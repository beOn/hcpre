"""
Routines for extracting data from Siemens DICOM files.

The simplest way to read a file is to call read(filename). If you like you
can also call lower level functions like read_data().

Except for the map of internal data types to numpy type strings (which
doesn't require an import of numpy), this code is deliberately ignorant of
numpy. It returns native Python types that are easy to convert into
numpy types.
"""

# Python modules
from __future__ import division
import struct
import exceptions
import math

# 3rd party modules
import dicom

# Our modules
import util_mrs_file
import constants

TYPE_NONE = 0
TYPE_IMAGE = 1
TYPE_SPECTROSCOPY = 2


# Change to True to enable the assert() statements sprinkled through the code
ASSERTIONS_ENABLED = False

# THese are some Siemens-specific tags
TAG_CONTENT_TYPE = (0x0029, 0x1008)
TAG_SPECTROSCOPY_DATA = (0x7fe1, 0x1010)

# I (Philip) ported much of the private tag parsing code from the IDL routines
# dicom_fill_rsp.pro and dicom_fill_util.pro, except for the CSA header
# parsing which is a port of C++ code in the GDCM project.

# Since a lot (all?) of the Siemens format is undocumented, there are magic
# numbers and logic in here that I can't explain. Sorry! Where appropriate
# I have copied or paraphrased comments from the IDL code; they're marked
# with [IDL]. Unmarked comments are mine. Where ambiguous, I labelled my
# comments with [PS] (Philip Semanchuk).


def read(filename, ignore_data=False):
    """ This is the simplest (and recommended) way for our code to read a
    Siemens DICOM file.
    
    It returns a tuple of (parameters, data). The parameters are a dict.
    The data is in a Python list.
    """
    # Since a DICOM file is params + data together, it's not so simple to
    # ignore the data part. The best we can do is tell PyDicom to apply
    # lazy evaluation which is probably less efficient in the long run.
    defer_size = 4096 if ignore_data else 0
    
    dataset = dicom.read_file(filename)
    
    params = read_parameters_from_dataset(dataset)
    data = read_data_from_dataset(dataset)
    
    return params, data
    
    
def read_parameters(filename):
    return read_parameters_from_dataset(dicom.read_file(filename))


def read_data(filename):
    return read_data_from_dataset(dicom.read_file(filename))

    
def read_data_from_dataset(dataset):
    """Given a PyDicom dataset, returns the data in the Siemens DICOM 
    spectroscopy data tag (0x7fe1, 0x1010) as a list of complex numbers.    
    """
    data = _get(dataset, TAG_SPECTROSCOPY_DATA)
    if data:
        # Big simplifying assumptions --
        # 1) Data is a series of complex numbers organized as ririri... 
        #    where r = real and i = imaginary.
        # 2) Each real & imaginary number is a 4 byte float.
        # 3) Data is little endian.
        data = struct.unpack("<%df" % (len(data) / 4), data)
        data = util_mrs_file.collapse_complexes(data)
    else:
        data = [ ]
        
    return data


def read_parameters_from_dataset(dataset):
    """Given a PyDicom dataset, returns a fairly extensive subset of the
    parameters therein as a dictionary.
    """
    params = { }
    
    # The code below refers to slice_index as a variable, but here it is 
    # hardcoded to one. It could vary, in theory, but in practice I don't
    # know how it would actually be used. How would the slice index or 
    # indices be passed? How would the data be returned? For now, I'll 
    # leave the slice code active but hardcode the index to 1.
    slice_index = 1
    
    # [PS] - Even after porting this code I still can't figure out what
    # ptag_img and ptag_ser stand for, so I left the names as is.
    ptag_img = { }
    ptag_ser = { }

    # (0x0029, 0x__10) is one of several possibilities
    # - SIEMENS CSA NON-IMAGE, CSA Data Info
    # - SIEMENS CSA HEADER, CSA Image Header Info
    # - SIEMENS CSA ENVELOPE, syngo Report Data
    # - SIEMENS MEDCOM HEADER, MedCom Header Info
    # - SIEMENS MEDCOM OOG, MedCom OOG Info (MEDCOM Object Oriented Graphics)
    # Pydicom identifies it as "CSA Image Header Info"
    for tag in ( (0x0029, 0x1010), (0x0029, 0x1210), (0x0029, 0x1110) ):
        tag_data = dataset.get(tag, None)
        if tag_data:
            break

    if tag_data:
        ptag_img = _parse_csa_header(tag_data.value)

    # [IDL] Access the SERIES Shadow Data
    # [PS] I don't know what makes this "shadow" data.
    for tag in ( (0x0029, 0x1020), (0x0029, 0x1220), (0x0029, 0x1120) ):
        tag_data = dataset.get(tag, None)
        if tag_data:
            break

    if tag_data:
        ptag_ser = _parse_csa_header(tag_data.value)

    # [IDL] "MrProtocol" (VA25) and "MrPhoenixProtocol" (VB13) are special
    # elements that contain many parameters.
    if ptag_ser.get("MrProtocol", ""):
        prot_ser = _parse_protocol_data(ptag_ser["MrProtocol"])

    if ptag_ser.get("MrPhoenixProtocol", ""):
        prot_ser = _parse_protocol_data(ptag_ser["MrPhoenixProtocol"])

    # [IDL] Determine if file is SVS,SI,EPSI, or OTHER
    # [PS] IDL code doesn't match comments. Possibilities appear to
    #      include EPSI, SVS, CSI, JPRESS and SVSLIP2. "OTHER" isn't
    #      considered.
    #      EPSI = Echo-Planar Spectroscopic Imaging
    #      SVS = Single voxel spectroscopy
    #      CSI = Chemical Shift Imaging
    #      JPRESS = J-resolved spectroscopy
    #      SVSLIP2 = No idea!
    is_epsi = False
    is_svs = False
    is_csi = False
    is_jpress = False
    is_svslip2 = False

    # [IDL] Protocol name
    parameter_filename = _extract_from_quotes(prot_ser.get("tProtocolName", ""))
    parameter_filename = parameter_filename.strip()

    # [IDL] Sequence file name
    sequence_filename = _extract_from_quotes(prot_ser.get("tSequenceFileName", ""))
    sequence_filename = sequence_filename.strip()

    sequence_filename2 = ptag_img.get("SequenceName", "")
    sequence_filename2 = sequence_filename2.strip()

    parameter_filename_lower = parameter_filename.lower()
    sequence_filename_lower = sequence_filename.lower()
    sequence_filename2_lower = sequence_filename2.lower()

    is_epsi = ("epsi" in (parameter_filename_lower, sequence_filename_lower))

    is_svs = ("svs" in (parameter_filename_lower, sequence_filename_lower,
                        sequence_filename2_lower))

    if "fid" in (parameter_filename_lower, sequence_filename_lower):
        if "csi" in (parameter_filename_lower, sequence_filename_lower):
            is_csi = True
        else:
            is_svs = True

    if "csi" in (parameter_filename_lower, sequence_filename_lower):
        is_csi = True

    is_jpress = ("jpress" in (parameter_filename_lower,
                              sequence_filename_lower))

    is_svslip2 = ("svs_li2" in (parameter_filename_lower,
                                sequence_filename2_lower))


    # Patient Info
    params["patient_name"] = _get(dataset, (0x0010, 0x0010), "")
    params["patient_id"] = _get(dataset, (0x0010, 0x0020))
    params["patient_birthdate"] = _get(dataset, (0x0010, 0x0030))
    params["patient_sex"] = _get(dataset, (0x0010, 0x0040), "")
    # [PS] Siemens stores the age as nnnY where 'n' is a digit, e.g. 042Y
    params["patient_age"] = \
                        int(_get(dataset, (0x0010, 0x1010), "000Y")[:3])
    params["patient_weight"] = round(_get(dataset, (0x0010, 0x1030), 0))
    params["study_code"] = _get(dataset, (0x0008, 0x1030), "")

    # Identification info
    params["bed_move_fraction"] = 0.0

    s = _get(dataset, (0x0008, 0x0080), "")
    if s:
        s = " " + s
    s += _get(dataset, (0x0008, 0x1090), "")
    params["institution_id"] = s
    params["parameter_filename"] = parameter_filename
    params["study_type"] = "spec"

    # DICOM date format is YYYYMMDD
    params["bed_move_date"] = _get(dataset, (0x0008, 0x0020), "")
    params["measure_date"] = params["bed_move_date"]
    # DICOM time format is hhmmss.fraction
    params["bed_move_time"] = _get(dataset, (0x0008, 0x0030), "")
    params["comment_1"] = _get(dataset, (0x0008, 0x0031), "")
    if not params["comment_1"]:
        params["comment_1"] = _get(dataset, (0x0020, 0x4000), "")
    # DICOM time format is hhmmss.fraction
    params["measure_time"] = _get(dataset, (0x0008, 0x0032), "")

    params["sequence_filename"] = ptag_img.get("SequenceName", "")
    params["sequence_type"] = ptag_img.get("SequenceName", "")

    # Measurement info
    params["echo_position"] = "0.0"
    params["image_contrast_mode"] = "unknown"
    params["kspace_mode"] = "unknown"
    params["measured_slices"] = "1"
    params["saturation_bands"] = "0"

    # Seems to me that a quantity called "NumberOfAverages" would be an
    # int, but it is stored as a float, e.g. "128.0000" which makes
    # Python's int() choke unless I run it through float() first.
    params["averages"] = int(_float(ptag_img.get("NumberOfAverages", "")))
    params["flip_angle"] = _float(ptag_img.get("FlipAngle", ""))
    # [PS] DICOM stores frequency as MHz, we store it as Hz. Mega = 1x10(6)
    params["frequency"] = float(ptag_img.get("ImagingFrequency", 0)) * 1e6

    inversion_time = float(ptag_img.get("InversionTime", 0))
    params["inversion_time_1"] = inversion_time
    params["number_inversions"] = 1 if inversion_time else 0

    params["measured_echoes"] = ptag_img.get("EchoTrainLength", "1")
    params["nucleus"] = ptag_img.get("ImagedNucleus", "")
    params["prescans"] = prot_ser.get("sSpecPara.lPreparingScans", 0)

    # Gain
    gain = prot_ser.get("sRXSPEC.lGain", None)
    if gain == 0:
        gain = "-20.0"
    elif gain == 1:
        gain = "0.0"
    else:
        gain = ""

    params["receiver_gain"] = gain

    params["ft_scale_factor"] = \
                float(prot_ser.get("sRXSPEC.aFFT_SCALE[0].flFactor", 0))

    # Receiver Coil
    coil = prot_ser.get("sCOIL_SELECT_MEAS.asList[0].sCoilElementID.tCoilID", "")
    params["receiver_coil"] = _extract_from_quotes(coil)

    # [IDL] differs in EPSI
    params["repetition_time_1"] = float(prot_ser.get("alTR[0]", 0)) * 0.001

    sweep_width = ""
    remove_oversample_flag = prot_ser.get("sSpecPara.ucRemoveOversampling", "")
    remove_oversample_flag = (remove_oversample_flag.strip() == "0x1")
    readout_os = float(ptag_ser.get("ReadoutOS", 1.0))

    dwelltime = float(ptag_img.get("RealDwellTime", 1.0)) * 1e-9

    if dwelltime:
        sweep_width = 1 / dwelltime
        if not remove_oversample_flag:
            sweep_width *= readout_os

        sweep_width = str(sweep_width)

    params["transmitter_voltage"] = \
        prot_ser.get("sTXSPEC.asNucleusInfo[0].flReferenceAmplitude", "0.0")
    params["total_duration"] = \
        prot_ser.get("lTotalScanTimeSec", "0.0")

    prefix = "sSliceArray.asSlice[%d]." % slice_index
    image_parameters = (
                        ("image_dimension_line",      "dPhaseFOV"),
                        ("image_dimension_column",    "dReadoutFOV"),
                        ("image_dimension_partition", "dThickness"),
                        ("image_position_sagittal",   "sPosition.dSag"),
                        ("image_position_coronal",    "sPosition.dCor"),
                        ("image_position_transverse", "sPosition.dTra"),
                       )
    for key, name in image_parameters:
        params[key] = float(prot_ser.get(prefix + name, "0.0"))

    # [IDL] Image Normal/Column
    image_orientation = ptag_img.get("ImageOrientationPatient", "")

    if not image_orientation:
        slice_orientation_pitch = ""
        slice_distance = ""
    else:
        # image_orientation is a list of strings, e.g. --
        # ['-1.00000000', '0.00000000', '0.00000000', '0.00000000', 
        #   '1.00000000', '0.00000000']

        # [IDL] If the data we are processing is a Single Voxel
        # Spectroscopy data, interchange rows and columns. Due to an error
        # in the protocol used.
        if is_svs:
            image_orientation = image_orientation[3:] + image_orientation[:3]

        # Convert the values to float and discard ones smaller than 1e-4
        f = lambda value: 0.0 if abs(value) < 1e-4 else value
        image_orientation = [f(float(value)) for value in image_orientation]

        row = image_orientation[:3]
        column = image_orientation[3:6]
        normal = ( ((row[1] * column[2]) - (row[2] * column[1])), 
                   ((row[2] * column[0]) - (row[0] * column[2])),
                   ((row[0] * column[1]) - (row[1] * column[0])),
                 )

        params["image_normal_sagittal"] = normal[0]
        params["image_normal_coronal"] = normal[1]
        params["image_normal_transverse"] = normal[2]
        params["image_column_sagittal"] = column[0]
        params["image_column_coronal"] = column[0]
        params["image_column_transverse"] = column[0]

        # Second part of the return tuple is orientation; we don't use it.
        slice_orientation_pitch, _ = _dicom_orientation_string(normal)

        # Slice distance
        # http://en.wikipedia.org/wiki/Dot_product
        keys = ("image_position_sagittal", "image_position_coronal",
                "image_position_transverse")
        a = [params[key] for key in keys]
        b = normal
        bb = math.sqrt(sum([value ** 2 for value in normal]))

        slice_distance = ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])) / bb


    params["slice_orientation_pitch"] = slice_orientation_pitch
    params["slice_distance"] = slice_distance

    regions = ( ("region_dimension_line",      "dPhaseFOV"),
                ("region_dimension_column",    "dReadoutFOV"),
                ("region_dimension_partition", "dThickness"),
                ("region_position_sagittal",   "sPosition.dSag"),
                ("region_position_coronal",    "sPosition.dCor"),
                ("region_position_transverse", "sPosition.dTra"),
              )

    for key, name in regions:
        name = "sSpecPara.sVoI." + name
        params[key] = float(prot_ser.get(name, 0))

    # 'DATA INFORMATION'

    params["measure_size_spectral"] = \
                                long(prot_ser.get('sSpecPara.lVectorSize', 0))

    params["slice_thickness"] = _float(ptag_img.get("SliceThickness", 0))
    params["current_slice"] = "1"
    params["number_echoes"] = "1"
    params["number_slices"] = "1"
    params["data_size_spectral"] = params["measure_size_spectral"]

    # ;------------------------------------------------------
    # [IDL] Sequence Specific Changes
    if not is_epsi:
        # [IDL] Echo time - JPRESS handling added by Dragan
        echo_time = 0.0

        if is_jpress:
            # [IDL] Yingjian saves echo time in a private 'echotime' field
            # [PS] The IDL code didn't use a dict to store these values
            # but instead did a brute force case-insensitive search over
            # an array of strings. In that context, key case didn't matter
            # but here it does.
            keys = prot_ser.keys()
            for key in keys:
                if key.upper() == "ECHOTIME":
                    echo_time = float(prot_ser[key])

        if is_svslip2:
            # [IDL] BJS found TE value set in ICE to be updated in
            #      'echotime' field
            # [PS] The IDL code didn't use a dict to store these values
            # but instead did a brute force case-insensitive search over
            # an array of strings. In that context, key case didn't matter
            # but here it does.
            keys = ptag_img.keys()
            for key in keys:
                if key.upper() == "ECHOTIME":
                    echo_time = float(ptag_img[key])

        if not echo_time:
            # [IDL] still no echo time - try std place
            echo_time = float(prot_ser.get('alTE[0]', 0.0))
            echo_time /= 1000

        params["echo_time"] = echo_time

        params["data_size_line"] = \
                    int(prot_ser.get('sSpecPara.lFinalMatrixSizePhase', 1))

        params["data_size_column"] = \
                        int(prot_ser.get('sSpecPara.lFinalMatrixSizeRead', 1))

        params["data_size_partition"] = \
                    int(prot_ser.get('sSpecPara.lFinalMatrixSizeSlice', 1))

        if is_svs:
            # [IDL] For Single Voxel Spectroscopy data (SVS) only
            params["image_dimension_line"] = \
                                        params["region_dimension_line"]
            params["image_dimension_column"] = \
                                        params["region_dimension_column"]
            params["image_dimension_partition"] = \
                                        params["region_dimension_partition"]


            # [IDL] For SVS data the following three parameters cannot be
            # anything other than 1
            params["measure_size_line"] = 1
            params["measure_size_column"] = 1
            params["measure_size_partition"] = 1
        else:
            # Not SVS
            #   ;--------------------------------------------------
            #   ; [IDL] For CSI or OTHER Spectroscopy data only
            #   ;--------------------------------------------------

            measure_size_line = int(prot_ser.get('sKSpace.lPhaseEncodingLines', 1))
            params["measure_size_line"] = str(measure_size_line)

            measure_size_column = int(prot_ser.get('sKSpace.lPhaseEncodingLines', 0))
            params["measure_size_column"] = str(measure_size_column)

            measure_size_partition = int(prot_ser.get('sKSpace.lPartitions', '0'))
            kspace_dimension = prot_ser.get('sKSpace.ucDimension', '')
            if kspace_dimension.strip() == "0x2":
                measure_size_partition = 1
                params["data_size_partition"] = 1
                data_size_partition = 1

            params["measure_size_partition"] = measure_size_partition


    if sequence_filename in ("svs_cp_press", "svs_se_ir", "svs_tavg"):
        # [IDL] Inversion Type 0-Volume,1-None
        s = prot_ser.get("SPREPPULSES.UCINVERSION", "")
        if s == "0x1":
            params["number_inversions"] = 1
        elif s == "0x2":
            params["number_inversions"] = 0
        # else:
        #     params["number_inversions"] doesn't get set at all.
        #     This matches the behavior of the IDL code. Note that
        #     params["number_inversions"] is also populated
        #     unconditionally in code many lines above.


    if sequence_filename in ("svs_se", "svs_st", "fid", "fid3", "fid_var",
                             "csi_se", "csi_st", "csi_fid", "csi_fidvar",
                             "epsi"):
         # [IDL] FOR EPSI Measure_size and Data_size parameters are the same
        params["region_dimension_line"] = \
                                params["image_dimension_line"]
        params["region_dimension_column"] = \
                                params["image_dimension_column"]
        params["ft_scale_factor"] = "1.0"

        params["data_size_line"] = \
                int(prot_ser.get('sKSpace.lPhaseEncodingLines', 0))
        params["data_size_column"] = \
                int(prot_ser.get('sKSpace.lBaseResolution', 0)) * readout_os
        params["data_size_partition"] = \
                int(prot_ser.get('sKSpace.lPartitions', 0))

        params["measure_size_line"] = params["data_size_line"]
        measure_size_column = params["data_size_column"]
        measure_size_partition = params["data_size_partition"]

        index = 0 if ((int(dataset.get("InstanceNumber", 0)) % 2) == 1) else 1
        echo_time =         float(prot_ser.get('alTE[%d]' % index, 0)) / 1000
        repetition_time_1 = float(prot_ser.get('alTR[%d]' % index, 0)) / 1000

        params["echo_time"] = str(echo_time)
        params["repetition_time_1"] = str(repetition_time_1)

        dwelltime = float(ptag_img.get("RealDwellTime", 0.0))
        if dwelltime and base_resolution:
            sweep_width = 1 / (dwelltime * base_resolution * readout_os)
        else:
            sweep_width = ""

        params["sweep_width"] = sweep_width

    # Added by BTA
    ip_rot = prot_ser.get("sSliceArray.asSlice[0].dInPlaneRot", None)
    pol_swap = prot_ser.get("sWipMemBlock.alFree[40]", None)
    if ip_rot:
        try:
            ip_rot = float(ip_rot)
            params["in_plane_rotation"] = ip_rot
        except Exception, e:
            pass
    if pol_swap:
        try:
            pol_swap = int(pol_swap)
            params["polarity_swap"] = pol_swap
        except Exception, e:
            raise e

    return params


def _my_assert(expression):
    if ASSERTIONS_ENABLED:
        assert(expression)


def _dicom_orientation_string(normal):
    """Given a 3-item list (or other iterable) that represents a normal vector
    to the "imaging" plane, this function determines the orientation of the
    vector in 3-dimensional space. It returns a tuple of (angle, orientation)
    in which angle is e.g. "Tra" or "Tra>Cor -6" or "Tra>Sag 14.1 >Cor 9.3"
    and orientation is e.g. "Sag" or "Cor-Tra".

    For double angulation, errors in secondary angle occur that may be due to
    rounding errors in internal Siemens software, which calculates row and
    column vectors.
    """
    # docstring paraphrases IDL comments
    TOLERANCE = 1.e-4
    orientations = ('Sag', 'Cor', 'Tra')

    final_angle = ""
    final_orientation = ""

    # [IDL] evaluate orientation of normal vector:
    #
    # Find principal direction of normal vector (i.e. axis with its largest
    # component)
    # Find secondary direction (second largest component)
    # Calc. angle btw. projection of normal vector into the plane that
    #     includes both principal and secondary directions on the one hand
    #     and the principal direction on the other hand ==> 1st angulation:
    #     "principal>secondary = angle"
    # Calc. angle btw. projection into plane perpendicular to principal
    #     direction on the one hand and secondary direction on the other
    #     hand ==> 2nd angulation: "secondary>third dir. = angle"


    # get principal, secondary and ternary directions
    sorted_normal = sorted(normal)

    for i, value in enumerate(normal):
        if value == sorted_normal[2]:
            # [IDL] index of principal direction
            principal = i
        if value == sorted_normal[1]:
            # [IDL] index of secondary direction
            secondary = i
        if value == sorted_normal[0]:
            # [IDL] index of ternary direction
            ternary = i

    # [IDL] calc. angle between projection into third plane (spawned by
    # principle & secondary directions) and principal direction:
    angle_1 = math.atan2(normal[secondary], normal[principal]) * \
                                                constants.RADIANS_TO_DEGREES

    # [IDL] calc. angle btw. projection on rotated principle direction and
    # secondary direction:
    # projection on rotated principle dir.
    new_normal_ip = math.sqrt((normal[principal] ** 2) + (normal[secondary] ** 2))

    angle_2 = math.atan2(normal[ternary], new_normal_ip) * \
                                                constants.RADIANS_TO_DEGREES

    # [IDL] SIEMENS notation requires modifications IF principal dir. indxs SAG !
    # [PS] In IDL, indxs is the name of the variable that is "secondary" here.
    #      Even with that substitution, I don't understand the comment above.
    if not principal:
        if abs(angle_1) > 0:
            sign1 = angle_1 / abs(angle_1)
        else:
            sign1 = 1.0

        angle_1 -= (sign1 * 180.0)
        angle_2 *= -1

    if (abs(angle_2) < TOLERANCE) or (abs(abs(angle_2) - 180) < TOLERANCE):
        if (abs(angle_1) < TOLERANCE) or (abs(abs(angle_1) - 180) < TOLERANCE):
            # [IDL] NON-OBLIQUE:
            final_angle = orientations[principal]
            final_orientation = ang
        else:
            # [IDL] SINGLE-OBLIQUE:
            final_angle = "%s>%s %.3f" % \
                    (orientations[principal], orientations[secondary],
                     (-1 * angle_1)
                    )
            final_orientation = orientations[principal] + '-' + orientations[secondary]
    else:
        # [IDL] DOUBLE-OBLIQUE:
        final_angle = "%s>%s %.3f >%s %f" % \
                (orientations[principal], orientations[secondary],
                 (-1 * angle_1), orientations[ternary], (-1 * angle_2))
        final_orientation = "%s-%s-%s" % \
                (orientations[principal], orientations[secondary],
                 orientations[ternary])

    return final_angle, final_orientation


def _float(value):
    """Attempts to return value as a float. No different from Python's 
    built-in float(), except that it accepts None and "" (for which it 
    returns 0.0).
    """
    return float(value) if value else 0.0
    
    
def _extract_from_quotes(s):
    """Given a string, returns the portion between the first and last
    double quote (ASCII 34). If there aren't at least two quote characters,
    the original string is returned."""
    start = s.find('"')
    end = s.rfind('"')

    if (start != -1) and (end != -1):
        s = s[start + 1 : end]

    return s


def _null_truncate(s):
    """Given a string, returns a version truncated at the first '\0' if
    there is one. If not, the original string is returned."""
    i = s.find(chr(0))
    if i != -1:
        s = s[:i]

    return s


def _scrub(item):
    """Given a string, returns a version truncated at the first '\0' and
    stripped of leading/trailing whitespace. If the param is not a string,
    it is returned unchanged."""
    if isinstance(item, basestring):
        return _null_truncate(item).strip()
    else:
        return item


def _get_chunks(tag, index, format, little_endian=True):
    """Given a CSA tag string, an index into that string, and a format
    specifier compatible with Python's struct module, returns a tuple
    of (size, chunks) where size is the number of bytes read and
    chunks are the data items returned by struct.unpack(). Strings in the
    list of chunks have been run through _scrub().
    """
    # The first character of the format string indicates endianness.
    format = ('<' if little_endian else '>') + format
    size = struct.calcsize(format)
    chunks = struct.unpack(format, tag[index:index + size])

    chunks = [_scrub(item) for item in chunks]

    return (size, chunks)


def _parse_protocol_data(protocol_data):
    """Returns a dictionary containing the name/value pairs inside the
    "ASCCONV" section of the MrProtocol or MrPhoenixProtocol elements
    of a Siemens CSA Header tag.
    """
    # Protocol_data is a large string (e.g. 32k) that lists a lot of
    # variables in a JSONish format with which I'm not familiar. Following
    # that there's another chunk of data delimited by the strings you see
    # below.
    # That chunk is a list of name=value pairs, INI file style. We
    # ignore everything outside of the ASCCONV delimiters. Everything inside
    # we parse and return as a dictionary.
    start = protocol_data.find("### ASCCONV BEGIN ###")
    end = protocol_data.find("### ASCCONV END ###")

    _my_assert(start != -1)
    _my_assert(end != -1)

    start += len("### ASCCONV BEGIN ###")
    protocol_data = protocol_data[start:end]

    lines = protocol_data.split('\n')

    # The two lines of code below turn the 'lines' list into a list of
    # (name, value) tuples in which name & value have been stripped and
    # all blank lines have been discarded.
    f = lambda pair: (pair[0].strip(), pair[1].strip())
    lines = [f(line.split('=')) for line in lines if '=' in line]

    return dict(lines)


def _get(dataset, tag, default=None):
    """Returns the value of a dataset tag, or the default if the tag isn't
    in the dataset.
    PyDicom datasets already have a .get() method, but it returns a
    dicom.DataElement object. In practice it's awkward to call dataset.get()
    and then figure out if the result is the default or a DataElement,
    and if it is the latter _get the .value attribute. This function allows
    me to avoid all that mess.
    It is also a workaround for this bug (which I submitted) which should be
    fixed in PyDicom > 0.9.3:
    http://code.google.com/p/pydicom/issues/detail?id=72
    """
    return default if tag not in dataset else dataset[tag].value


def _parse_csa_header(tag, little_endian = True):
    """The CSA header is a Siemens private tag that should be passed as
    a string. Any of the following tags should work: (0x0029, 0x1010),
    (0x0029, 0x1210), (0x0029, 0x1110), (0x0029, 0x1020), (0x0029, 0x1220),
    (0x0029, 0x1120).

    The function returns a dictionary keyed by element name.
    """
    # Let's have a bit of fun, shall we? A Siemens CSA header is a mix of
    # binary glop, ASCII, binary masquerading as ASCII, and noise masquerading
    # as signal. It's also undocumented, so there's no specification to which
    # to refer.

    # The format is a good one to show to anyone who complains about XML being
    # verbose or hard to read. Spend an afternoon with this and XML will
    # look terse and read like a Shakespearean sonnet.

    # The algorithm below is a translation of the GDCM project's
    # CSAHeader::LoadFromDataElement() inside gdcmCSAHeader.cxx. I don't know
    # how that code's author figured out what's in a CSA header, but the
    # code works.

    # I added comments and observations, but they're inferences. I might
    # be wrong. YMMV.

    # Some observations --
    # - If you need to debug this code, a hexdump of the tag data will be
    #   your best friend.
    # - The data in the tag is a list of elements, each of which contains
    #   zero or more subelements. The subelements can't be further divided
    #   and are either empty or contain a string.
    # - Everything begins on four byte boundaries.
    # - This code will break on big endian data. I don't know if this data
    #   can be big endian, and if that's possible I don't know what flag to
    #   read to indicate that. However, it's easy to pass an endianness flag
    #   to _get_chunks() should the need to parse big endian data arise.
    # - Delimiters are thrown in here and there; they are 0x4d = 77 which is
    #   ASCII 'M' and 0xcd = 205 which has no ASCII representation.
    # - Strings in the data are C-style NULL terminated.

    # I sometimes read delimiters as strings and sometimes as longs.
    DELIMITERS = ("M", "\xcd", 0x4d, 0xcd)

    # This dictionary of elements is what this function returns
    elements = { }

    # I march through the tag data byte by byte (actually a minimum of four
    # bytes at a time), and current points to my current position in the tag
    # data.
    current = 0

    # The data starts with "SV10" followed by 0x04, 0x03, 0x02, 0x01.
    # It's meaningless to me, so after reading it, I discard it.
    size, chunks = _get_chunks(tag, current, "4s4s")
    current += size

    _my_assert(chunks[0] == "SV10")
    _my_assert(chunks[1] == "\4\3\2\1")

    # get the number of elements in the outer list
    size, chunks = _get_chunks(tag, current, "L")
    current += size
    element_count = chunks[0]

    # Eat a delimiter (should be 0x77)
    size, chunks = _get_chunks(tag, current, "4s")
    current += size
    _my_assert(chunks[0] in DELIMITERS)

    for i in range(element_count):
        # Each element looks like this:
        # - (64 bytes) Element name, e.g. ImagedNucleus, NumberOfFrames,
        #   VariableFlipAngleFlag, MrProtocol, etc. Only the data up to the
        #   first 0x00 is important. The rest is helpfully populated with
        #   noise that has enough pattern to make it look like something
        #   other than the garbage that it is.
        # - (4 bytes) VM
        # - (4 bytes) VR
        # - (4 bytes) syngo_dt
        # - (4 bytes) # of subelements in this element (often zero)
        # - (4 bytes) a delimiter (0x4d or 0xcd)
        size, chunks = _get_chunks(tag, current,
                                  "64s" + "4s" + "4s" + "4s" + "L" + "4s")
        current += size

        name, vm, vr, syngo_dt, subelement_count, delimiter = chunks
        _my_assert(delimiter in DELIMITERS)

        # The subelements hold zero or more strings. Those strings are stored
        # temporarily in the values list.
        values = [ ]

        for j in range(subelement_count):
            # Each subelement looks like this:
            # - (4 x 4 = 16 bytes) Call these four bytes A, B, C and D. For
            #   some strange reason, C is always a delimiter, while A, B and
            #   D are always equal to one another. They represent the length
            #   of the associated data string.
            # - (n bytes) String data, the length of which is defined by
            #   A (and A == B == D).
            # - (m bytes) Padding if length is not an even multiple of four.
            size, chunks = _get_chunks(tag, current, "4L")
            current += size

            _my_assert(chunks[0] == chunks[1])
            _my_assert(chunks[1] == chunks[3])
            _my_assert(chunks[2] in DELIMITERS)
            length = chunks[0]

            # get a chunk-o-stuff, length indicated by code above.
            # Note that length can be 0.
            size, chunks = _get_chunks(tag, current, "%ds" % length)
            current += size
            if chunks[0]:
                values.append(chunks[0])

            # If we're not at a 4 byte boundary, move.
            # Clever modulus code below swiped from GDCM
            current += (4 - (length % 4)) % 4

        # The value becomes a single string item (possibly "") or a list
        # of strings
        if len(values) == 0:
            values = ""
        if len(values) == 1:
            values = values[0]

        _my_assert(name not in elements)
        elements[name] = values

    return elements
    
