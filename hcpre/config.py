import os
import sys
from multiprocessing import Pool, cpu_count
from configobj import ConfigObj

from hcpre.util import *

SCAN_TYPES = [
    "bold",
    "bold_sbref",
    "fieldmap_magnitude",
    "fieldmap_phase",
    "fieldmap_ap",
    "fieldmap_lr",
    "fieldmap_pa",
    "fieldmap_rl",
    "t1",
    "t2",
    "polarity_swapped"]

TEMPL_KEYS = [
    "templates_dir",
    "t1_template",
    "t1_template_brain",
    "t1_template_2mm",
    "t2_template",
    "t2_template_brain",
    "t2_template_2mm",
    "template_mask",
    "template_2mm_mask",
    "surf_atlas_dir",
    "grayordinates_dir",
    "ref_myelin_maps",
    ]

CONF_FILE_KEYS = [
    "fnirt_config",
    "grad_distort_coeffs",
    "subcort_gray_labels",
    "freesurfer_labels",
    "top_up_config",]

CONF_FILE_DEFAULTS = [
    "%(hcp_dir)s/global/config/T1_2_MNI152_2mm.cnf",
    "%(hcp_dir)s/global/config/coeff_SC72C_Skyra.grad",
    "%(hcp_dir)s/global/config/FreeSurferSubcorticalLabelTableLut.txt",
    "%(hcp_dir)s/global/config/FreeSurferAllLut.txt",
    "%(hcp_dir)s/global/config/b02b0.cnf",]

POS_FIELDMAPS = ["fieldmap_rl", "fieldmap_pa"]
NEG_FIELDMAPS = ["fieldmap_lr", "fieldmap_ap"]
YES_WORDS = [1,"1","y","Y","yes","Yes","YES", "True", "true"]
NO_WORDS = [0,"0","n","N","no","No","NO", "False", "false"]
EPFM_SELECTION_OPTS = ["first","most_recent"]

def setup_conf():
    import os
    # get name for conf file
    name = raw_input("New config file name [hcp.conf]: ")
    name = name.strip()
    name = name if name else "hcp.conf"
    name = name if (len(name) > 5 and name[-5:] == ".conf") else name + ".conf"
    # make sure that file doesn't already exist
    if os.path.exists(name):
        print "File already exists."
        return
    # get the directory containing all data for all subjects
    print "\nThe subjects directory should contain all raw data for all subjects."
    subs_dir = raw_input("Subjects Directory [./]: ")
    subs_dir = subs_dir.strip()
    subs_dir = subs_dir if subs_dir else "./"
    subs_dir = os.path.abspath(subs_dir)
    # get the template used for getting form the subject dir to the dicoms
    print "\nThe DICOM template should be a format string for a glob which, " + \
          "when combined with an individual subject ID, will get us all of " + \
          "the subject's DICOM files."
    dcm_temp = raw_input("DICOM template [data/raw_dicom/%s/*.dcm]: ")
    dcm_temp = dcm_temp.strip()
    dcm_temp = dcm_temp if dcm_temp else "data/raw_dicom/%s/*.dcm"
    dcm_temp = dcm_temp if ".dcm" in dcm_temp else os.path.join(dcm_temp, "*.dcm")
    # get a list of subjects
    print "\nSubjects should be a comma separated list of subject ids."
    subs = raw_input("Subject list ['']: ")
    subs = subs.strip()
    # set up config obj
    config = ConfigObj(name, unrepr=True)
    # get the basics
    config["general"] = {}
    config["general"]["subjects"] = [s.strip() for s in subs.split(",")]
    config["general"]["subject_dir"] = subs_dir
    config["general"]["dicom_template"] = dcm_temp
    # write the config file
    config.write()
    update_conf(name)

def get_series_desc(dicom_path):
    import dicom
    d = dicom.read_file(dicom_path, stop_before_pixels=True)
    sd = getattr(d, "SeriesDescription", None)
    return sd

def templ_defaults_for_rez(rez):
    if rez in [.7, .8]:
        rez = "%0.1f" % rez
    elif rez in [1]:
        rez = "%1.f" % rez
    else:
        return None
    return [
        "%(hcp_dir)s/global/templates",
        "%%(templates_dir)s/MNI152_T1_%smm.nii.gz" % rez,
        "%%(templates_dir)s/MNI152_T1_%smm_brain.nii.gz" % rez,
        "%(templates_dir)s/MNI152_T1_2mm.nii.gz",
        "%%(templates_dir)s/MNI152_T2_%smm.nii.gz" % rez,
        "%%(templates_dir)s/MNI152_T2_%smm_brain.nii.gz" % rez,
        "%(templates_dir)s/MNI152_T2_2mm.nii.gz",
        "%%(templates_dir)s/MNI152_T1_%smm_brain_mask.nii.gz" % rez,
        "%(templates_dir)s/MNI152_T1_2mm_brain_mask_dil.nii.gz",
        "%(templates_dir)s/standard_mesh_atlases",
        "%(templates_dir)s/91282_Greyordinates",
        "%(templates_dir)s/standard_mesh_atlases/Conte69.MyelinMap_BC.164k_fs_LR.dscalar.nii"]

def update_conf(conf_path):
    import os
    import sys
    from glob import glob
    import dicom
    from time import sleep
    from numpy import unique
    # TODO: make sure conf is there
    config = ConfigObj(conf_path, unrepr=True)
    # update the series map if needed or desired
    needs_smap = "series" not in config
    if not needs_smap:
        needs_smap = raw_input("Do you want to re-map the series descriptions y/[n]? ").strip() in YES_WORDS
    if needs_smap:
        # get list of all sequence descriptions
        s = config["general"]["subject_dir"]
        t = config["general"]["dicom_template"] % "*"
        dicoms = glob(os.path.join(s,t))
        message = "\rChecking series names (this may several minutes) %s"
        dcm_count = len(dicoms)
        pool = Pool(processes=min(15, int(round(cpu_count() * .75))))
        result = pool.map_async(get_series_desc, dicoms)
        while not result.ready():
            sleep(.5)
            # perc = float(dcm_count - result._number_left) / float(dcm_count)
            # perc = "(%d%%)..." % int(round(perc * 100))
            perc = "(%d chunks remaining)..." % result._number_left
            m = message % perc
            sys.stdout.write(m)
            sys.stdout.flush()
        series = unique(result.get()).tolist()
        if None in series:
            series.remove(None)
        found = "\nFound %d unique series descriptions.\n" % len(series)
        sys.stdout.write(found)
        # message if we didn't find anything
        if not series:
            raise ValueError("couldn't find any dicoms! please double check your paths and templates...")
        # update the series confg
        series.sort()
        print "-------\nSeries:\n-------"
        print "\n".join(["%d:\t%s" % (i, ser) for i,ser in enumerate(series)]) + "\n"
        type_matches = dict(zip(SCAN_TYPES, [[] for s in SCAN_TYPES]))
        i = 0
        series_count = len(series)
        while i < len(SCAN_TYPES):
            ft = SCAN_TYPES[i]
            m = raw_input("\nWhich series do you use for '%s'?\n[None] or comma separated values 0-%d: " % (ft, len(series)-1)).strip()
            if not m:
                i += 1
                continue
            try:
                m = [int(n) for n in m.split(",")]
                if any([n < 0 or n >= series_count for n in m]):
                    raise ValueError()
            except Exception, e:
                print "Invalid Selection..."
                continue
            i += 1
            for n in m:
                type_matches[ft].append(series[n])
        config["series"] = {}
        for key in sorted(type_matches.keys()):
            config["series"][key] = type_matches[key]
    config["DEFAULT"] = {}
    # freesurfer home (defaults)
    d_val = os.environ.get("FREESURFER_HOME", "")
    fsh = raw_input("\nPath for FREESURFER_HOME [%s]: " % d_val).strip()
    fsh = fsh if fsh else d_val
    config["DEFAULT"]["freesurfer_home"] = fsh
    # fsl dir (defaults)
    d_val = os.environ.get("FSLDIR", "")
    fslh = raw_input("\nPath for FSLDIR [%s]: " % d_val).strip()
    fslh = fslh if fslh else d_val
    config["DEFAULT"]["fsl_dir"] = fslh
    # hcp dir (defaults)
    hcph = raw_input("\nPath to HCP Pipelines dir [ ]: ").strip()
    config["DEFAULT"]["hcp_dir"] = hcph
    # template file locations (templates)
    rez = raw_input("\nWhat is your structural image resolution (mm)?\n[Skip] or one of (.7, .8, 1): ").strip()
    rez = float_or_none(rez)
    rez_temps = None
    def_t = False
    if rez:
        rez_temps = templ_defaults_for_rez(rez)
    if rez_temps:
        def_t = raw_input("\nUse default template files for resolution %0.1f [y]/n? " % rez).strip() not in NO_WORDS
        if def_t:
            t_vals = rez_temps
    if not t_vals:
        t_vals = ['' for k in TEMPL_KEYS]
    config["templates"] = {}
    for i, key in enumerate(TEMPL_KEYS):
        config["templates"][key] = t_vals[i]
    if not def_t:
        print "\nWhen finished, please open your config file to set the tamplate locations manually.\n"
    # config file locations (config_files)
    def_c = raw_input("\nUse default config files [y]/n?").strip() not in NO_WORDS
    c_vals = CONF_FILE_DEFAULTS if def_c else ['' for k in CONF_FILE_KEYS]
    config["config_files"] = {}
    for i, key in enumerate(CONF_FILE_KEYS):
        config["config_files"][key] = c_vals[i]
    if not def_c:
        print "\nWhen finished, please open your config file to set the other config file locations manually.\n"
    config["nifti_wrangler"] = {}
    # ep unwarp policy
    print "\nIf you have more than one ep fieldmap set, you may either\nwant to use the first, or always use the most recent."
    ep_sel = numbered_choice(EPFM_SELECTION_OPTS, prompt="Which policy would you like to use:")
    config["nifti_wrangler"]["ep_fieldmap_selection"] = ep_sel
    # struct averaging block
    print "\nIf you collect multiple t1 or t2 images, and averaging them yields warped\nresults, try blocking structural image averaging."
    block_avg = raw_input("\nBlock averaging of structural images [y]/n? ").strip() not in NO_WORDS
    config["nifti_wrangler"]["block_struct_averaging"] = block_avg
    # a terrible hack, but take a guess at the unwarp direction
    uwd = "x"
    ser = config["series"]
    if ser.get("fieldmap_ap",None) and ser.get("fieldmap_pa",None):
        uwd = "y"
    elif ser.get("fieldmap_rl",None) and ser.get("fieldmap_lr",None):
        uwd = "x"
    config["nifti_wrangler"]["ep_unwarp_dir"] = uwd
    print "\nVery weak guess that your primary unwarp direction is %s.\nDid I mention this is a GUESS?" % uwd
    print "\nWhen finished, please open your config file check the value for ep_unwarp_dir.\n"
    # output level selection
    print """
    By default, we only output the MNINonLinear results folder (recommended, for
    the sake of your drive capacity). If you want to check out the full outputs,
    including all of the precursor files, set output_mni_only to False in the
    config file."""
    config["output_select"] = {"output_mni_only":True}
    # write the file!
    config.write()

def select_conf():
    # if there"s only one conf around, select it. otherwise, offer a choice.
    from glob import glob
    confs = glob("./*.conf")
    if not confs:
        raise ValueError("Could not find any .conf files in current directory.")
    if len(confs) == 1:
        return confs[0]
    return numbered_choice(confs, prompt = "There were multiple config files. Please choose.")

def validate_config(conf_dict):
    """ validates the config dict
    return: True if it passes, False if it fails
    """
    # TODO: more of this
    # make sure we can get an env out of it
    try:
        get_hcp_env_for_config(conf_dict)
        get_hcp_commands_for_config(conf_dict)
    except Exception, e:
        return False
    return True

def numbered_choice(choices, prompt=None, allow_multiple=False, allow_none=False):
    import numpy as np
    # TODO: print a choice, return the chosen value
    if not prompt:
        prompt = "Options:\n" + "-"*len("Options:")
    print prompt + "\n"
    for i, val in enumerate(choices):
        print "%d: %s" % (i, val)
    # print "\n"
    nstr = "[None] or " if allow_none else ""
    mstr = "comma separated values " if allow_multiple else ""
    selection = None
    done = False
    while not done:
        try:
            m = raw_input("\nSelect %s%s0-%d: " % (nstr, mstr, len(choices)-1)).strip()
            if not m:
                if allow_none:
                    selection = None
                    done = True
                else:
                    raise ValueError()
            else:
                m = [int(n) for n in m.split(",")]
                if len(m) > 1 and not allow_multiple:
                    raise ValueError()
                if any([n < 0 or n >= len(choices) for n in m]):
                    raise ValueError()
                m = np.unique(m)
                selection = [choices[n] for n in m]
                if not allow_multiple:
                    selection = selection[0]
                done = True
        except ValueError, e:
            print "Invalid Selection..."
            continue
    return selection

def get_config_dict(conf_path):
    config = ConfigObj(conf_path, unrepr=True)
    return config

def get_hcp_env_for_config(conf_dict):
    # exceptions will be raised if config isn't good.
    # hint: we call this from the validation method!
    d = conf_dict["DEFAULT"]
    dp = lambda x: os.path.join(d["hcp_dir"], x)
    new_d = {
        "FSLDIR":d["fsl_dir"],
        "FREESURFER":d["freesurfer_home"],
        "HCPPIPEDIR":d["hcp_dir"],
        "CARET7DIR":dp("global/binaries/caret7/bin_rh_linux64"),
        "HCPPIPEDIR_Templates":dp("global/templates"),
        "HCPPIPEDIR_Bin":dp("global/binaries"),
        "HCPPIPEDIR_Config":dp("global/config"),
        "HCPPIPEDIR_PreFS":dp("PreFreeSurfer/scripts"),
        "HCPPIPEDIR_FS":dp("FreeSurfer/scripts"),
        "HCPPIPEDIR_PostFS":dp("PostFreeSurfer/scripts"),
        "HCPPIPEDIR_fMRISurf":dp("fMRISurface/scripts"),
        "HCPPIPEDIR_fMRIVol":dp("fMRIVolume/scripts"),
        "HCPPIPEDIR_tfMRI":dp("tfMRI/scripts"),
        "HCPPIPEDIR_dMRI":dp("DiffusionPreprocessing/scripts"),
        "HCPPIPEDIR_dMRITract":dp("DiffusionTractography/scripts"),
        "HCPPIPEDIR_Global":dp("global/scripts"),
        "HCPPIPEDIR_tfMRIAnalysis":dp("TaskfMRIAnalysis/scripts"),
        "MSMBin":dp("MSMBinaries"),
        }
    new_d.update(conf_dict.get("env",{}))
    return new_d


def get_hcp_commands_for_config(conf_dict):
    # exceptions will be raised if config isn't good.
    # hint: we call this from the validation method!
    dp = lambda x: os.path.join(conf_dict["DEFAULT"]["hcp_dir"], x)
    return [
        dp("PreFreeSurfer/PreFreeSurferPipeline.sh"),
        dp("FreeSurfer/FreeSurferPipeline.sh"),
        dp("PostFreeSurfer/PostFreeSurferPipeline.sh"),
        dp("fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh"),
        dp("fMRISurface/GenericfMRISurfaceProcessingPipeline.sh"),
        ]

def apply_dict_to_obj(the_d, obj, skip_names=[]):
    if not the_d:
        return
    for name, val in the_d.iteritems():
        if name in skip_names or "traits" not in dir(obj) or name not in obj.traits().keys():
            continue
        setattr(obj, name, val)
