def orientation_from_dcm_header(header):
    # borrowed from hcp xnat dicom juggling java
    # get and check the base values
    if not header:
        raise ValueError("didn't get a header")
    o = getattr(header, "ImageOrientationPatient", None)
    o = [float(a) for a in o]
    if not o:
        raise ValueError("couldn't find ImageOrientationPatient in header")
    if len(o) != 6:
        raise ValueError("cannot be translated to cosine vectors")
    # consistency checks
    epsilon = 0.001
    if abs(o[0] * o[3] + o[1] * o[4] + o[2] * o[5]) > 0.001:
        raise ValueError("cosine vectors not orthogonal")
    if abs(1.0 - o[0] * o[0] - o[1] * o[1] - o[2] * o[2]) > epsilon:
        raise ValueError("cosine vectors not normal")
    # looks like we're good to go. derive the value
    absNormalX = abs(o[1] * o[5] - o[2] * o[4])
    absNormalY = abs(o[2] * o[3] - o[0] * o[5])
    absNormalZ = abs(o[0] * o[4] - o[1] * o[3])
    if absNormalX > absNormalY:
        return "sagittal" if absNormalX > absNormalZ else "transverse"
    else:
        return "coronal" if absNormalY > absNormalZ else "transverse"

def numberfy(s):
    n = s
    try:
        n = float(n)
        return n
    except Exception:
        return s

def float_or_none(s):
    n = s
    try:
        n = float(n)
        return n
    except Exception:
        return None

def int_or_none(s):
    n = s
    try:
        n = int(n)
        return n
    except ValueError:
        return None