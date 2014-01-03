#!/usr/bin/env python

import sys
import os

# put this module at the front of the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.split(__file__)[0], "../../")))

from hcp_prep.interfaces import *

out_file = os.path.join(os.path.split(__file__)[0], "interface_docs.txt")

ints = [HCDcm2nii, DicomInfo, NiiWrangler, HCPCommand, PreFS, FS, PostFS, VolumeProcessing, SurfaceProcessing]

with open(out_file, 'w') as f:
    for c in ints:
        print >>f, "%s\n%s\n%s" % (c.__name__, "="*len(c.__name__), c.help(returnhelp=True))