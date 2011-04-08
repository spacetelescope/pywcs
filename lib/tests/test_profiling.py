import glob
import os
import sys

import numpy as np
from numpy.testing import assert_array_almost_equal

import pywcs

ROOT_DIR = None
def setup():
    global ROOT_DIR
    ROOT_DIR = os.path.dirname(__file__)

def test_maps():
    def test_map(filename):
        fd = open(filename, 'rb')
        header = fd.read()
        fd.close()
        wcs = pywcs.WCS(header)

        x = np.random.rand(2 ** 16, wcs.wcs.naxis)
        world = wcs.wcs_pix2sky(x, 1)
        pix = wcs.wcs_sky2pix(x, 1)

    for filename in glob.glob(os.path.join(ROOT_DIR, "maps", "*.hdr")):
        yield test_map, filename

def test_spectra():
    def test_spectrum(filename):
        fd = open(filename, 'rb')
        header = fd.read()
        fd.close()
        wcs = pywcs.WCS(header)

        x = np.random.rand(2 ** 16, wcs.wcs.naxis)
        world = wcs.wcs_pix2sky(x, 1)
        pix = wcs.wcs_sky2pix(x, 1)

    for filename in glob.glob(os.path.join(ROOT_DIR, "spectra", "*.hdr")):
        yield test_spectrum, filename
