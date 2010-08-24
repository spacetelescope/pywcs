import glob
import os

import numpy as np
import pyfits
from numpy.testing import assert_array_almost_equal

import pywcs

ROOT_DIR = None
def setup():
    global ROOT_DIR
    ROOT_DIR = os.path.dirname(__file__)

def test_maps():
    def test_map(filename):
        hdulist = pyfits.open(filename)
        wcs = pywcs.WCS(hdulist[0].header)

        world = wcs.wcs_pix2sky([[97, 97]], 1)

        assert_array_almost_equal(world, [[285.0, -66.25]], decimal=1)

    for filename in glob.glob(os.path.join(ROOT_DIR, "maps", "*.fits")):
        yield test_map, filename

def test_spectra():
    def test_spectrum(filename):
        hdulist = pyfits.open(filename)
        wcs = pywcs.WCS(hdulist[0].header)

        all = pywcs.find_all_wcs(hdulist[0].header)
        assert len(all) == 9

    for filename in glob.glob(os.path.join(ROOT_DIR, "spectra", "*.fits")):
        yield test_spectrum, filename

def test_units():
    u = pywcs.UnitConverter("log(MHz)", "ln(Hz)")
    print u.convert([1,2,3,4])
