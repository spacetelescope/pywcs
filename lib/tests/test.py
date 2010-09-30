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

basic_units = "m s g rad sr K A mol cd".split()
derived_units = "Hz J W V N Pa C Ohm ohm S F Wb T H lm lx".split()
add_all_units = "eV Jy R G barn".split()
add_sup_units = "a yr pc bit byte Byte".split()
add_sub_units = "mag".split()
general_units = "deg arcmin arcsec mas d h min erg Ry u D".split()
astro_units = "Angstrom angstrom AU lyr beam solRad solMass solLum Sun".split()
device_units = "adu bin chan count ct photon ph pixel pix voxel".split()
sub_prefixes = "y z a f p n u m c d".split()
sup_prefixes = "da h k M G T P E Z Y".split()
def test_all_units():
    def test_self(x):
        try:
            u = pywcs.UnitConverter(x, x)
        except ValueError, e:
            if str(e) == "Potentially unsafe translation" and \
                    x in ("S", "H", "D"):
                return
            else:
                raise
        assert u.scale == 1.0
        assert u.offset == 0.0
        assert u.power == 1.0

    for unit in (basic_units + derived_units + add_all_units + add_sup_units +
                 add_sub_units + general_units + astro_units + device_units):
        yield test_self, unit

def test_unit_prefixes():
    def test_self(x, p):
        unit = p + x
        try:
            u = pywcs.UnitConverter(unit, unit)
        except ValueError, e:
            if str(e) == "Potentially unsafe translation" and \
                    x in ("S", "H", "D"):
                return
            else:
                raise
        assert u.scale == 1.0
        assert u.offset == 0.0
        assert u.power == 1.0

    for unit in (basic_units + derived_units + add_all_units):
        for prefix in (sub_prefixes + sup_prefixes):
            yield test_self, unit, prefix

    for unit in add_sup_units:
        for prefix in sup_prefixes:
            yield test_self, unit, prefix

    for unit in add_sub_units:
        for prefix in sub_prefixes:
            yield test_self, unit, prefix

