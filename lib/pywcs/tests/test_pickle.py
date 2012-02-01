# TODO: Test that this works for subclasses

import os
try:
    import cPickle as pickle
except ImportError:
    import pickle

import numpy as np
from numpy.testing import assert_array_almost_equal

import pywcs


ROOT_DIR = None
def setup():
    global ROOT_DIR
    ROOT_DIR = os.path.join(os.path.dirname(pywcs.__file__), "tests")


def test_basic():
    wcs = pywcs.WCS()
    s = pickle.dumps(wcs)
    wcs2 = pickle.loads(s)


def test_dist():
    try:
        import pyfits
    except ImportError:
        pass

    hdulist = pyfits.open(os.path.join(ROOT_DIR, "data", "dist.fits"))
    wcs1 = pywcs.WCS(hdulist[0].header, hdulist)
    assert wcs1.det2im2 is not None
    s = pickle.dumps(wcs1)
    wcs2 = pickle.loads(s)

    x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
    world1 = wcs1.all_pix2sky(x, 1)
    world2 = wcs2.all_pix2sky(x, 1)

    assert_array_almost_equal(world1, world2)


def test_sip():
    try:
        import pyfits
    except ImportError:
        pass

    hdulist = pyfits.open(os.path.join(ROOT_DIR, "data", "sip.fits"),
                          ignore_missing_end=True)
    wcs1 = pywcs.WCS(hdulist[0].header)
    assert wcs1.sip is not None
    s = pickle.dumps(wcs1)
    wcs2 = pickle.loads(s)

    x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
    world1 = wcs1.all_pix2sky(x, 1)
    world2 = wcs2.all_pix2sky(x, 1)

    assert_array_almost_equal(world1, world2)


def test_sip2():
    try:
        import pyfits
    except ImportError:
        pass

    hdulist = pyfits.open(os.path.join(ROOT_DIR, "data", "sip2.fits"),
                          ignore_missing_end=True)
    wcs1 = pywcs.WCS(hdulist[0].header)
    assert wcs1.sip is not None
    s = pickle.dumps(wcs1)
    wcs2 = pickle.loads(s)

    x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
    world1 = wcs1.all_pix2sky(x, 1)
    world2 = wcs2.all_pix2sky(x, 1)

    assert_array_almost_equal(world1, world2)


def test_wcs():
    try:
        import pyfits
    except ImportError:
        pass

    filename = os.path.join(ROOT_DIR, "data", "outside_sky.hdr")
    fd = open(filename, 'rb')
    header = fd.read()
    fd.close()

    wcs1 = pywcs.WCS(header)
    s = pickle.dumps(wcs1)
    wcs2 = pickle.loads(s)

    x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
    world1 = wcs1.all_pix2sky(x, 1)
    world2 = wcs2.all_pix2sky(x, 1)

    assert_array_almost_equal(world1, world2)


class Sub(pywcs.WCS):
    def __init__(self, *args, **kwargs):
        self.foo = 42


def test_subclass():
    wcs = Sub()
    s = pickle.dumps(wcs)
    wcs2 = pickle.loads(s)

    assert isinstance(wcs2, Sub)
    assert wcs.foo == 42
    assert wcs2.foo == 42
    assert wcs2.wcs is not None
