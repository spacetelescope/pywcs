#!/usr/bin/env python

# This is really more of a smoke-test than a unit test.
# It just checks that a set of test FITS files load and
# that all methods can be called on the resulting objects.
# There is no verification of the values.  That should be
# done eventually, but also keep in mind that the wrapper
# itself does no math.

import glob
import os
import sys

import numpy

import pywcs
import pyfits

members = """alt cd cdelt cname colax colnum crder crota crpix crval csyer ctype cubeface cunit dateavg dateobs equinox lat latpole lng lonpole mjdavg mjdobs name naxis obsgeo pc radesys restfrq restwav spec specsys ssysobs ssyssrc velosys velangl zsource""".split()

def test_file(path):
    print "=" * 75
    print path
    hdulist = pyfits.open(path)
    wcs = pywcs.WCS(hdulist[0].header)

    data1 = numpy.array([0,2,4,6])
    data2 = numpy.array([1,3,5,7])
    data3 = numpy.array([[0,1],[2,3],[4,5],[6,7]], numpy.float_)

    wcs.fix()
    wcs.set()
    print "p2s: %s" % wcs.p2s(data3)
    print "pixel2world: %s %s" % wcs.pixel2world(data1, data2)
    print "s2p: %s" % wcs.s2p(data3)

#         try:
#             print "mix: %s" % wcs.mix(1, 1, (-120,120), 0.0, 10, data1, data2)
#         except ValueError, e:
#             print "mix: %s" % e
    print "has_cdi_ja: %s" % wcs.has_cdi_ja()
    print "has_crotaia: %s" % wcs.has_crotaia()
    print "has_pci_ja: %s" % wcs.has_pci_ja()
    for member in members:
        try:
            val = getattr(wcs, member)
        except Exception, e:
            print "wcs.%s: EXCEPTION: %s" % (member, e)
        else:
            print "wcs.%s: %s" % (member, val)
    wcs.print_contents()


def run_directory(directory):
    for filepath in glob.glob(os.path.join(directory, "*.fits")):
        test_file(filepath)

if __name__ == '__main__':
    directory = sys.argv[-1]
    if not os.path.exists(directory):
        print "usage: test.py [directory of fits files]"
        sys.exit(1)

    run_directory(directory)
