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

members = """cd cdelt crota crpix crval ctype cubeface cunit lat latpole lng lonpole naxis pc ps pv restfrq restwav spec""".split()

def test_file(path):
    print "=" * 75
    print path
    hdulist = pyfits.open(path)
    wcslist = pywcs.parse_hdulist(hdulist)
    data1 = numpy.array([100,200])
    data2 = numpy.array([200,300])
    data3 = numpy.array([[0,1],[2,3],[4,5],[6,7]], numpy.float_)
    for wcs in wcslist:
        wcs.fix()
        wcs.set()
        print "p2s: %s" % wcs.p2s(data3)
        print "s2p: %s" % wcs.s2p(data3)
        try:
            print "mix: %s" % wcs.mix(1, 1, (-120,120), 0.0, 10, data1, data2)
        except ValueError, e:
            print "mix: %s" % e
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
