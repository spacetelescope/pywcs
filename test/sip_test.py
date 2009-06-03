from __future__ import division # confidence high

import pywcs
import pyfits
import numpy
import sys

hdulist = pyfits.open(sys.argv[-1])

data1 = numpy.array([0,2,4,6])
data2 = numpy.array([1,3,5,7])

header = hdulist[1].header
wcs = pywcs.WCS(header)
assert wcs.sip is not None

print data1, data2
focal = wcs.sip_pix2foc(data1, data2)
print focal

focal2 = wcs.pix2foc(data1, data2)
print focal
