from __future__ import division # confidence high

import pywcs
import pyfits
import numpy
import sys

hdulist = pyfits.open(sys.argv[-1])

data1 = numpy.array([0,2,4,6])
data2 = numpy.array([1,3,5,7])

header = hdulist[1].header
wcs = pywcs.WCS(header, hdulist)
print wcs.cpdis1, wcs.cpdis2
wcs.cpdis2 = None
wcs.cpdis1.data = numpy.array([[1.0], [2.0], [2.0]], numpy.float32) # numpy.ones((1, 1), numpy.float32) * 5.0
print "crval %s crpix %s cdelt %s" % (wcs.cpdis1.crval, wcs.cpdis1.crpix, wcs.cpdis1.cdelt)
wcs.cpdis1.crval = (1.0, 0.0)
wcs.cpdis1.crpix = (1.0, 1.0)

print "Just P4"
focal = wcs.p4_pix2foc(data1, data2, 1)
print focal

print "Just SIP"
focal2 = wcs.sip_pix2foc(data1, data2, 0)
print focal2

print "Both"
focal3 = wcs.pix2foc(data1, data2, 0)
print focal3
