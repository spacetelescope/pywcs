# Set the WCS information manually by setting properties of the WCS
# object.

from __future__ import division # confidence high

import numpy
import pywcs
import pyfits
import sys

# Create a new WCS object.  The number of axes must be set
# from the start
wcs = pywcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
wcs.wcs.crpix = [-234.75, 8.3393]
wcs.wcs.cdelt = numpy.array([-0.066667, 0.066667])
wcs.wcs.crval = [0, -90]
wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
wcs.wcs.set_pv([(2, 1, 45.0)])

# Print out the "name" of the WCS, as defined in the FITS header
print wcs.wcs.name

wcs.wcs.print_contents()

# Some pixel coordinates of interest.
pixcrd = numpy.array([[0,0],[24,38],[45,98]], numpy.float_)

# Convert pixel coordinates to world coordinates
world = wcs.wcs_pix2sky(pixcrd, 1)
print world

# Convert the same coordinates back to pixel coordinates.
pixcrd2 = wcs.wcs_sky2pix(world, 1)
print pixcrd2

# These should be the same as the original pixel coordinates, modulo
# some floating-point error.
assert numpy.max(numpy.abs(pixcrd - pixcrd2)) < 1e-6

# Now, write out the WCS object as a FITS header
header = wcs.to_header()

# header is a PyFITS Header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
hdu = pyfits.PrimaryHDU(header=header)
hdu.writeto('test.fits')
