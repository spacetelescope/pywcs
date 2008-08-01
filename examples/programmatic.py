# Set the WCS information manually by setting properties of the WCS
# object.

import numpy
import pywcs
import pyfits
import sys

# Create a new WCS object.  The number of axes must be set
# from the start
wcs = pywcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
wcs.crpix = [-234.75, 8.3393]
wcs.cdelt = numpy.array([-0.066667, 0.066667])
wcs.crval = [0, -90]
wcs.ctype = ["RA---AIR", "DEC--AIR"]
wcs.set_pv([(2, 1, 45.0)])

# Print out the "name" of the WCS, as defined in the FITS header
print wcs.name

wcs.print_contents()

# Some pixel coordinates of interest.  These are 0-based coordinates
pixcrd = numpy.array([[0,0],[24,38],[45,98]], numpy.float_)

# Convert pixel coordinates to world coordinates
world = wcs.pixel2world(pixcrd)
print world

# Convert the same coordinates back to pixel coordinates.
pixcrd2 = wcs.world2pixel(world)
print pixcrd2

# These should be the same as the original pixel coordinates, modulo
# some floating-point error.
assert numpy.max(numpy.abs(pixcrd - pixcrd2)) < 1e-6

# Now, write out the WCS object as a FITS header
header = wcs.to_header()

# header is a PyFITS header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
hdu = pyfits.PrimaryHDU(header=header)
hdu.writeto('test.fits')
