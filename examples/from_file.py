# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.

import numpy
import pywcs
import pyfits
import sys

hdulist = pyfits.open(sys.argv[-1])

# Parse the WCS keywords in the primary HDU
wcs = pywcs.WCS(hdulist[0].header)

# Print out the "name" of the WCS, as defined in the FITS header
print wcs.name
print wcs.get_pv()

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
