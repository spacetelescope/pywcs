# Copyright (C) 2008 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

__docformat__ = "epytext"

import numpy

import _pywcs
try:
    import pyfits
    _has_pyfits = True
except ImportError:
    _has_pyfits = False

# This is here for the sake of epydoc
WCSBase = _pywcs._WCS

# A wrapper around the C WCS type
class WCS(WCSBase):
    """%s""" % _pywcs._WCS.__doc__

    def __init__(self, header, key=' ', relax=False):
        if _has_pyfits:
            if isinstance(header, pyfits.NP_pyfits.Header):
                header = str(header.ascardlist())
        WCSBase.__init__(self, header, key, relax)

    def pixel2world(self, x, y=None):
        """
        pixel2world(x, y=None) -> world

        Transforms world coordinates to pixel coordinates.
        L{pixel2world} is a convenience wrapper around L{p2s}.  If
        intermediate values from the transform are needed, use the
        more complex L{p2s} directly.

        @param x: Array of I{x} pixel coordinates, or if L{y} is not
            provided, a 2-dimensional array of both I{x}- and
            I{y}-coordinates.  B{Pixel coordinates are zero based.}
        @type x: array of double

        @param y: Array of I{y} pixel coordinates.  If provided, L{x}
            must be 1-dimensional and the same length as L{y}.
        @type y: array of double or None

        @return: If both L{x} and L{y} were provided, a 2-tuple of
            arrays, where the first element is latitude coordinates
            and the second element is longitude coordinates.
            Otherwise, a single 2D array containing both latitude and
            longitude.
        """
        if y is None:
            return self.p2s(x)['world']
        else:
            assert len(x) == len(y)
            length = len(x)
            xy = numpy.hstack((x.reshape((length, 1)), y.reshape((length, 1))))
            world = self.p2s(xy)['world']
            return world[:, 0], world[:, 1]

    def world2pixel(self, ra, dec=None):
        """
        world2pixel(ra, dec=None) -> pixel

        Transforms world coordinates to pixel coordinates.
        L{world2pixel} is a convenience wrapper around L{s2p}.  If
        intermediate values from the transform are needed, use the
        more complex L{s2p} directly.

        @param ra: Array of I{ra} world coordinates, or if L{dec} is
            not provided, a 2-dimensional array of both I{ra}- and
            I{dec}-coordinates.  B{World coordinates are in decimal
            degrees.}
        @type ra: array of double

        @param dec: Array of I{dec} world coordinates.  If provided,
            L{ra} must be 1-dimensional and the same length as L{dec}.
        @type dec: array of double or None

        @return: If both L{ra} and L{dec} were provided, a 2-tuple of
            arrays, where the first element is I{x} coordinates and
            the second element is I{y} coordinates.  Otherwise, a
            single 2D array containing both I{x} and I{y}.  B{Pixel
            coordinates are zero-based.}
        """
        if dec is None:
            return self.p2s(x)['pixcrd']
        else:
            assert len(ra) == len(dec)
            length = len(ra)
            radec = numpy.hstack((ra.reshape((length, 1)), dec.reshape((length, 1))))
            pixcrd = self.p2s(radec)['pixcrd']
            return pixcrd[:, 0], pixcrd[:, 1]
