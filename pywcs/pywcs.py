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

    def __init__(self, header=None, key=' ', relax=False, naxis=2):
        if _has_pyfits:
            if isinstance(header, pyfits.NP_pyfits.Header):
                header = str(header.ascardlist())

        WCSBase.__init__(self, header=header, key=key, relax=relax, naxis=naxis)

        if header is None:
            # Set some reasonable defaults.
            self.crpix = numpy.zeros((self.naxis,), numpy.double)
            self.crval = numpy.zeros((self.naxis,), numpy.double)
            self.ctype = ['RA---TAN', 'DEC--TAN']

    def pixel2world(self, *args):
        """
        pixel2world(*args) -> world

        Transforms world coordinates to pixel coordinates.
        L{pixel2world} is a convenience wrapper around L{p2s}.  If
        intermediate values from the transform are needed, use the
        more complex L{p2s} directly.

        Either one or two arguments may be provided.

          - one argument: An Nx2 array of I{x}- and I{y}-coordinates.

          - two arguments: Two one-dimensional arrays of I{x} and I{y}
            coordinates.

        @return: Returns the world coordinates.  If the input was a
            single array, a single array is returned, otherwise a
            tuple of arrays is returned.
        """
        if len(args) == 1:
            return self.p2s(args[0])['world']
        elif len(args) == 2:
            x, y = args
            assert len(x) == len(y)
            length = len(x)
            xy = numpy.hstack((x.reshape((length, 1)),
                               y.reshape((length, 1))))
            world = self.p2s(xy)['world']
            return [world[:, i] for i in range(world.shape[1])]
        raise TypeError("Expected 1 or 2 arguments, %d given" % len(args))

    def world2pixel(self, *args):
        """
        world2pixel(*args) -> pixel

        Transforms world coordinates to pixel coordinates.
        L{world2pixel} is a convenience wrapper around L{s2p}.  If
        intermediate values from the transform are needed, use the
        more complex L{s2p} directly.

        Either one or L{naxis} arguments may be provided.

          - one argument: An NxL{naxis} array of world coordinates.

          - L{naxis} arguments: L{naxis} one-dimensional arrays of
            world coordinates.

        @return: Returns the pixel coordinates.  If the input was a
            single array, a single array is returned, otherwise a
            tuple of arrays is returned.
        """
        if len(args) == 1:
            return self.s2p(args[0])['pixcrd']
        elif len(args) == self.naxis:
            length = len(args[0])
            combined = numpy.hstack([x.reshape((length, 1)) for x in args])
            pixcrd = self.s2p(combined)['pixcrd']
            return pixcrd[:, 0], pixcrd[:, 1]
        raise TypeError("Expected 1 or %d arguments, %d given" %
                        (self.naxis, len(args)))
