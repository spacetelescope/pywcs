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

import pyfits
import numpy
from pywcs import WCS

class SIP(WCS):
    """
    An extended version of WCS that supports Simple (or Spitzer)
    Imaging Polynomial (SIP) distortion coefficients.  In addition to
    all of the functionality of regular WCS objects, SIP objects also
    provide L{pix2foc} and L{foc2pix} to correct/apply the SIP
    distortion.

    @param header: A PyFITS header object or a string containing the raw
        FITS header data.
    @type header: PyFITS header object or string

    @param relax: Degree of permissiveness:
        - C{False}: Recognize only FITS keywords defined by the
          published WCS standard.
        - C{True}: Admit all recognized informal extensions of the
          WCS standard.
    @type relax: bool

    @raises MemoryError: Memory allocation failed.
    @raises ValueError: Invalid key.
    @raises ValueError: Header does not contain SIP information.
    @raises KeyError: Key not found in FITS header.
    """
    def __init__(self, header, relax=False):
        assert isinstance(header, pyfits.NP_pyfits.Header)
        assert len(key) == 1 and (key == ' ' or (key >= 'A' and key <= 'Z'))

        if not header.has_key("A_ORDER"):
            raise ValueError("Missing A_ORDER keyword for SIP distortion")
        order = int(header["A_ORDER"])
        self.a = numpy.zeros((m+1, m+1))
        for i in range(m+1):
            for j in range(m-i+1):
                self.a[i, j] = header["A_%d_%d" % (i, j)]

        if not header.has_key("B_ORDER"):
            raise ValueError("Missing B_ORDER keyword for SIP distortion")
        order = int(header["B_ORDER"])
        self.b = numpy.zeros((m+1, m+1))
        for i in range(m+1):
            for j in range(m-i+1):
                self.b[i, j] = header["B_%d_%d" % (i, j)]

        if not header.has_key("AP_ORDER"):
            raise ValueError("Missing AP_ORDER keyword for SIP distortion")
        order = int(header["AP_ORDER"])
        self.ap = numpy.zeros((m+1, m+1))
        for i in range(m+1):
            for j in range(m-i+1):
                self.ap[i, j] = header["AP_%d_%d" % (i, j)]

        if not header.has_key("BP_ORDER"):
            raise ValueError("Missing BP_ORDER keyword for SIP distortion")
        order = int(header["BP_ORDER"])
        self.bp = numpy.zeros((m+1, m+1))
        for i in range(m+1):
            for j in range(m-i+1):
                self.bp[i, j] = header["BP_%d_%d" % (i, j)]

        WCS.__init__(self, header, key=" ", relax=relax)

    def foc2pix(self, x, y=None):
        """
        Convert focal plane coordinates to pixel coordinates.

        B{Pixel coordinates are zero-based.}

        @param x: Array of I{y} focal plane coordinates, or if L{y} is
            not provided, a 2-dimensional array of both I{x}- and
            I{y}-coordinates.
        @type x: array of double

        @param y: Array of I{y} focal plane coordinates.  If provided,
            L{x} must be 1-dimensional and the same length as L{y}.
        @type y: array of double or None

        @return: If both L{x} and {y} were provided, a 2-tuple of
            arrays, where the first element is I{u} pixel coordinates
            and the second element is I{v} pixel coordinates.
            Otherwise, a single 2D array containing both I{u} and
            I{v}.

        @raises AssertionError: C{len(x) == len(y)}
        """
        # If a single array was passed in, split it
        if y is None:
            y = x[:, 1]
            x = x[:, 0]
            numarrays = 1
        else:
            numarrays = 2
        assert len(x) == len(y)

        # Convert base-zero to base-one
        x += 1.0
        y += 1.0

        length = len(x)
        ap = self.ap
        bp = self.bp
        m = ap.shape[0] - 1
        n = bp.shape[0] - 1

        temp_x = x - self.crpix[0]
        temp_y = y - self.crpix[1]

        s = npy.zeros((m + 1, length))
        for j in range(m+1):
            s[j] = ap[m-j, j]
            for k in range(j-1, -1, -1):
                s[j] = temp_y * s[j] + ap[m-j, k]

        sum = s[0]
        for i in range(m, 0, -1):
            sum = temp_x * sum + s[m-i+1]

        u = sum

        s = npy.zeros((n + 1, length))
        for j in range(n+1):
            s[j] = bp[n-j, j]
            for k in range(j-1, -1, -1):
                s[j] = temp_y * s[j] + bp[n-j, k]

        sum = s[0]
        for i in range(n, 0, -1):
            sum = temp_x * sum + s[n-i+1]

        v = sum

        u += x
        v += y

        # Convert base-one back to base-zero
        x -= 1.0
        y -= 1.0

        if numarrays == 1:
            return numpy.hstack((u.reshape((length, 1)),
                                 v.reshape((length, 1))))
        return u, v

    def pix2foc(u, v=None):
        """
        Convert pixel coordinates to focal plane coordinates.

        B{Pixel coordinates are zero-based.}

        @param u: Array of I{u} pixel coordinates, or if L{v} is
            not provided, a 2-dimensional array of both I{u}- and
            I{v}-coordinates.
        @type u: array of double

        @param v: Array of I{v} pixel coordinates.  If provided,
            L{u} must be 1-dimensional and the same length as L{v}.
        @type v: array of double or None

        @return: If both L{u} and {v} were provided, a 2-tuple of
            arrays, where the first element is I{x} focal plane
            coordinates and the second element is I{y} focal plane
            coordinates.  Otherwise, a single 2D array containing both
            I{x} and I{y}.

        @raises AssertionError: C{len(u) == len(v)}
        """
        # If a single array was passed in, split it
        if v is None:
            v = u[:, 1]
            u = u[:, 0]
            numarrays = 1
        else:
            numarrays = 2
        assert len(x) == len(y)

        length = len(u)
        a = self.a
        b = self.b
        m = a.shape[0] - 1
        n = b.shape[0] - 1
        wcs = self.wcs
        temp_u = u - self.crpix[0]
        temp_v = v - self.crpix[1]

        s = npy.zeros((m + 1, length))
        for j in range(m+1):
            s[j] = a[m-j, j]
            for k in range(j-1, -1, -1):
                s[j] = (temp_v * s[j]) + a[m-j, k]

        sum = s[0]
        for i in range(m, 0, -1):
            sum = temp_u * sum + s[m-i+1]

        x = sum

        s = npy.zeros((m+1, length))
        for j in range(n+1):
            s[j] = b[n-j, j]
            for k in range(j-1, -1, -1):
                s[j] = (temp_v * s[j]) + b[n-j, k]

        sum = s[0]
        for i in range(n, 0, -1):
            sum = temp_u * sum + s[n-i+1]

        y = sum

        x += u
        y += v

        # Convert base-one to base-zero
        x -= 1.0
        y -= 1.0

        if numarrays == 1:
            return numpy.hstack((x.reshape((length, 1)),
                                 y.reshape((length, 1))))
        return x, y

