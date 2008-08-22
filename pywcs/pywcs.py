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

if hasattr(_pywcs, "Distortion"):
    Distortion = _pywcs.Distortion
DistortionLookupTable = _pywcs.DistortionLookupTable

# A wrapper around the C WCS type
class WCS(WCSBase):
    """%s""" % _pywcs._WCS.__doc__

    def __init__(self, header=None, fobj=None, key=' ', relax=False, naxis=2):
        if isinstance(header, pyfits.NP_pyfits.Header):
            WCSBase.__init__(self, header=str(header.ascardlist()), key=key, relax=relax, naxis=naxis)
            self.addDistortion(header, fobj)
        else:        
            WCSBase.__init__(self, header=header, key=key, relax=relax, naxis=naxis)
        
        if header is None:
            # Set some reasonable defaults.
            self.crpix = numpy.zeros((self.naxis,), numpy.double)
            self.crval = numpy.zeros((self.naxis,), numpy.double)
            self.ctype = ['RA---TAN', 'DEC--TAN']
            self.cd = numpy.array([[1.0, 0.0], [0.0, 1.0]], numpy.double)
            
    def addDistortion(self, header=None, fobj=None):
        """
        Adds distortion information as per WCS paper IV.
        Only lookup table distortion is omplemented sofar.
        """
        assert isinstance(header, pyfits.NP_pyfits.Header)
        d_keys = []
        for i in range(1+self.naxis+1):
            d_keys.append('CPDIS'+str(i))
            d_keys.append('CQDIS'+str(i))
        d_val = numpy.array([header.has_key(k) for k in d_keys])
        if d_val.any():
        
            for dist in ['CPDIS', 'CQDIS']:
                self.readDistortionKw(header, fobj, dist)
            
    def readDistortionKw(self, header, fobj, dist=""):
        """
        Reads paper IV distortion keywords and data.
        """
        if dist == 'CPDIS':
            d_kw= 'DP'
        else:
            d_kw = 'DQ'
            
        for i in range(1, self.naxis+1):
            distortion = dist+str(i)
            if header.has_key(distortion):
                dis = header[distortion].lower()
                if dis == 'lookup':
                    assert isinstance(fobj, pyfits.NP_pyfits.HDUList), 'A pyfits HDUList is required for Lookup table distrtion.'
                    dp = d_kw+str(i)
                    d_extver = header[dp+'.EXTVER']
                    d_data = fobj['WCSDVARR', d_extver].data
                    d_header = fobj['WCSDVARR', d_extver].header
                    d_crpix = (d_header['CRPIX1'], d_header['CRPIX2'])
                    d_crval = (d_header['CRVAL1'], d_header['CRVAL2'])
                    d_cdelt = (d_header['CDELT1'], d_header['CDELT2'])
                    d_lookup = DistortionLookupTable(d_data, d_crpix, d_crval, d_cdelt)
                    setattr(self, distortion.lower(), d_lookup)
                else:
                    print 'Polynomial distortion is not implemented.\n'
            else:
                pass
                
    
    def _pixel2world_generic(self, func, *args):
        if len(args) == 1:
            return func(args[0])['world']
        elif len(args) == 2:
            x, y = args
            assert len(x) == len(y)
            length = len(x)
            xy = numpy.hstack((x.reshape((length, 1)),
                               y.reshape((length, 1))))
            world = func(xy)['world']
            return [world[:, i] for i in range(world.shape[1])]
        raise TypeError("Expected 1 or 2 arguments, %d given" % len(args))

    def pixel2world(self, *args):
        """
        pixel2world(*args) -> world

        Transforms world coordinates to pixel coordinates.

        L{pixel2world} is a convenience wrapper around L{p2s}.  If
        intermediate values from the transform are needed, use the
        more complex L{p2s} directly.

        B{The pixel coordinates given are 0-based (like array indices
        in C and Python).  If your pixel coordinates are 1-based (like
        array indices in Fortran), use L{pixel2world_fits} instead.}

        Either one or two arguments may be provided.

          - one argument: An Nx2 array of I{x}- and I{y}-coordinates.

          - two arguments: Two one-dimensional arrays of I{x} and I{y}
            coordinates.

        @return: Returns the world coordinates.  If the input was a
            single array, a single array is returned, otherwise a
            tuple of arrays is returned.
        """
        return self._pixel2world_generic(self.p2s, *args)

    def pixel2world_fits(self, *args):
        """
        pixel2world_fits(*args) -> world

        Identical to L{pixel2world}, except pixel coordinates are
        1-based (like array indices in Fortran), instead of 0-based
        (like array indices C and Python).
        """
        return self._pixel2world_generic(self.p2s_fits, *args)

    def _world2pixel_generic(self, func, *args):
        if len(args) == 1:
            return func(args[0])['pixcrd']
        elif len(args) == self.naxis:
            length = len(args[0])
            combined = numpy.hstack([x.reshape((length, 1)) for x in args])
            pixcrd = func(combined)['pixcrd']
            return pixcrd[:, 0], pixcrd[:, 1]
        raise TypeError("Expected 1 or %d arguments, %d given" %
                        (self.naxis, len(args)))

    def world2pixel(self, *args):
        """
        world2pixel(*args) -> pixel

        Transforms world coordinates to pixel coordinates.
        L{world2pixel} is a convenience wrapper around L{s2p}.  If
        intermediate values from the transform are needed, use the
        more complex L{s2p} directly.

        B{The pixel coordinates returned are 0-based (like array
        indices in C and Python).  If you require pixel coordinates to
        be 1-based (like array indices in Fortran), use
        L{world2pixel_fits} instead.}

        Either one or L{naxis} arguments may be provided.

          - one argument: An NxL{naxis} array of world coordinates.

          - L{naxis} arguments: L{naxis} one-dimensional arrays of
            world coordinates.

        @return: Returns the pixel coordinates.  If the input was a
            single array, a single array is returned, otherwise a
            tuple of arrays is returned.
        """
        return self._world2pixel_generic(self.s2p, *args)

    def world2pixel_fits(self, *args):
        """
        pixel2world_fits(*args) -> world

        Identical to L{world2pixel}, except pixel coordinates are
        1-based (like array indices in Fortran), instead of 0-based
        (like array indices C and Python).
        """
        return self._world2pixel_generic(self.s2p_fits, *args)

    if _has_pyfits:
        def to_header(self, relax=False):
            header_string = WCSBase.to_header(self, relax)
            cards = pyfits.CardList()
            for i in range(0, len(header_string), 80):
                card_string = header_string[i:i+80]
                card = pyfits.Card()
                card.fromstring(card_string)
                cards.append(card)
            return pyfits.Header(cards)
    else:
        def to_header(self, *args, **kwargs):
            raise NotImplementedError(
                "PyFITS must be installed to generate a FITS header")

    to_header.__doc__ = """
    Generate a PyFITS header object with the WCS information stored in
    this object.

    The output header will almost certainly differ from the input in a
    number of respects:

      1. The output header only contains WCS-related keywords.  In
         particular, it does not contain syntactically-required
         keywords such as C{SIMPLE}, C{NAXIS}, C{BITPIX}, or C{END}.

      2. Deprecated (e.g. C{CROTAn}) or non-standard usage will be
         translated to standard (this is partially dependent on
         whether L{fix} was applied).

      3. Quantities will be converted to the units used internally,
         basically SI with the addition of degrees.

      4. Floating-point quantities may be given to a different decimal
         precision.

      5. Elements of the C{PCi_j} matrix will be written if and only
         if they differ from the unit matrix.  Thus, if the matrix is
         unity then no elements will be written.

      6. Additional keywords such as C{WCSAXES}, C{CUNITia},
         C{LONPOLEa} and C{LATPOLEa} may appear.

      7. The original keycomments will be lost, although L{to_header}
         tries hard to write meaningful comments.

      8. Keyword order may be changed.

    @param relax: Degree of permissiveness:

      - C{False}: Recognize only FITS keywords defined by the
        published WCS standard.

      - C{True}: Admit all recognized informal extensions of the WCS
        standard.

    @type relax: bool

    @return: A PyFITS Header object.
    """

    def to_header_string(self, relax=False):
        return WCSBase.to_header(self, relax)
    to_header_string.__doc__ = WCSBase.to_header.__doc__
