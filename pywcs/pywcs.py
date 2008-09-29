# Copyright (C) 2008 Association of Universities for Research in
# Astronomy (AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#     1. Redistributions of source code must retain the above
#       copyright notice, this list of conditions and the following
#       disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials
#       provided with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
"""
Pywcs provides transformations following the SIP conventions, Paper IV
table lookup distortion, and the core WCS functionality provided by
wcslib.  Each of these transformations can be used independently or
together in a standard pipeline.

The basic workflow is as follows:

    1. C{import pywcs}

    2. Call the C{pywcs.WCS} constructor with a PyFITS header object.

    3. Optionally, if the FITS file uses any deprecated or
       non-standard features, you may need to call one of the C{fix}
       methods on the object.

    4. Use one of the following transformation methods:

       - all_pix2sky: Perform all three transformations from pixel to
         sky coords.

       - wcs_pix2sky: Perform just the core WCS transformation from
         pixel to sky coords.

       - wcs_sky2pix: Perform just the core WCS transformation from
         sky to pixel coords.

       - sip_pix2foc: Convert from pixel to focal plane coords using
         the SIP polynomial coefficients.

       - sip_foc2pix: Convert from focal plane to pixel coords using
         the SIP polynomial coefficients.

       - p4_pix2foc: Convert from pixel to focal plane coords using
         the table lookup distortion method described in Paper IV.
"""

__docformat__ = "epytext"

import numpy

import _docutil as __
import _pywcs
import pyfits

# This is here for the sake of epydoc
WCSBase = _pywcs._Wcs
DistortionLookupTable = _pywcs.DistortionLookupTable
Sip = _pywcs.Sip
Wcsprm = _pywcs._Wcsprm

# A wrapper around the C WCS type
class WCS(WCSBase):
    """
    WCS objects correct for SIP and Paper IV table-lookup distortions
    and transform pixel to/from sky coordinates, based on the WCS
    keywords and data in a FITS file.
    """

    def __init__(self, header=None, fobj=None, key=' ', relax=False, naxis=2):
        """
        WCS(header=None, fobj=None, key=' ', relax=False, naxis=2)

        @param header: A PyFITS header object.  If header is not
            provided, the object will be initialized to default
            values.
        @type header: PyFITS header object

        @param fobj: A PyFITS file object. It is needed when header
            keywords point to a Lookup table distortion stored in a
            different extension as per WCS paper IV.
        @type fobj: PyFITS HDUList object

        @param key: The key referring to a particular WCS transform in
           the header.  This may be either C{' '} or C{'A'}-C{'Z'} and
           corresponds to the C{"a"} part of C{"CTYPEia"}.  (C{key}
           may only be provided if C{header} is also provided.)
        @type key: string

        @param relax: Degree of permissiveness:
            - C{False}: Recognize only FITS keywords defined by the
              published WCS standard.
            - C{True}: Admit all recognized informal extensions of the
              WCS standard.
        @type relax: bool

        @param naxis: The number of sky coordinates axes for the
            object.  (C{naxis} may only be provided if C{header} is
            C{None}.)  The only number of axis currently supported is 2.
        @type naxis: int

        @raises MemoryError: Memory allocation failed.
        @raises ValueError: Invalid key.
        @raises KeyError: Key not found in FITS header.
        @raises AssertionError: Lookup table distortion present in the
            header but fobj not provided.
        """
        if naxis != 2:
            raise ValueError("Only 2 axes are supported")
        self.naxis = naxis

        if header is None:
            wcsprm = _pywcs._Wcsprm(header=None, key=key,
                                    relax=relax, naxis=naxis)
            # Set some reasonable defaults.
            wcsprm.crpix = numpy.zeros((self.naxis,), numpy.double)
            wcsprm.crval = numpy.zeros((self.naxis,), numpy.double)
            wcsprm.ctype = ['RA---TAN', 'DEC--TAN']
            wcsprm.cd = numpy.array([[1.0, 0.0], [0.0, 1.0]], numpy.double)
            cpdis = (None, None)
        else:
            wcsprm = _pywcs._Wcsprm(header=str(header.ascardlist()), key=key,
                                    relax=relax, naxis=naxis)
            cpdis = self._read_distortion_kw(header, fobj, dist='CPDIS')
            sip = self._read_sip_kw(header)

        WCSBase.__init__(self, sip, cpdis, wcsprm)

    def _read_distortion_kw(self, header, fobj, dist='CPDIS'):
        """
        Reads paper IV table-lookup distortion keywords and data, and
        returns a 2-tuple of DistortionLookupTable objects.

        If no Paper IV distortion keywords are found, (None, None) is
        returned.
        """
        if dist == 'CPDIS':
            d_kw = 'DP'
        else:
            d_kw = 'DQ'

        tables = {}
        for i in range(1, self.naxis+1):
            distortion = dist+str(i)
            if header.has_key(distortion):
                dis = header[distortion].lower()
                if dis == 'lookup':
                    assert isinstance(fobj, pyfits.NP_pyfits.HDUList), \
                    'A pyfits HDUList is required for Lookup table distortion.'
                    dp = d_kw+str(i)
                    d_extver = header[dp+'.EXTVER']
                    d_data = fobj['WCSDVARR', d_extver].data
                    d_header = fobj['WCSDVARR', d_extver].header
                    d_crpix = (d_header['CRPIX1'], d_header['CRPIX2'])
                    d_crval = (d_header['CRVAL1'], d_header['CRVAL2'])
                    d_cdelt = (d_header['CDELT1'], d_header['CDELT2'])
                    d_lookup = DistortionLookupTable(d_data, d_crpix,
                                                     d_crval, d_cdelt)
                    tables[i] = d_lookup
                else:
                    print 'Polynomial distortion is not implemented.\n'
            else:
                pass

        if len(tables) == 2:
            return (tables[0], tables[1])
        return (None, None)

    def _read_sip_kw(self, header):
        """
        Reads SIP header keywords and returns a Sip object.

        If no SIP header keywords are found, None is returned.
        """
        if header.has_key("A_ORDER"):
            if not header.has_key("B_ORDER"):
                raise ValueError(
                    "A_ORDER provided without corresponding B_ORDER "
                    "keyword for SIP distortion")

            m = int(header["A_ORDER"])
            a = numpy.zeros((m+1, m+1))
            for i in range(m+1):
                for j in range(m-i+1):
                    a[i, j] = header.get("A_%d_%d" % (i, j), 0.0)

            m = int(header["B_ORDER"])
            b = numpy.zeros((m+1, m+1))
            for i in range(m+1):
                for j in range(m-i+1):
                    b[i, j] = header.get("B_%d_%d" % (i, j), 0.0)
        elif header.has_key("B_ORDER"):
            raise ValueError(
                "B_ORDER provided without corresponding A_ORDER "
                "keyword for SIP distortion")
        else:
            a = None
            b = None

        if header.has_key("AP_ORDER"):
            if not header.has_key("BP_ORDER"):
                raise ValueError(
                    "AP_ORDER provided without corresponding BP_ORDER "
                    "keyword for SIP distortion")

            m = int(header["AP_ORDER"])
            ap = numpy.zeros((m+1, m+1))
            for i in range(m+1):
                for j in range(m-i+1):
                    ap[i, j] = header.get("AP_%d_%d" % (i, j), 0.0)

            m = int(header["BP_ORDER"])
            bp = numpy.zeros((m+1, m+1))
            for i in range(m+1):
                for j in range(m-i+1):
                    bp[i, j] = header.get("BP_%d_%d" % (i, j), 0.0)
        elif header.has_key("BP_ORDER"):
            raise ValueError(
                "BP_ORDER provided without corresponding AP_ORDER "
                "keyword for SIP distortion")
        else:
            ap = None
            bp = None

        if a is None and b is None and ap is None and bp is None:
            return None

        return Sip(a, b, ap, bp, self.wcs.crpix)

    def _array_converter(self, func, *args):
        if len(args) == 1:
            return func(args[0])
        elif len(args) == 2:
            x, y = args
            if len(x) != len(y):
                raise ValueError("x and y arrays are not the same size")
            length = len(x)
            xy = numpy.hstack((x.reshape((length, 1)),
                               y.reshape((length, 1))))
            sky = func(xy)
            return [sky[:, i] for i in range(sky.shape[1])]
        raise TypeError("Expected 1 or 2 arguments, %d given" % len(args))

    def all_pix2sky(self, *args):
        return self._array_converter(self._all_pix2sky, *args)
    all_pix2sky.__doc__ = """
        all_pix2sky(*args) -> sky

        Transforms pixel coordinates to sky coordinates by doing all
        of the following in order:

            - SIP distortion correction (optionally)

            - Paper IV table-lookup distortion correction (optionally)

            - wcslib WCS transformation

        %s

        %s

        @raises MemoryError: Memory allocation failed.
        @raises SingularMatrixError: Linear transformation matrix is
            singular.
        @raises InconsistentAxisTypesError: Inconsistent or
            unrecognized coordinate axis types.
        @raises ValueError: Invalid parameter value.
        @raises ValueError: Invalid coordinate transformation
            parameters.
        @raises ValueError: x- and y-coordinate arrays are not the
            same size.
        @raises InvalidTransformError: Invalid coordinate
            transformation parameters.
        @raises InvalidTransformError: Ill-conditioned coordinate
            transformation parameters.
        """ % (__.FITS_EQUIVALENT('all_pix2sky', 8),
               __.ONE_OR_TWO_ARGS('pixel', 8))

    def all_pix2sky_fits(self, *args):
        return self._array_converter(self._all_pix2sky_fits, *args)
    all_pix2sky_fits.__doc__ = """
        all_pix2sky_fits(*args) -> sky

        %s
        """ % (__.NON_FITS_EQUIVALENT('all_pix2sky', 8))

    def wcs_pix2sky(self, *args):
        if self.wcs is None:
            raise ValueError("No basic WCS settings were created.")
        return self._array_converter(lambda x: self.wcs.p2s(x)['world'], *args)
    wcs_pix2sky.__doc__ = """
        wcs_pix2sky(*args) -> sky

        Transforms pixel coordinates to sky coordinates by doing only
        the basic wcslib transformation.  No SIP or Paper IV table
        lookup distortion correction is applied.  To perform
        distortion correction, see L{all_pix2sky}, L{sip_pix2foc},
        L{p4_pix2foc}, and L{pix2foc}.

        %s

        %s

        @raises MemoryError: Memory allocation failed.
        @raises SingularMatrixError: Linear transformation matrix is
            singular.
        @raises InconsistentAxisTypesError: Inconsistent or
            unrecognized coordinate axis types.
        @raises ValueError: Invalid parameter value.
        @raises ValueError: Invalid coordinate transformation
            parameters.
        @raises ValueError: x- and y-coordinate arrays are not the
            same size.
        @raises InvalidTransformError: Invalid coordinate
            transformation parameters.
        @raises InvalidTransformError: Ill-conditioned coordinate
            transformation parameters.
        """ % (__.FITS_EQUIVALENT('wcs_pix2sky', 8),
               __.ONE_OR_TWO_ARGS('sky', 8))

    def wcs_pix2sky_fits(self, *args):
        return self._array_converter(lambda x: self.wcs.p2s(x)['world'], *args)
    wcs_pix2sky_fits.__doc__ = """
        wcs_pix2sky_fits(*args) -> sky

        %s
        """ % (__.NON_FITS_EQUIVALENT('wcs_pix2sky', 8))

    def wcs_sky2pix(self, *args):
        if self.wcs is None:
            raise ValueError("No basic WCS settings were created.")
        return self._array_converter(lambda x: self.wcs.s2p(x)['pixel'], *args)
    wcs_sky2pix.__doc__ = """
        wcs_sky2pix(*args) -> pixel

        Transforms sky coordinates to pixel coordinates, using only
        the basic libwcs WCS transformation.  No SIP or Paper IV table
        lookup distortion is applied.

        %s

        %s

        @raises MemoryError: Memory allocation failed.
        @raises SingularMatrixError: Linear transformation matrix is
            singular.
        @raises InconsistentAxisTypesError: Inconsistent or
            unrecognized coordinate axis types.
        @raises ValueError: Invalid parameter value.
        @raises InvalidTransformError: Invalid coordinate
            transformation parameters.
        @raises InvalidTransformError: Ill-conditioned coordinate
            transformation parameters.
        """ % (__.FITS_EQUIVALENT('wcs_sky2pix', 8),
               __.ONE_OR_TWO_ARGS('pixel', 8))

    def wcs_sky2pix_fits(self, *args):
        if self.wcs is None:
            raise ValueError("No basic WCS settings were created.")
        return self._array_converter(lambda x: self.wcs.s2p_fits(x)['pixel'],
                                     *args)
    wcs_sky2pix_fits.__doc__ = """
        wcs_sky2pix_fits(*args) -> pixel

        %s
        """ % (__.NON_FITS_EQUIVALENT('wcs_sky2pix', 8))

    def pix2foc(self, *args):
        return self._array_converter(self._pix2foc, *args)
    pix2foc.__doc__ = """
        pix2foc(*args) -> focal plane

        Convert pixel coordinates to focal plane coordinates using the
        SIP polynomial distortion convention and Paper IV table-lookup
        distortion correction.

        %s

        %s

        @raises MemoryError: Memory allocation failed.
        @raises ValueError: Invalid coordinate transformation parameters.
        """ % (__.FITS_EQUIVALENT('pix2foc', 8),
               __.ONE_OR_TWO_ARGS('pixel', 8))

    def pix2foc_fits(self, *args):
        return self._array_converter(self._pix2foc_fits, *args)
    pix2foc_fits.__doc__ = """
        pix2foc_fits(*args) -> focal plane

        %s
        """ % (__.NON_FITS_EQUIVALENT('pix2foc', 8))

    def p4_pix2foc(self, *args):
        return self._array_converter(self._p4_pix2foc, *args)
    p4_pix2foc.__doc__ = """
        p4_pix2foc(*args) -> focal plane

        Convert pixel coordinates to focal plane coordinates using
        Paper IV table-lookup distortion correction.

        %s

        %s

        @raises MemoryError: Memory allocation failed.
        @raises ValueError: Invalid coordinate transformation parameters.
        """ % (__.FITS_EQUIVALENT('p4_pix2foc', 8),
               __.ONE_OR_TWO_ARGS('pixel', 8))

    def p4_pix2foc_fits(self, *args):
        return self._array_converter(self._p4_pix2foc_fits, *args)
    p4_pix2foc_fits.__doc__ = """
        p4_pix2foc_fits(*args) -> focal plane

        %s
        """ % (__.NON_FITS_EQUIVALENT('p4_pix2foc', 8))

    def sip_pix2foc(self, *args):
        return self._array_converter(self.sip.pix2foc, *args)
    sip_pix2foc.__doc__ = """
        sip_pix2foc(*args) -> focal plane

        Convert pixel coordinates to focal plane coordinates using the
        SIP polynomial distortion convention.

        Paper IV table lookup distortion correction is not applied,
        even if that information existed in the FITS file that
        initialized this L{WCS} object.  To correct for that use
        L{pix2foc} or L{p4_pix2foc}.

        %s

        %s

        @raises MemoryError: Memory allocation failed.
        @raises ValueError: Invalid coordinate transformation parameters.
        """ % (__.FITS_EQUIVALENT('sip_pix2foc', 8),
               __.ONE_OR_TWO_ARGS('pixel', 8))

    def sip_pix2foc_fits(self, *args):
        return self._array_converter(self.sip.pix2foc_fits, *args)
    sip_pix2foc_fits.__doc__ = """
        sip_pix2foc_fits(*args) -> focal plane

        %s
        """ % (__.NON_FITS_EQUIVALENT('sip_pix2foc', 8))

    def sip_foc2pix(self, *args):
        return self._array_converter(self.sip.foc2pix, *args)
    sip_foc2pix.__doc__ = """
        sip_foc2pix(*args) -> pixel

        Convert focal plane coordinates to pixel coordinates using the
        SIP polynomial distortion convention.

        Paper IV table lookup distortion correction is not applied,
        even if that information existed in the FITS file that
        initialized this L{WCS} object.

        %s

        %s

        @raises MemoryError: Memory allocation failed.
        @raises ValueError: Invalid coordinate transformation parameters.
        """ % (__.FITS_EQUIVALENT('sip_foc2pix', 8),
               __.ONE_OR_TWO_ARGS('focal plane', 8))

    def sip_foc2pix_fits(self, *args):
        return self._array_converter(self.sip.foc2pix_fits, *args)
    sip_foc2pix_fits.__doc__ = """
        sip_foc2pix_fits(*args) -> pixels

        %s
        """ % (__.NON_FITS_EQUIVALENT('sip_foc2pix', 8))

    def to_header(self, relax=False):
        """
        Generate a PyFITS header object with the WCS information stored in
        this object.

        B{This function does not write out SIP or Paper IV distortion
        keywords, yet, only the core WCS support by wcslib.}

        The output header will almost certainly differ from the input in a
        number of respects:

          1. The output header only contains WCS-related keywords.  In
             particular, it does not contain syntactically-required
             keywords such as C{SIMPLE}, C{NAXIS}, C{BITPIX}, or C{END}.

          2. Deprecated (e.g. C{CROTAn}) or non-standard usage will be
             translated to standard (this is partially dependent on
             whether C{fix} was applied).

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
        header_string = self.wcs.to_header(relax)
        cards = pyfits.CardList()
        for i in range(0, len(header_string), 80):
            card_string = header_string[i:i+80]
            card = pyfits.Card()
            card.fromstring(card_string)
            cards.append(card)
        return pyfits.Header(cards)

    def to_header_string(self, relax=False):
        return self.to_header(self, relax).to_string()
    to_header_string.__doc__ = to_header.__doc__
