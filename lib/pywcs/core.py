# Copyright (C) 2008-2012 Association of Universities for Research in
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
Under the hood, there are 3 separate classes that perform different
parts of the transformation:

   - `~pywcs.Wcsprm`: Is a direct wrapper of the core WCS
     functionality in `wcslib`_.

   - `~pywcs.Sip`: Handles polynomial distortion as defined in the
     `SIP`_ convention.

   - `~pywcs.DistortionLookupTable`: Handles `Paper IV`_ distortion
     lookup tables.

Additionally, the class `WCS` aggregates all of these transformations
together in a pipeline:

   - Detector to image plane correction (by a pair of
     `~pywcs.DistortionLookupTable` objects).

   - `SIP`_ distortion correction (by an underlying `~pywcs.Sip`
     object)

   - `Paper IV`_ table-lookup distortion correction (by a pair of
     `~pywcs.DistortionLookupTable` objects).

   - `wcslib`_ WCS transformation (by a `~pywcs.Wcsprm` object)
"""

from __future__ import division # confidence high

# stdlib
import os
import sys
import copy
import warnings
if sys.version_info[0:2] >= (2, 6):
    from io import BytesIO
else:
    from cStringIO import StringIO as BytesIO

# third-party
import numpy as np
try:
    import pyfits
    HAS_PYFITS = True
except ImportError:
    HAS_PYFITS = False

# local
if sys.version_info[0] >= 3:
    from . import _docutil as __
    from . import _pywcs
else:
    import _docutil as __
    import _pywcs

assert _pywcs._sanity_check(), \
    """PyWcs did not pass its sanity check for your build on your platform.
Please send details about your build and platform to mdroe@stsci.edu"""

if sys.version_info[0] >= 3:
    string_types = (bytes,)
else:
    string_types = (str, unicode)

# bit flags for relax keyword
WCSHDO_SIP =  0x40

# This is here for the sake of epydoc
WCSBase = _pywcs._Wcs
DistortionLookupTable = _pywcs.DistortionLookupTable
Sip = _pywcs.Sip
UnitConverter = _pywcs.UnitConverter
class Wcsprm(_pywcs._Wcsprm): pass
# Copy all the constants from the C extension into this module's namespace
for key, val in _pywcs.__dict__.items():
    if (key.startswith('WCSSUB') or
        key.startswith('WCSHDR') or
        key.startswith('WCSHDO')):
        locals()[key] = val

# A wrapper around the C WCS type

def _parse_keysel(keysel):
    keysel_flags = 0
    if keysel is not None:
        for element in keysel:
            if element.lower() == 'image':
                keysel_flags |= _pywcs.WCSHDR_IMGHEAD
            elif element.lower() == 'binary':
                keysel_flags |= _pywcs.WCSHDR_BIMGARR
            elif element.lower() == 'pixel':
                keysel_flags |= _pywcs.WCSHDR_PIXLIST
            else:
                raise ValueError(
                    "keysel must be a list of 'image', 'binary' and/or 'pixel'")
    else:
        keysel_flags = -1

    return keysel_flags


class WCS(WCSBase):
    """
    WCS objects perform standard WCS transformations, and correct for
    `SIP`_ and `Paper IV`_ table-lookup distortions, based on the WCS
    keywords and supplementary data read from a FITS file.
    """

    def __init__(self, header=None, fobj=None, key=' ', minerr=0.0,
                 relax=False, naxis=None, keysel=None, colsel=None,
                 fix=True):
        """
        - *header*: A string containing the header content, or a
          PyFITS header object.  If *header* is not provided or None,
          the object will be initialized to default values.

        - *fobj*: A PyFITS file (hdulist) object. It is needed when
          header keywords point to a `Paper IV`_ Lookup table
          distortion stored in a different extension.

        - *key*: A string.  The name of a particular WCS transform to
          use.  This may be either ``' '`` or ``'A'``-``'Z'`` and
          corresponds to the ``"a"`` part of the ``CTYPEia`` cards.
          *key* may only be provided if *header* is also provided.

        - *minerr*: A floating-point value.  The minimum value a
          distortion correction must have in order to be applied. If
          the value of ``CQERRja`` is smaller than *minerr*, the
          corresponding distortion is not applied.

        - *relax*: Degree of permissiveness:

            - `False`: Recognize only FITS keywords defined by the
              published WCS standard.

            - `True`: Admit all recognized informal extensions of the
              WCS standard.

            - `int`: a bit field selecting specific extensions to
              accept.  See :ref:`relaxread` for details.

        - *naxis*: int or sequence.  Extracts specific coordinate axes
          using :meth:`~pywcs.Wcsprm.sub`.  If a header is provided,
          and *naxis* is not ``None``, *naxis* will be passed to
          :meth:`~pywcs.Wcsprm.sub` in order to select specific axes
          from the header.  See :meth:`~pywcs.Wcsprm.sub` for more
          details about this parameter.

        - *keysel*: A list of flags used to select the keyword types
          considered by wcslib.  When ``None``, only the standard
          image header keywords are considered (and the underlying
          wcspih() C function is called).  To use binary table image
          array or pixel list keywords, *keysel* must be set.

          Each element in the list should be one of the following
          strings:

            - 'image': Image header keywords

            - 'binary': Binary table image array keywords

            - 'pixel': Pixel list keywords

          Keywords such as ``EQUIna`` or ``RFRQna`` that are common to
          binary table image arrays and pixel lists (including
          ``WCSNna`` and ``TWCSna``) are selected by both 'binary' and
          'pixel'.

        - *colsel*: A sequence of table column numbers used
          to restrict the WCS transformations considered to only those
          pertaining to the specified columns.  If `None`, there is no
          restriction.

        - *fix* : When `True` (default), call `fix` on the resulting
          object to fix any non-standard uses in the header.

        .. warning::

          pywcs supports arbitrary *n* dimensions for the core WCS
          (the transformations handled by WCSLIB).  However, the Paper
          IV lookup table and SIP distortions must be two dimensional.
          Therefore, if you try to create a WCS object where the core
          WCS has a different number of dimensions than 2 and that
          object also contains a Paper IV lookup table or SIP
          distortion, a `ValueError` exception will be raised.  To
          avoid this, consider using the *naxis* kwarg to select two
          dimensions from the core WCS.

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `ValueError`: Invalid key.

        - `KeyError`: Key not found in FITS header.

        - `AssertionError`: Lookup table distortion present in the
          header but fobj not provided.
        """
        if header is None:
            if naxis is None:
                naxis = 2
            wcsprm = _pywcs._Wcsprm(header=None, key=key,
                                    relax=relax, naxis=naxis)
            self.naxis = wcsprm.naxis
            # Set some reasonable defaults.
            det2im = (None, None)
            cpdis = (None, None)
            sip = None
        else:
            keysel_flags = _parse_keysel(keysel)

            if isinstance(header, string_types):
                if HAS_PYFITS and os.path.exists(header):
                    hdulist = pyfits.open(header)
                    header = hdulist[0].header
                    header_string = header.tostring()
                else:
                    header_string = header
            elif HAS_PYFITS:
                assert isinstance(header, pyfits.Header)
                header_string = header.tostring()
            else:
                raise TypeError(
                    "header must be a string or a pyfits.Header object")
            try:
                wcsprm = _pywcs._Wcsprm(header=header_string, key=key,
                                        relax=relax, keysel=keysel_flags,
                                        colsel=colsel)
            except _pywcs.NoWcsKeywordsFoundError:
                # The header may have SIP or distortions, but no core
                # WCS.  That isn't an error -- we want a "default"
                # (identity) core Wcs transformation in that case.
                if colsel is None:
                    wcsprm = _pywcs._Wcsprm(header=None, key=key,
                                            relax=relax, keysel=keysel_flags,
                                            colsel=colsel)
                else:
                    raise

            if naxis is not None:
                wcsprm = wcsprm.sub(naxis)
            self.naxis = wcsprm.naxis

            det2im = self._read_det2im_kw(header, fobj, err=minerr)
            cpdis = self._read_distortion_kw(
                header, fobj, dist='CPDIS', err=minerr)
            sip = self._read_sip_kw(header)
            if (wcsprm.naxis != 2 and
                (det2im[0] or det2im[1] or cpdis[0] or cpdis[1] or sip)):
                raise ValueError(
                    """
Paper IV lookup tables and SIP distortions only work in 2 dimensions.
However, WCSLIB has detected %d dimensions in the core WCS keywords.
To use core WCS in conjunction with Paper IV lookup tables or SIP
distortion, you must select or reduce these to 2 dimensions using the
naxis kwarg.
""" %
                    wcsprm.naxis)

        if fix:
            wcsprm.fix()

        self.get_naxis(header)
        WCSBase.__init__(self, sip, cpdis, wcsprm, det2im)

    def __copy__(self):
        new_copy = self.__class__()
        WCSBase.__init__(new_copy, self.sip,
                         (self.cpdis1, self.cpdis2),
                         self.wcs,
                         (self.det2im1, self.det2im2))
        new_copy.__dict__.update(self.__dict__)
        return new_copy

    def __deepcopy__(self, memo):
        new_copy = self.__class__()
        new_copy.naxis = copy.deepcopy(self.naxis, memo)
        WCSBase.__init__(new_copy, copy.deepcopy(self.sip, memo),
                         (copy.deepcopy(self.cpdis1, memo),
                          copy.deepcopy(self.cpdis2, memo)),
                         copy.deepcopy(self.wcs, memo),
                         (copy.deepcopy(self.det2im1, memo),
                          copy.deepcopy(self.det2im2, memo)))
        for key in self.__dict__:
            val = self.__dict__[key]
            new_copy.__dict__[key] = copy.deepcopy(val, memo)
        return new_copy

    def copy(self):
        """
        Return a shallow copy of the object.

        Convenience method so user doesn't have to import the :mod:`copy`
        stdlib module.
        """
        return copy.copy(self)

    def deepcopy(self):
        """
        Return a deep copy of the object.

        Convenience method so user doesn't have to import the :mod:`copy`
        stdlib module.
        """
        return copy.deepcopy(self)

    def sub(self, axes=None):
        copy = self.deepcopy()
        copy.wcs = self.wcs.sub(axes)
        copy.naxis = copy.wcs.naxis
        return copy
    sub.__doc__ = _pywcs._Wcsprm.sub.__doc__

    def calcFootprint(self, header=None, undistort=True):
        """
        Calculates the footprint of the image on the sky.

        A footprint is defined as the positions of the corners of the
        image on the sky after all available distortions have been
        applied.

        Returns a (4, 2) array of (*x*, *y*) coordinates.
        """
        if header is None:
            try:
                # classes that inherit from WCS and define naxis1/2
                # do not require a header parameter
                naxis1 = self.naxis1
                naxis2 = self.naxis2
            except AttributeError :
                print("Need a valid header in order to calculate footprint\n")
                return None
        else:
            naxis1 = header.get('NAXIS1', None)
            naxis2 = header.get('NAXIS2', None)

        corners = np.zeros(shape=(4,2),dtype=np.float64)
        if naxis1 is None or naxis2 is None:
            return None

        corners[0,0] = 1.
        corners[0,1] = 1.
        corners[1,0] = 1.
        corners[1,1] = naxis2
        corners[2,0] = naxis1
        corners[2,1] = naxis2
        corners[3,0] = naxis1
        corners[3,1] = 1.
        if undistort:
            return self.all_pix2sky(corners, 1)
        else:
            return self.wcs_pix2sky(corners,1)

    def _read_det2im_kw(self, header, fobj, err=0.0):
        """
        Create a `Paper IV`_ type lookup table for detector to image
        plane correction.
        """

        
        if fobj is None:
            return (None, None)

        if not HAS_PYFITS:
            raise ImportError(
                "pyfits is required to use Paper IV lookup tables")

        if not isinstance(fobj, pyfits.HDUList):
            return (None, None)

        try:            
            axiscorr = header['AXISCORR']
            d2imdis = self._old_style_d2im(header, fobj, axiscorr)
            return d2imdis
        except KeyError:
            pass
        
        dist = 'D2IMDIS'
        d_kw = 'D2IM'
        err_kw = 'D2IMERR'
        tables = {}
        for i in range(1, self.naxis+1):
            d_error = header.get(err_kw+str(i), 0.0)
            if d_error < err:
                tables[i] = None
                continue
            distortion = dist+str(i)
            if distortion in header:
                dis = header[distortion].lower()
                if dis == 'lookup':
                    if fobj is not None and not HAS_PYFITS:
                        raise ImportError(
                            "pyfits is required to use Paper IV lookup tables")

                    assert isinstance(fobj, pyfits.HDUList), \
                        'A pyfits HDUList is required for Lookup table distortion.'
                    dp = (d_kw+str(i)).strip()
                    d_extver = header.get(dp+'.EXTVER', 1)
                    if i == header[dp+'.AXIS.%s'%i]:
                        d_data = fobj['D2IMARR', d_extver].data
                    else:
                        d_data = (fobj['D2IMARR', d_extver].data).transpose()
                    d_header = fobj['D2IMARR', d_extver].header
                    d_crpix = (d_header.get('CRPIX1', 0.0), d_header.get('CRPIX2', 0.0))
                    d_crval = (d_header.get('CRVAL1', 0.0), d_header.get('CRVAL2', 0.0))
                    d_cdelt = (d_header.get('CDELT1', 1.0), d_header.get('CDELT2', 1.0))
                    d_lookup = DistortionLookupTable(d_data, d_crpix,
                                                     d_crval, d_cdelt)
                    tables[i] = d_lookup
                else:
                    print('Polynomial distortion is not implemented.\n')
            else:
                tables[i] = None
        if not tables:
            return (None, None)
        else:
            return (tables.get(1), tables.get(2))
        
    def _old_style_d2im(self, header, fobj, axiscorr):
        warnings.warn("The use of ``AXISCORR`` for D2IM correction has been deprecated."
                      "The new style of this correction is described at"
                      ""
                      "PyWCS will read in files with ``AXISCORR`` but to_fits() will write"
                      "out files in the new style",
                      DeprecationWarning)
        cpdis = [None, None]
        crpix = [0.,0.]
        crval = [0.,0.]
        cdelt = [1.,1.]
        try:
            d2im_data = fobj[('D2IMARR', 1)].data
        except KeyError:
            return (None, None)
        except AttributeError:
            return (None, None)
        
        d2im_data = np.array([d2im_data])
        d2im_hdr = fobj[('D2IMARR', 1)].header
        naxis = d2im_hdr['NAXIS']

        for i in range(1,naxis+1):
            crpix[i-1] = d2im_hdr.get('CRPIX'+str(i), 0.0)
            crval[i-1] = d2im_hdr.get('CRVAL'+str(i), 0.0)
            cdelt[i-1] = d2im_hdr.get('CDELT'+str(i), 1.0)

        cpdis = DistortionLookupTable(d2im_data, crpix, crval, cdelt)

        if axiscorr == 1:
            return (cpdis, None)
        elif axiscorr == 2:
            return (None, cpdis)
        else:
            print "Expected AXISCORR to be 1 or 2"
            return (None, None)

    def _write_det2im(self, hdulist):
        """
        Write out Detector to Image Correction to the given
        `pyfits.HDUList`.
        """
        if self.det2im1 is None and self.det2im2 is None:
            return
        
        dist = 'D2IMDIS'
        d_kw = 'D2IM'
        err_kw = 'D2IMERR'
            
        def write_d2i(num, det2im):
            if det2im is None:
                return

            hdulist[0].header.update('%s%d' % (dist, num), value='LOOKUP',
                                     comment='Detector to image correction type')
            hdulist[0].header.update('%s%d.EXTVER' % (d_kw, num), value=num,
                                     comment='Version number of WCSDVARR extension')
            hdulist[0].header.update('%s%d.NAXES' % (d_kw, num), value=len(det2im.data.shape),
                            comment='Number of independent variables in d2im function')
            for i in range(det2im.data.ndim):
                hdulist[0].header.update('%s%d.AXIS.%d' % (d_kw, num, i+1), value=i+1,
                                comment='Axis number of the jth independent variable in a d2im function')

            image = pyfits.ImageHDU(det2im.data, name='D2IMARR')
            header = image.header

            header.update('CRPIX1', value=det2im.crpix[0],
                                comment='Coordinate system reference pixel')
            header.update('CRPIX2', value=det2im.crpix[1],
                                comment='Coordinate system reference pixel')
            header.update('CRVAL1', value=det2im.crval[0],
                                comment='Coordinate system value at reference pixel')
            header.update('CRVAL2', value=det2im.crval[1],
                                comment='Coordinate system value at reference pixel')
            header.update('CDELT1', value=det2im.cdelt[0],
                                comment='Coordinate increment along axis')
            header.update('CDELT2', value=det2im.cdelt[1],
                                comment='Coordinate increment along axis')
            image.update_ext_version(int(hdulist[0].header['%s%d.EXTVER' % (d_kw, num)]))
            hdulist.append(image)

        write_d2i(1, self.det2im1)
        write_d2i(2, self.det2im2)

    def _read_distortion_kw(self, header, fobj, dist='CPDIS', err=0.0):
        """
        Reads `Paper IV`_ table-lookup distortion keywords and data,
        and returns a 2-tuple of `~pywcs.DistortionLookupTable`
        objects.

        If no `Paper IV`_ distortion keywords are found, ``(None,
        None)`` is returned.
        """
        if isinstance(header, string_types):
            return (None, None)

        if dist == 'CPDIS':
            d_kw = 'DP'
            err_kw = 'CPERR'
        else:
            d_kw = 'DQ'
            err_kw = 'CQERR'

        tables = {}
        for i in range(1, self.naxis+1):
            d_error = header.get(err_kw+str(i), 0.0)
            if d_error < err:
                tables[i] = None
                continue
            distortion = dist+str(i)
            if distortion in header:
                dis = header[distortion].lower()
                if dis == 'lookup':
                    if fobj is not None and not HAS_PYFITS:
                        raise ImportError(
                            "pyfits is required to use Paper IV lookup tables")

                    assert isinstance(fobj, pyfits.HDUList), \
                        'A pyfits HDUList is required for Lookup table distortion.'
                    dp = (d_kw+str(i)).strip()
                    d_extver = header.get(dp+'.EXTVER', 1)
                    if i == header[dp+'.AXIS.%s'%i]:
                        d_data = fobj['WCSDVARR', d_extver].data
                    else:
                        d_data = (fobj['WCSDVARR', d_extver].data).transpose()
                    d_header = fobj['WCSDVARR', d_extver].header
                    d_crpix = (d_header.get('CRPIX1', 0.0), d_header.get('CRPIX2', 0.0))
                    d_crval = (d_header.get('CRVAL1', 0.0), d_header.get('CRVAL2', 0.0))
                    d_cdelt = (d_header.get('CDELT1', 1.0), d_header.get('CDELT2', 1.0))
                    d_lookup = DistortionLookupTable(d_data, d_crpix,
                                                     d_crval, d_cdelt)
                    tables[i] = d_lookup
                else:
                    print('Polynomial distortion is not implemented.\n')
            else:
                tables[i] = None

        if not tables:
            return (None, None)
        else:
            return (tables.get(1), tables.get(2))

    def _write_distortion_kw(self, hdulist, dist='CPDIS'):
        """
        Write out Paper IV distortion keywords to the given
        `pyfits.HDUList`.
        """
        if self.cpdis1 is None and self.cpdis2 is None:
            return

        if dist == 'CPDIS':
            d_kw = 'DP'
            err_kw = 'CPERR'
        else:
            d_kw = 'DQ'
            err_kw = 'CQERR'

        def write_dist(num, cpdis):
            if cpdis is None:
                return

            hdulist[0].header.update('%s%d' % (dist, num), value='LOOKUP',
                                     comment='Prior distortion function type')
            hdulist[0].header.update('%s%d.EXTVER' % (d_kw, num), value=num,
                                     comment='Version number of WCSDVARR extension')
            hdulist[0].header.update('%s%d.NAXES' % (d_kw, num), value=len(cpdis.data.shape),
                            comment='Number of independent variables in distortion function')
            for i in range(cpdis.data.ndim):
                hdulist[0].header.update('%s%d.AXIS.%d' % (d_kw, num, i+1), value=i+1,
                                comment='Axis number of the jth independent variable in a distortion function')

            image = pyfits.ImageHDU(cpdis.data, name='WCSDVARR')
            header = image.header

            header.update('CRPIX1', value=cpdis.crpix[0],
                                comment='Coordinate system reference pixel')
            header.update('CRPIX2', value=cpdis.crpix[1],
                                comment='Coordinate system reference pixel')
            header.update('CRVAL1', value=cpdis.crval[0],
                                comment='Coordinate system value at reference pixel')
            header.update('CRVAL2', value=cpdis.crval[1],
                                comment='Coordinate system value at reference pixel')
            header.update('CDELT1', value=cpdis.cdelt[0],
                                comment='Coordinate increment along axis')
            header.update('CDELT2', value=cpdis.cdelt[1],
                                comment='Coordinate increment along axis')
            image.update_ext_version(int(hdulist[0].header['%s%d.EXTVER' % (d_kw, num)]))
            hdulist.append(image)

        write_dist(1, self.cpdis1)
        write_dist(2, self.cpdis2)

    def _read_sip_kw(self, header):
        """
        Reads `SIP`_ header keywords and returns a `~pywcs.Sip`
        object.

        If no `SIP`_ header keywords are found, ``None`` is returned.
        """
        if isinstance(header, string_types):
            # TODO: Parse SIP from a string without pyfits around
            return None

        if "A_ORDER" in header and header['A_ORDER'] > 1:
            if "B_ORDER" not in header:
                raise ValueError(
                    "A_ORDER provided without corresponding B_ORDER "
                    "keyword for SIP distortion")

            m = int(header["A_ORDER"])
            a = np.zeros((m+1, m+1), np.double)
            for i in range(m+1):
                for j in range(m-i+1):
                    a[i, j] = header.get(("A_%d_%d" % (i, j)), 0.0)

            m = int(header["B_ORDER"])
            if m > 1:
                b = np.zeros((m+1, m+1), np.double)
                for i in range(m+1):
                    for j in range(m-i+1):
                        b[i, j] = header.get(("B_%d_%d" % (i, j)), 0.0)
            else:
                a = None
                b = None
        elif "B_ORDER" in header and header['B_ORDER'] > 1:
            raise ValueError(
                "B_ORDER provided without corresponding A_ORDER "
                "keyword for SIP distortion")
        else:
            a = None
            b = None

        if "AP_ORDER" in header and header['AP_ORDER'] > 1:
            if "BP_ORDER" not in header:
                raise ValueError(
                    "AP_ORDER provided without corresponding BP_ORDER "
                    "keyword for SIP distortion")

            m = int(header["AP_ORDER"])
            ap = np.zeros((m+1, m+1), np.double)
            for i in range(m+1):
                for j in range(m-i+1):
                    ap[i, j] = header.get("AP_%d_%d" % (i, j), 0.0)

            m = int(header["BP_ORDER"])
            if m > 1:
                bp = np.zeros((m+1, m+1), np.double)
                for i in range(m+1):
                    for j in range(m-i+1):
                        bp[i, j] = header.get("BP_%d_%d" % (i, j), 0.0)
            else:
                ap = None
                bp = None
        elif "BP_ORDER" in header and header['BP_ORDER'] > 1:
            raise ValueError(
                "BP_ORDER provided without corresponding AP_ORDER "
                "keyword for SIP distortion")
        else:
            ap = None
            bp = None

        if a is None and b is None and ap is None and bp is None:
            return None

        if "CRPIX1" not in header or "CRPIX2" not in header:
            raise ValueError(
                "Header has SIP keywords without CRPIX keywords")

        crpix1 = header.get("CRPIX1")
        crpix2 = header.get("CRPIX2")

        return Sip(a, b, ap, bp, (crpix1, crpix2))

    def _write_sip_kw(self):
        """
        Write out SIP keywords.  Returns a dictionary of key-value
        pairs.
        """
        if self.sip is None:
            return {}

        keywords = {}

        def write_array(name, a):
            if a is None:
                return
            size = a.shape[0]
            keywords['%s_ORDER' % name] = size - 1
            for i in range(size):
                for j in range(size - i):
                    if a[i, j] != 0.0:
                        keywords['%s_%d_%d' % (name, i, j)] = a[i, j]

        write_array('A', self.sip.a)
        write_array('B', self.sip.b)
        write_array('AP', self.sip.ap)
        write_array('BP', self.sip.bp)

        return keywords

    def _denormalize_sky(self, sky):
        if self.wcs.lngtyp != 'RA':
            raise ValueError(
                "WCS does not have longitude type of 'RA', therefore " +
                "(ra, dec) data can not be used as input")
        if self.wcs.lattype != 'DEC':
            raise ValueError(
                "WCS does not have longitude type of 'DEC', therefore " +
                "(ra, dec) data can not be used as input")
        if self.wcs.naxis == 2:
            if self.wcs.lng == 0 and self.wcs.lat == 1:
                return sky
            elif self.wcs.lng == 1 and self.wcs.lat == 0:
                # Reverse the order of the columns
                return sky[:,::-1]
            else:
                raise ValueError(
                    "WCS does not have longitude and latitude celestial " +
                    "axes, therefore (ra, dec) data can not be used as input")
        else:
            if self.wcs.lng < 0 or self.wcs.lat < 0:
                raise ValueError(
                    "WCS does not have both longitude and latitude celestial " +
                    "axes, therefore (ra, dec) data can not be used as input")
            out = np.zeros((sky.shape[0], self.wcs.naxis))
            out[:,self.wcs.lng] = sky[:,0]
            out[:,self.wcs.lat] = sky[:,1]
            return out

    def _normalize_sky(self, sky):
        if self.wcs.lngtyp != 'RA':
            raise ValueError(
                "WCS does not have longitude type of 'RA', therefore " +
                "(ra, dec) data can not be returned")
        if self.wcs.lattype != 'DEC':
            raise ValueError(
                "WCS does not have longitude type of 'DEC', therefore " +
                "(ra, dec) data can not be returned")
        if self.wcs.naxis == 2:
            if self.wcs.lng == 0 and self.wcs.lat == 1:
                return sky
            elif self.wcs.lng == 1 and self.wcs.lat == 0:
                # Reverse the order of the columns
                return sky[:,::-1]
            else:
                raise ValueError(
                    "WCS does not have longitude and latitude celestial "
                    "axes, therefore (ra, dec) data can not be returned")
        else:
            if self.wcs.lng < 0 or self.wcs.lat < 0:
                raise ValueError(
                    "WCS does not have both longitude and latitude celestial "
                    "axes, therefore (ra, dec) data can not be returned")
            out = np.empty((sky.shape[0], 2))
            out[:,0] = sky[:,self.wcs.lng]
            out[:,1] = sky[:,self.wcs.lat]
            return out

    def _array_converter(self, func, sky, *args, **kwargs):
        """
        A helper function to support reading either a pair of arrays
        or a single Nx2 array.
        """
        ra_dec_order = kwargs.get('ra_dec_order')
        if len(args) == 2:
            xy, origin = args
            try:
                xy = np.asarray(xy)
                origin = int(origin)
            except:
                raise TypeError(
                    "When providing two arguments, they must be (xy, origin)")
            if ra_dec_order and sky == 'input':
                xy = self._denormalize_sky(xy)
            result = func(xy, origin)
            if ra_dec_order and sky == 'output':
                result = self._normalize_sky(result)
            return result
        elif len(args) == 3:
            x, y, origin = args
            try:
                x = np.asarray(x)
                y = np.asarray(y)
                origin = int(origin)
            except:
                raise TypeError(
                    "When providing three arguments, they must be (x, y, origin)")
            if x.size != y.size:
                raise ValueError("x and y arrays are not the same size")
            length = x.size
            xy = np.hstack((x.reshape((length, 1)),
                            y.reshape((length, 1))))
            if ra_dec_order and sky == 'input':
                xy = self._denormalize_sky(xy)
            sky = func(xy, origin)
            if ra_dec_order and sky == 'output':
                sky = self._normalize_sky_output(sky)
                return sky[:, 0], sky[:, 1]
            return [sky[:, i] for i in range(sky.shape[1])]
        raise TypeError("Expected 2 or 3 arguments, %d given" % len(args))

    def all_pix2sky(self, *args, **kwargs):
        return self._array_converter(self._all_pix2sky, 'output', *args, **kwargs)
    all_pix2sky.__doc__ = """
        Transforms pixel coordinates to sky coordinates by doing all
        of the following in order:

            - Detector to image plane correction (optionally)

            - `SIP`_ distortion correction (optionally)

            - `Paper IV`_ table-lookup distortion correction (optionally)

            - `wcslib`_ WCS transformation

        %s

        %s

        For a transformation that is not two-dimensional, the
        two-argument form must be used.

        .. note::

            The order of the axes for the result is determined by the
            `CTYPEia` keywords in the FITS header, therefore it may
            not always be of the form (*ra*, *dec*).  The
            `~pywcs.Wcsprm.lat`, `~pywcs.Wcsprm.lng`,
            `~pywcs.Wcsprm.lattyp` and `~pywcs.Wcsprm.lngtyp` members
            can be used to determine the order of the axes.

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `SingularMatrixError`: Linear transformation matrix is
          singular.

        - `InconsistentAxisTypesError`: Inconsistent or unrecognized
          coordinate axis types.

        - `ValueError`: Invalid parameter value.

        - `ValueError`: Invalid coordinate transformation parameters.

        - `ValueError`: x- and y-coordinate arrays are not the same
          size.

        - `InvalidTransformError`: Invalid coordinate transformation
          parameters.

        - `InvalidTransformError`: Ill-conditioned coordinate
          transformation parameters.
        """ % (__.TWO_OR_THREE_ARGS(
            'sky coordinates, in degrees', 'naxis', 8),
               __.RA_DEC_ORDER(8))

    def wcs_pix2sky(self, *args, **kwargs):
        if self.wcs is None:
            raise ValueError("No basic WCS settings were created.")
        return self._array_converter(lambda xy, o: self.wcs.p2s(xy, o)['world'],
                                     'output', *args, **kwargs)
    wcs_pix2sky.__doc__ = """
        Transforms pixel coordinates to sky coordinates by doing only
        the basic `wcslib`_ transformation.  No `SIP`_ or `Paper IV`_
        table lookup distortion correction is applied.  To perform
        distortion correction, see `~pywcs.WCS.all_pix2sky`,
        `~pywcs.WCS.sip_pix2foc`, `~pywcs.WCS.p4_pix2foc`, or
        `~pywcs.WCS.pix2foc`.

        %s

        %s

        For a transformation that is not two-dimensional, the
        two-argument form must be used.

        .. note::

            The order of the axes for the result is determined by the
            `CTYPEia` keywords in the FITS header, therefore it may
            not always be of the form (*ra*, *dec*).  The
            `~pywcs.Wcsprm.lat`, `~pywcs.Wcsprm.lng`,
            `~pywcs.Wcsprm.lattyp` and `~pywcs.Wcsprm.lngtyp` members
            can be used to determine the order of the axes.

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `SingularMatrixError`: Linear transformation matrix is
          singular.

        - `InconsistentAxisTypesError`: Inconsistent or unrecognized
          coordinate axis types.

        - `ValueError`: Invalid parameter value.

        - `ValueError`: Invalid coordinate transformation parameters.

        - `ValueError`: x- and y-coordinate arrays are not the same
          size.

        - `InvalidTransformError`: Invalid coordinate transformation
          parameters.

        - `InvalidTransformError`: Ill-conditioned coordinate
          transformation parameters.
        """ % (__.TWO_OR_THREE_ARGS('sky coordinates, in degrees.', 'naxis', 8),
               __.RA_DEC_ORDER(8))

    def wcs_sky2pix(self, *args, **kwargs):
        if self.wcs is None:
            raise ValueError("No basic WCS settings were created.")
        return self._array_converter(lambda xy, o: self.wcs.s2p(xy, o)['pixcrd'],
                                     'input', *args, **kwargs)
    wcs_sky2pix.__doc__ = """
        Transforms sky coordinates to pixel coordinates, using only
        the basic `wcslib`_ WCS transformation.  No `SIP`_ or `Paper
        IV`_ table lookup distortion is applied.

        %s

        %s

        For a transformation that is not two-dimensional, the
        two-argument form must be used.

        .. note::

            The order of the axes for the input sky array is
            determined by the `CTYPEia` keywords in the FITS header,
            therefore it may not always be of the form (*ra*, *dec*).
            The `~pywcs.Wcsprm.lat`, `~pywcs.Wcsprm.lng`,
            `~pywcs.Wcsprm.lattyp` and `~pywcs.Wcsprm.lngtyp` members
            can be used to determine the order of the axes.

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `SingularMatrixError`: Linear transformation matrix is
          singular.

        - `InconsistentAxisTypesError`: Inconsistent or unrecognized
          coordinate axis types.

        - `ValueError`: Invalid parameter value.

        - `InvalidTransformError`: Invalid coordinate transformation
          parameters.

        - `InvalidTransformError`: Ill-conditioned coordinate
          transformation parameters.
        """ % (__.TWO_OR_THREE_ARGS('pixel coordinates', 'naxis', 8),
               __.RA_DEC_ORDER(8))

    def pix2foc(self, *args, **kwargs):
        return self._array_converter(self._pix2foc, None, *args, **kwargs)
    pix2foc.__doc__ = """
        Convert pixel coordinates to focal plane coordinates using the
        `SIP`_ polynomial distortion convention and `Paper IV`_
        table-lookup distortion correction.

        %s

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `ValueError`: Invalid coordinate transformation parameters.
        """ % (__.TWO_OR_THREE_ARGS('focal coordinates', '2', 8))

    def p4_pix2foc(self, *args, **kwargs):
        return self._array_converter(self._p4_pix2foc, None, *args, **kwargs)
    p4_pix2foc.__doc__ = """
        Convert pixel coordinates to focal plane coordinates using
        `Paper IV`_ table-lookup distortion correction.

        %s

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `ValueError`: Invalid coordinate transformation parameters.
        """ % (__.TWO_OR_THREE_ARGS('focal coordinates', '2', 8))

    def det2im(self, *args, **kwargs):
        return self._array_converter(self._det2im, None, *args, **kwargs)
    det2im.__doc__ = """
        Convert detector coordinates to image plane coordinates using
        `Paper IV`_ table-lookup distortion correction.

        %s

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `ValueError`: Invalid coordinate transformation parameters.
        """ % (__.TWO_OR_THREE_ARGS('pixel coordinates', '2', 8))

    def sip_pix2foc(self, *args, **kwargs):
        if self.sip is None:
            if len(args) == 2:
                return args[0]
            elif len(args) == 3:
                return args[:2]
            else:
                raise TypeError("Wrong number of arguments")
        return self._array_converter(self.sip.pix2foc, None, *args, **kwargs)
    sip_pix2foc.__doc__ = """
        Convert pixel coordinates to focal plane coordinates using the
        `SIP`_ polynomial distortion convention.

        `Paper IV`_ table lookup distortion correction is not applied,
        even if that information existed in the FITS file that
        initialized this :class:`~pywcs.WCS` object.  To correct for that,
        use `~pywcs.WCS.pix2foc` or `~pywcs.WCS.p4_pix2foc`.

        %s

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `ValueError`: Invalid coordinate transformation parameters.
        """ % (__.TWO_OR_THREE_ARGS('focal coordinates', '2', 8))

    def sip_foc2pix(self, *args, **kwargs):
        if self.sip is None:
            if len(args) == 2:
                return args[0]
            elif len(args) == 3:
                return args[:2]
            else:
                raise TypeError("Wrong number of arguments")
        return self._array_converter(self.sip.foc2pix, None, *args, **kwargs)
    sip_foc2pix.__doc__ = """
        Convert focal plane coordinates to pixel coordinates using the
        `SIP`_ polynomial distortion convention.

        `Paper IV`_ table lookup distortion correction is not applied,
        even if that information existed in the FITS file that
        initialized this `~pywcs.WCS` object.

        %s

        **Exceptions:**

        - `MemoryError`: Memory allocation failed.

        - `ValueError`: Invalid coordinate transformation parameters.
        """ % (__.TWO_OR_THREE_ARGS('pixel coordinates', '2', 8))

    def to_fits(self, relax=False, wkey=None):
        """
        Generate a `pyfits.HDUList` object with all of the information
        stored in this object.  This should be logically identical to
        the input FITS file, but it will be normalized in a number of
        ways.

        See `WCS.to_header` for some warnings about the output
        produced.

        - *relax*: Degree of permissiveness:

          - `False`: Recognize only FITS keywords defined by the
            published WCS standard.

          - `True`: Admit all recognized informal extensions of the WCS
            standard.

          - `int`: a bit field selecting specific extensions to write.
            See :ref:`relaxwrite` for details.

        Returns a `pyfits.HDUList` object.
        """
        if not HAS_PYFITS:
            raise ImportError(
                "pyfits is required to generate a FITS file")

        header = self.to_header(relax=relax, wkey=wkey)

        hdu = pyfits.PrimaryHDU(header=header)
        hdulist = pyfits.HDUList(hdu)

        self._write_det2im(hdulist)
        self._write_distortion_kw(hdulist)

        return hdulist

    def to_header(self, relax=False, wkey=None):
        """
        Generate a `pyfits.Header` object with the basic WCS and SIP
        information stored in this object.  This should be logically
        identical to the input FITS file, but it will be normalized in
        a number of ways.

        .. warning::

          This function does not write out Paper IV distortion
          information, since that requires multiple FITS header data
          units.  To get a full representation of everything in this
          object, use `to_fits`.

        The output header will almost certainly differ from the input in a
        number of respects:

          1. The output header only contains WCS-related keywords.  In
             particular, it does not contain syntactically-required
             keywords such as ``SIMPLE``, ``NAXIS``, ``BITPIX``, or
             ``END``.

          2. Deprecated (e.g. ``CROTAn``) or non-standard usage will
             be translated to standard (this is partially dependent on
             whether `fix` was applied).

          3. Quantities will be converted to the units used internally,
             basically SI with the addition of degrees.

          4. Floating-point quantities may be given to a different decimal
             precision.

          5. Elements of the ``PCi_j`` matrix will be written if and
             only if they differ from the unit matrix.  Thus, if the
             matrix is unity then no elements will be written.

          6. Additional keywords such as ``WCSAXES``, ``CUNITia``,
             ``LONPOLEa`` and ``LATPOLEa`` may appear.

          7. The original keycomments will be lost, although
             `to_header` tries hard to write meaningful comments.

          8. Keyword order may be changed.

        - *relax*: Degree of permissiveness:

          - `False`: Recognize only FITS keywords defined by the
            published WCS standard.

          - `True`: Admit all recognized informal extensions of the WCS
            standard.

          - `int`: a bit field selecting specific extensions to write.
            See :ref:`relaxwrite` for details.
            
        - *wkey*: A string.  The name of a particular WCS transform to
          use.  This may be either ``' '`` or ``'A'``-``'Z'`` and
          corresponds to the ``"a"`` part of the ``CTYPEia`` cards.

        Returns a `pyfits.Header` object.
        """
        if wkey:
            self.wcs.alt = wkey
        
        if not HAS_PYFITS:
            raise ImportError(
                "pyfits is required to generate a FITS header")

        dosip = False
        if (relax == True or relax == WCSHDO_all or
            (relax & WCSHDO_SIP)):
            dosip = True
        if relax not in (True, False):
            relax &= ~WCSHDO_SIP

        if self.wcs is not None:
            header_string = self.wcs.to_header(relax)
            cards = pyfits.CardList()
            for i in range(0, len(header_string), 80):
                card_string = header_string[i:i+80]
                if pyfits.__version__[0] >= '3':
                    card = pyfits.Card.fromstring(card_string)
                else:
                    card = pyfits.Card()
                    card.fromstring(card_string)
                
                cards.append(card)
            header = pyfits.Header(cards)
        else:
            header = pyfits.Header()

        if dosip and self.sip is not None:
            for key, val in self._write_sip_kw().items():
                header.update(key, val)

        return header

    def to_header_string(self, relax=False):
        """
        Identical to `to_header`, but returns a string containing the
        header cards.
        """
        return str(self.to_header(relax))

    def __reduce__(self):
        """
        Support pickling of WCS objects.  This is done by serializing
        to an in-memory FITS file and dumping that as a string.
        """
        if not HAS_PYFITS:
            raise ImportError(
                "pyfits is required for pickling support")

        hdulist = self.to_fits(relax=True)

        buffer = BytesIO()
        hdulist.writeto(buffer)

        return (_unpickle_wcs, (self.__class__, self.__dict__, buffer.getvalue(),))

    def footprint_to_file(self, filename=None, color='green', width=2):
        """
        Writes out a `ds9`_ style regions file. It can be loaded
        directly by `ds9`_.

        - *filename*: string.  Output file name - default is
          ``'footprint.reg'``

        - *color*: string.  Color to use when plotting the line.

        - *width*: int.  Width of the region line.
        """
        if not filename:
            filename = 'footprint.reg'
        comments = '# Region file format: DS9 version 4.0 \n'
        comments += '# global color=green font="helvetica 12 bold select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n'

        f = open(filename, 'a')
        f.write(comments)
        f.write('linear\n')
        f.write('polygon(')
        self.footprint.tofile(f, sep=',')
        f.write(') # color=%s, width=%d \n' % (color, width))
        f.close()

    def get_naxis(self, header=None):
        self.naxis1 = 0.0
        self.naxis2 = 0.0
        if header != None and not isinstance(header, string_types):
            self.naxis1 = header.get('NAXIS1', 0.0)
            self.naxis2 = header.get('NAXIS2', 0.0)

    def rotateCD(self, theta):
        _theta = DEGTORAD(theta)
        _mrot = np.zeros(shape=(2,2),dtype=np.double)
        _mrot[0] = (np.cos(_theta),np.sin(_theta))
        _mrot[1] = (-np.sin(_theta),np.cos(_theta))
        new_cd = np.dot(self.wcs.cd, _mrot)
        self.wcs.cd = new_cd

    def printwcs(self):
        """
        Temporary function for internal use.
        """
        print('WCS Keywords\n')
        if hasattr(self.wcs, 'cd'):
            print('CD_11  CD_12: %r %r' % (self.wcs.cd[0,0],  self.wcs.cd[0,1]))
            print('CD_21  CD_22: %r %r' % (self.wcs.cd[1,0],  self.wcs.cd[1,1]))
        print('CRVAL    : %r %r' % (self.wcs.crval[0], self.wcs.crval[1]))
        print('CRPIX    : %r %r' % (self.wcs.crpix[0], self.wcs.crpix[1]))
        print('NAXIS    : %r %r' % (self.naxis1, self.naxis2))

    def get_axis_types(self):
        """
        ``list of dicts``

        Similar to `self.wcsprm.axis_types <_pywcs._Wcsprm.axis_types>`
        but provides the information in a more Python-friendly format.

        Returns a list of dictionaries, one for each axis, each
        containing attributes about the type of that axis.

        Each dictionary has the following keys:

        - 'coordinate_type':

          - None: Non-specific coordinate type.

          - 'stokes': Stokes coordinate.

          - 'celestial': Celestial coordinate (including ``CUBEFACE``).

          - 'spectral': Spectral coordinate.

        - 'scale':

          - 'linear': Linear axis.

          - 'quantized': Quantized axis (``STOKES``, ``CUBEFACE``).

          - 'non-linear celestial': Non-linear celestial axis.

          - 'non-linear spectral': Non-linear spectral axis.

          - 'logarithmic': Logarithmic axis.

          - 'tabular': Tabular axis.

        - 'group'

          - Group number, e.g. lookup table number

        - 'number'

          - For celestial axes:

            - 0: Longitude coordinate.

            - 1: Latitude coordinate.

            - 2: ``CUBEFACE`` number.

          - For lookup tables:

            - the axis number in a multidimensional table.

        ``CTYPEia`` in ``"4-3"`` form with unrecognized algorithm code will
        generate an error.
        """
        if self.wcs is None:
            raise AttributeError(
                "This WCS object does not have a wcsprm object.")

        coordinate_type_map = {
            0: None,
            1: 'stokes',
            2: 'celestial',
            3: 'spectral'
            }

        scale_map = {
            0: 'linear',
            1: 'quantized',
            2: 'non-linear celestial',
            3: 'non-linear spectral',
            4: 'logarithmic',
            5: 'tabular'
            }

        result = []
        for axis_type in self.wcs.axis_types:
            subresult = {}

            coordinate_type = (axis_type // 1000) % 10
            subresult['coordinate_type'] = coordinate_type_map[coordinate_type]

            scale = (axis_type // 100) % 10
            subresult['scale'] = scale_map[scale]

            group = (axis_type // 10) % 10
            subresult['group'] = group

            number = axis_type % 10
            subresult['number'] = number

            result.append(subresult)

        return result


def DEGTORAD(deg):
    return (deg * np.pi / 180.)

def RADTODEG(rad):
    return (rad * 180. / np.pi)


def find_all_wcs(header, relax=False, keysel=None):
    """
    Find all the WCS transformations in the given header.

    - *header*: A string or PyFITS header object.

    - *relax*: Degree of permissiveness:

        - `False`: Recognize only FITS keywords defined by the
          published WCS standard.

        - `True`: Admit all recognized informal extensions of the
          WCS standard.

        - `int`: a bit field selecting specific extensions to accept.
          See :ref:`relaxread` for details.

    - *keysel*: A list of flags used to select the keyword types
      considered by wcslib.  When ``None``, only the standard image
      header keywords are considered (and the underlying wcspih() C
      function is called).  To use binary table image array or pixel
      list keywords, *keysel* must be set.

      Each element in the list should be one of the following strings:

        - 'image': Image header keywords

        - 'binary': Binary table image array keywords

        - 'pixel': Pixel list keywords

      Keywords such as ``EQUIna`` or ``RFRQna`` that are common to
      binary table image arrays and pixel lists (including ``WCSNna``
      and ``TWCSna``) are selected by both 'binary' and 'pixel'.

    Returns a list of `WCS` objects.
    """
    if isinstance(header, string_types):
        header_string = header
    elif HAS_PYFITS:
        assert isinstance(header, pyfits.Header)
        header_string = repr(header.ascard)
    else:
        raise TypeError(
            "header must be a string or pyfits.Header object")

    keysel_flags = _parse_keysel(keysel)

    wcsprms = _pywcs.find_all_wcs(header_string, relax, keysel_flags)

    result = []
    for wcsprm in wcsprms:
        subresult = WCS()
        subresult.wcs = wcsprm
        result.append(subresult)

    return result

def _unpickle_wcs(cls, dct, fits_data):
    """
    Unpickles a WCS object from a serialized FITS string.
    """
    if not HAS_PYFITS:
        raise ImportError(
            "pyfits is required for pickling support")

    self = cls.__new__(cls)
    self.__dict__.update(dct)

    buffer = BytesIO(fits_data)
    hdulist = pyfits.open(buffer)

    WCS.__init__(self, hdulist[0].header, hdulist)

    return self
