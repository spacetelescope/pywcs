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

from _pywcs import *
import _pywcs
try:
    import pyfits
    _has_pyfits = True
except ImportError:
    _has_pyfits = False

if _has_pyfits:
    def parse_hdulist(hdulist, relax=False):
        """parse_hdulist(hdulist, relax=False) -> list of L{Wcs} objects

Parses a FITS image header, either that of a primary HDU or of an image
extension.  All WCS keywords defined in Papers I, II, and III are
recognized, and also those used by the AIPS convention and certain
other keywords that existed in early drafts of the WCS papers.

Given a string containing a FITS image header, C{parse_image_header()}
identifies and reads all WCS keywords for the primary coordinate
representation and up to 26 alternate representations.  It returns
this information as a list of C{Wcs} objects.

C{parse_image_header()} fills in information associated with
coordinate lookup tables.

C{parse_hdulist} determines the number of coordinate axes independently for
each alternate coordinate representation (denoted by the C{"a"} value in
keywords like C{CTYPEia}) from the higher of
    - C{NAXIS}
    - C{WCSAXES}
    - The highest axis number in any parameterized WCS keyword.  The
      keyvalue, as well as the keyword, must be syntactically valid
      otherwise it will not be considered.

If none of these keyword types is present, i.e. if the header only
contains auxiliary WCS keywords for a particular coordinate
representation, then no coordinate description is constructed for it.

C{parse_hdulist} enforces correct FITS C{"keyword = value"} syntax
with regard to C{"= "} occurring in columns 9 and 10.  However, it
does recognize free-format character (NOST 100-2.0, Sect. 5.2.1),
integer (Sect. 5.2.3), and floating-point values (Sect. 5.2.4) for all
keywords.

@param hdulist: A PyFITS hdulist
@type hdulist: string
@param relax: Degree of permissiveness:
    - C{False}: Recognize only FITS keywords defined by the
      published WCS standard.
    - C{True}: Admit all recognized informal extensions of the
      WCS standard.
@type relax: bool
    """
        return parse_image_header(str(hdulist[0].header.ascardlist()), relax)

# This is a hack so epydoc will document this method
def parse_image_header(header, relax=False):
    return _pywcs.parse_image_header(header, relax)
parse_image_header.__doc__ = _pywcs.parse_image_header.__doc__
