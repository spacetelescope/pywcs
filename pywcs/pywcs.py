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

# $Id: ocumentation_Guidelines.html,v 1.2 2007/09/04 20:44:02 dencheva Exp $

__docformat__ = "epytext"

from _pywcs import *
import _pywcs
try:
    import pyfits
    has_pyfits = True
except ImportError:
    has_pyfits = False

if has_pyfits:
    # TODO: update formatting in this docstring
    def parse_wcs(hdulist):
        """Parses the WCS information in a PyFITS hdulist object.

All WCS keywords defined in Papers I, II, and III are
recognized, and also those used by the AIPS convention and certain
other keywords that existed in early drafts of the WCS papers.

Given a string containing a FITS image header, parse_wcs()
identifies and reads all WCS keywords for the primary coordinate
representation and up to 26 alternate representations.  It returns
this information as a list of Wcs objects.

parse_wcs() fills in information associated with coordinate lookup
tables.

Parameters
----------

header : hdulist
    A PyFITS hdulist parsed from a FITS file.

    TODO: Deal with this next paragraph in a Pythonic way
    For negative values of ctrl (see below), header[] is modified so
    that WCS keyrecords processed by wcspih() are removed from it.

Returns
-------

wcs : list
    A list of Wcs objects, containing up to 27 coordinate
    representations.

Other parameters
----------------

relax : string
    Degree of permissiveness:
        'standard': Recognize only FITS keywords defined by the
            published WCS standard.
        'all': Admit all recognized informal extensions of the
            WCS standard.

TODO: handle the ctrl/nreject arguments

Exceptions
----------

MemoryError
    """
        return parse_image_header(str(hdulist[0].header.ascardlist()))

# This is a hack so epydoc will document this method
def parse_image_header(header, relax=False):
    return _pywcs.parse_image_header(header, relax)
parse_image_header.__doc__ = _pywcs.parse_image_header.__doc__
