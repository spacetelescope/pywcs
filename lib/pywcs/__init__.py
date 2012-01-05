# Copyright (C) 2008-2012 Association of Universities for Research in Astronomy (AURA)

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

"""
.. _wcslib: http://www.atnf.csiro.au/~mcalabre/WCS/
.. _pyfits: http://www.stsci.edu/resources/software_hardware/pyfits
.. _Paper IV: http://www.atnf.csiro.au/people/mcalabre/WCS/index.html
.. _SIP: http://ssc.spitzer.caltech.edu/postbcd/doc/shupeADASS.pdf
.. _ds9: http://hea-www.harvard.edu/RD/ds9/

Pywcs provides transformations following the `SIP`_ conventions,
`Paper IV`_ table lookup distortion, and the core WCS functionality
provided by `wcslib`_.  Each of these transformations can be used
independently or together in a standard pipeline.

The basic workflow is as follows:

    1. ``import pywcs``

    2. Call the `pywcs.WCS` constructor with a `pyfits`_ header
       and/or hdulist object.

    3. Optionally, if the FITS file uses any deprecated or
       non-standard features, you may need to call one of the
       `~pywcs.WCS.fix` methods on the object.

    4. Use one of the following transformation methods:

       - `~WCS.all_pix2sky`: Perform all three transformations from
         pixel to sky coordinates.

       - `~WCS.wcs_pix2sky`: Perform just the core WCS transformation
         from pixel to sky coordinates.

       - `~WCS.wcs_sky2pix`: Perform just the core WCS transformation
         from sky to pixel coordinates.

       - `~WCS.sip_pix2foc`: Convert from pixel to focal plane
         coordinates using the `SIP`_ polynomial coefficients.

       - `~WCS.sip_foc2pix`: Convert from focal plane to pixel
         coordinates using the `SIP`_ polynomial coefficients.

       - `~WCS.p4_pix2foc`: Convert from pixel to focal plane
         coordinates using the table lookup distortion method
         described in `Paper IV`_.

       - `~WCS.det2im`: Convert from detector coordinates to image
         coordinates.  Commonly used for narrow column correction.
"""

from __future__ import division # confidence high

import sys
if sys.version_info[0] >= 3:
    exec("from .pywcs import *")
else:
    from pywcs import *

__version__ = "1.11-4.7"
