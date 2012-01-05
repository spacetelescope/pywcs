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
pywcs-specific utilities for generating boilerplate in docstrings.
"""

from __future__ import division # confidence high

def _fix(content, indent=0):
    lines = content.split('\n')
    indent = '\n' + ' ' * indent
    return indent.join(lines)

def TWO_OR_THREE_ARGS(out_type, naxis, indent=0):
    return _fix(
"""Either two or three arguments may be provided.

    - 2 arguments: An *N* x *%s* array of *x*- and *y*-coordinates, and
      an *origin*.

    - 3 arguments: 2 one-dimensional arrays of *x* and *y*
      coordinates, and an *origin*.

Here, *origin* is the coordinate in the upper left corner of the
image.  In FITS and Fortran standards, this is 1.  In Numpy and C
standards this is 0.

Returns the %s.  If the input was a single array and
origin, a single array is returned, otherwise a tuple of arrays is
returned.""" % (naxis, out_type), indent)

def ORIGIN(indent=0):
    return _fix(
"""
- *origin*: int. Specifies the origin of pixel values.  The Fortran and
  FITS standards use an origin of 1.  Numpy and C use array indexing
  with origin at 0.
""", indent)

def RA_DEC_ORDER(indent=0):
    return _fix(
"""
An optional keyword argument, *ra_dec_order*, may be provided, that
when `True` will ensure that sky coordinates are always given and
returned in as (*ra*, *dec*) pairs, regardless of the order of the
axes specified by the in the ``CTYPE`` keywords.
""", indent)

