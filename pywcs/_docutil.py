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
pywcs-specific utilities for generating boilerplate in docstrings.
"""

def _fix(content, indent=0):
    lines = content.split('\n')
    indent = '\n' + ' ' * indent
    return indent.join(lines)

def TWO_OR_THREE_ARGS(out_type, indent=0):
    return _fix(
"""Either two or three arguments may be provided.

    - two arguments: An Nx2 array of I{x}- and I{y}-coordinates, and
      an origin

    - three arguments: Two one-dimensional arrays of I{x} and I{y}
      coordinates, and an origin.

Here, origin is the coordinate in the upper left corner of the image.
In FITS/Fortran standards, this is 1.  In Numpy/C standards this is 0.

@return: Returns the %s coordinates.  If the input was a single array
     and origin, a single array is returned, otherwise a tuple of
     arrays is returned.""" % out_type, indent)

def ORIGIN(indent=0):
    return _fix(
"""
@param origin: Specifies the origin of pixel values.  The Fortran and
FITS standards use an origin of 1.  Numpy and C use array indexing
with origin at 0.

@type origin: int
""", indent)

