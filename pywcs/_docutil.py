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

def ONE_OR_TWO_ARGS(out_type, indent=0):
    return _fix("""Either one or two arguments may be provided.

    - one argument: An Nx2 array of I{x}- and I{y}-coordinates.

    - two arguments: Two one-dimensional arrays of I{x} and I{y}
          coordinates.

@return: Returns the %s coordinates.  If the input was a
     single array, a single array is returned, otherwise a tuple of
     arrays is returned.""" % out_type, indent)

def FITS_EQUIVALENT(method_name, indent=0):
    return _fix("""B{The pixel coordinates are 0-based (like array indices in C and
Python).  If your pixel coordinates are 1-based (like array
indices in Fortran), use L{%s_fits} instead.}""" % method_name, indent)

def NON_FITS_EQUIVALENT(method_name, indent=0):
    return _fix("""Identical to L{%s}, except pixel coordinates are 1-based (like array
indices in Fortran), instead of 0-based (like array indices C and
Python).""" % method_name, indent)

