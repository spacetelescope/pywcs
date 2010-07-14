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

###########################################################################

from __future__ import division # confidence high
del division
# We don't want the "division" symbol in the namespace, since it
# should have only docstrings

# It gets to be really tedious to type long docstrings in ANSI C
# syntax (since multi-line string literals are not valid).
# Therefore, the docstrings are written here in doc/docstrings.py,
# which are then converted by setup.py into docstrings.h, which is
# included by pywcs.c

import _docutil as __

a = """
``double array[a_order+1][a_order+1]``

The `SIP`_ ``A_i_j`` matrix used for pixel to focal plane
transformation.

Its values may be changed in place, but it may not be resized, without
creating a new `~pywcs.Sip` object.
"""

a_order = """
``int`` (read-only)

The order of the polynomial in the `SIP`_ ``A_i_j`` array (``A_ORDER``).
"""

all_pix2sky = """
all_pix2sky(pixcrd, origin) -> ``double array[ncoord][nelem]``

Transforms pixel coordinates to sky coordinates by doing all of the
following:

    - Detector to image plane correction (optionally)

    - SIP distortion correction (optionally)

    - Paper IV distortion correction (optionally)

    - wcslib WCS transformation

The first three (the distortion corrections) are done in parallel.

- *pixcrd*: double array[ncoord][nelem].  Array of pixel coordinates.

%s

Returns an array of sky coordinates.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `SingularMatrixError`: Linear transformation matrix is singular.

- `InconsistentAxisTypesError`: Inconsistent or unrecognized
  coordinate axis types.

- `ValueError`: Invalid parameter value.

- `ValueError`: Invalid coordinate transformation parameters.

- `ValueError`: x- and y-coordinate arrays are not the same size.

- `InvalidTransformError`: Invalid coordinate transformation.

- `InvalidTransformError`: Ill-conditioned coordinate transformation
  parameters.
""" % __.ORIGIN()

alt = """
``str``

Character code for alternate coordinate descriptions.  For example,
the ``"a"`` in keyword names such as ``CTYPEia``.  This is a space
character for the primary coordinate description, or one of the 26
upper-case letters, A-Z.
"""

ap = """
``double array[ap_order+1][ap_order+1]``

The `SIP`_ ``AP_i_j`` matrix used for focal plane to pixel
transformation.  Its values may be changed in place, but it may not be
resized, without creating a new `~pywcs.Sip` object.
"""

ap_order = """
``int`` (read-only)

The order of the polynomial in the `SIP`_ ``AP_i_j`` array
(``AP_ORDER``).
"""

axis_types = """
``int array[naxis]``

An array of four-digit type codes for each axis.

- First digit (i.e. 1000s):

  - 0: Non-specific coordinate type.

  - 1: Stokes coordinate.

  - 2: Celestial coordinate (including ``CUBEFACE``).

  - 3: Spectral coordinate.

- Second digit (i.e. 100s):

  - 0: Linear axis.

  - 1: Quantized axis (``STOKES``, ``CUBEFACE``).

  - 2: Non-linear celestial axis.

  - 3: Non-linear spectral axis.

  - 4: Logarithmic axis.

  - 5: Tabular axis.

- Third digit (i.e. 10s):

  - 0: Group number, e.g. lookup table number

- The fourth digit is used as a qualifier depending on the axis type.

  - For celestial axes:

    - 0: Longitude coordinate.

    - 1: Latitude coordinate.

    - 2: ``CUBEFACE`` number.

  - For lookup tables: the axis number in a multidimensional table.

``CTYPEia`` in ``"4-3"`` form with unrecognized algorithm code will
have its type set to -1 and generate an error.
"""

b = """
``double array[b_order+1][b_order+1]``

The `SIP`_ ``B_i_j`` matrix used for pixel to focal plane
transformation.  Its values may be changed in place, but it may not be
resized, without creating a new `~pywcs.Sip` object.
"""

b_order = """
``int`` (read-only)

The order of the polynomial in the `SIP`_ ``B_i_j`` array
(``B_ORDER``).
"""

bp = """
``double array[bp_order+1][bp_order+1]``

The `SIP`_ ``BP_i_j`` matrix used for focal plane to pixel
transformation.  Its values may be changed in place, but it may not be
resized, without creating a new `~pywcs.Sip` object.
"""

bp_order = """
``int`` (read-only)

The order of the polynomial in the `SIP`_ ``BP_i_j`` array
(``BP_ORDER``).
"""

cd = """
``double array[2][2]``

The ``CDi_ja`` linear transformation matrix.

For historical compatibility, two alternate specifications of the
``CDi_ja`` and ``CROTAia`` keywords are supported.  Although these may
not formally co-exist with ``PCi_ja``, the approach here is simply to
ignore them if given in conjunction with ``PCi_ja``.

`~pywcs.Wcsprm.has_pci_ja`, `~pywcs.Wcsprm.has_cdi_ja` and
`~pywcs.Wcsprm.has_crotaia` can be used to determine which of these
alternatives are present in the header.

``CDi_ja`` and ``CROTAia`` keywords, if found, are to be stored in the
`~pywcs.Wcsprm.cd` and `~pywcs.Wcsprm.crota` arrays which are
dimensioned similarly to `~pywcs.Wcsprm.pc` and `~pywcs.Wcsprm.cdelt`.

These alternate specifications of the linear transformation matrix are
translated immediately to ``PCi_ja`` by `~pywcs.Wcsprm.set` and are
nowhere visible to the lower-level routines.  In particular,
`~pywcs.Wcsprm.set` resets `~pywcs.Wcsprm.cdelt` to unity if
``CDi_ja`` is present (and no ``PCi_ja``).  If no ``CROTAia`` is
associated with the latitude axis, `set` reverts to a unity ``PCi_ja``
matrix.
"""

cdelt = """
``double array[naxis]``

Coordinate increments (``CDELTia``) for each coord axis.

If a ``CDi_ja`` linear transformation matrix is present, a warning is
raised and `~pywcs.Wcsprm.cdelt` is ignored.  The ``CDi_ja`` matrix
may be deleted by::

  del wcs.wcs.cd

An undefined value is represented by NaN.
"""

cel_offset = """
``boolean``

If ``True``, an offset will be applied to ``(x, y)`` to force ``(x,y)
= (0,0)`` at the fiducial point.
"""

celfix = """
Translates AIPS-convention celestial projection types, ``-NCP`` and
``-GLS``.

Returns ``0`` for success; ``-1`` if no change required.
"""

cname = """
``list of strings``

A list of the coordinate axis names, from ``CNAMEia``.
"""

colax = """
``int array[naxis]``

An array recording the column numbers for each axis in a pixel list.
"""

colnum = """
``int``

Where the coordinate representation is associated with an image-array
column in a FITS binary table, this property may be used to record the
relevant column number.

It should be set to zero for an image header or pixel list.
"""

copy = """
Creates a deep copy of the WCS object.
"""

cpdis1 = """
`~pywcs.DistortionLookupTable`

The pre-linear transformation distortion lookup table, ``CPDIS1``.
"""

cpdis2 = """
`~pywcs.DistortionLookupTable`

The pre-linear transformation distortion lookup table, ``CPDIS2``.
"""

crder = """
``double array[naxis]``

The random error in each coordinate axis, ``CRDERia``.

An undefined value is represented by NaN.
"""

crota = """
``double array[naxis]``

``CROTAia`` keyvalues for each coordinate axis.

``CROTAia`` is an alternate specification of the linear transformation
matrix, maintained for historical compatibility.  See `cd` for the
current method.
"""

crpix = """
``double array[naxis]``

Coordinate reference pixels (``CRPIXja``) for each pixel axis.
"""

crval = """
``double array[naxis]``

Coordinate reference values (``CRVALia``) for each coordinate axis.
"""

csyer = """
``double array[naxis]``

The systematic error in the coordinate value axes, ``CSYERia``.

An undefined value is represented by NaN.
"""

ctype = """
``list of strings[naxis]``

List of ``CTYPEia`` keyvalues.

The `~pywcs.Wcsprm.ctype` keyword values must be in upper case and
there must be zero or one pair of matched celestial axis types, and
zero or one spectral axis.
"""

cubeface = """
``int``

Index into the ``pixcrd`` (pixel coordinate) array for the
``CUBEFACE`` axis.  This is used for quadcube projections where the
cube faces are stored on a separate axis.

The quadcube projections (``TSC``, ``CSC``, ``QSC``) may be represented
in FITS in either of two ways:

    - The six faces may be laid out in one plane and numbered as
      follows::


                                       0

                              4  3  2  1  4  3  2

                                       5

      Faces 2, 3 and 4 may appear on one side or the other (or both).
      The sky-to-pixel routines map faces 2, 3 and 4 to the left but
      the pixel-to-sky routines accept them on either side.

    - The ``COBE`` convention in which the six faces are stored in a
      three-dimensional structure using a ``CUBEFACE`` axis indexed
      from 0 to 5 as above.

These routines support both methods; `~pywcs.Wcsprm.set` determines
which is being used by the presence or absence of a ``CUBEFACE`` axis
in `~pywcs.Wcsprm.ctype`.  `~pywcs.Wcsprm.p2s` and `~pywcs.Wcsprm.s2p`
translate the ``CUBEFACE`` axis representation to the single plane
representation understood by the lower-level projection routines.
"""

cunit = """
``list of strings[naxis]``

List of ``CUNITia`` keyvalues which define the units of measurement of
the ``CRVALia``, ``CDELTia`` and ``CDi_ja`` keywords.

As ``CUNITia`` is an optional header keyword, `~pywcs.Wcsprm.cunit`
may be left blank but otherwise is expected to contain a standard
units specification as defined by WCS Paper I.  Utility function
`wcsutrn`, (not currently wrapped for Python) is available to
translate commonly used non-standard units specifications but this
must be done as a separate step before invoking `~pywcs.Wcsprm.set`.

For celestial axes, if `~pywcs.Wcsprm.cunit` is not blank,
`~pywcs.Wcsprm.set` uses `wcsunits` to parse it and scale
`~pywcs.Wcsprm.cdelt`, `~pywcs.Wcsprm.crval`, and `~pywcs.Wcsprm.cd`
to decimal degrees.  It then resets `~pywcs.Wcsprm.cunit` to
``"deg"``.

For spectral axes, if `~pywcs.Wcsprm.cunit` is not blank,
`~pywcs.Wcsprm.set` uses `wcsunits` to parse it and scale
`~pywcs.Wcsprm.cdelt`, `~pywcs.Wcsprm.crval`, and `~pywcs.Wcsprm.cd`
to SI units.  It then resets `~pywcs.Wcsprm.cunit` accordingly.

`~pywcs.Wcsprm.set` ignores `~pywcs.Wcsprm.cunit` for other coordinate
types; `~pywcs.Wcsprm.cunit` may be used to label coordinate values.
"""

cylfix = """
Fixes WCS keyvalues for malformed cylindrical projections.

Returns ``0`` for success; ``-1`` if no change required.
"""

data = """
``float array``

The array data for the `~pywcs.DistortionLookupTable`.
"""

dateavg = """
``string``

Representative mid-point of the date of observation in ISO format,
``yyyy-mm-ddThh:mm:ss``.

.. seealso::

   `~pywcs.Wcsprm.dateobs`
"""

dateobs = """
``string``

Start of the date of observation in ISO format,
``yyyy-mm-ddThh:mm:ss``.

.. seealso::

   `~pywcs.Wcsprm.dateavg`
"""

datfix = """
Translates the old ``DATE-OBS`` date format to year-2000 standard form
``(yyyy-mm-ddThh:mm:ss)`` and derives ``MJD-OBS`` from it if not
already set.  Alternatively, if `~pywcs.Wcsprm.mjdobs` is set and
`~pywcs.Wcsprm.dateobs` isn't, then `~pywcs.Wcsprm.datfix` derives
`~pywcs.Wcsprm.dateobs` from it.  If both are set but disagree by more
than half a day then `ValueError` is raised.

Returns ``0`` for success; ``-1`` if no change required.
"""

det2im = """
Convert detector coordinates to image plane coordinates.
"""

det2im1 = """
A `~pywcs.DistortionLookupTable` object for detector to image plane
correction in the *x*-axis.
"""

det2im2 = """
A `~pywcs.DistortionLookupTable` object for detector to image plane
correction in the *y*-axis.
"""

DistortionLookupTable = """
DistortionLookupTable(*table*, *crpix*, *crval*, *cdelt*)

- *table*: 2-dimensional array for the distortion lookup table.

- *crpix*: the distortion array reference pixel (a 2-tuple)

- *crval*: is the image array pixel coordinate (a 2-tuple)

- *cdelt*: is the grid step size (a 2-tuple)

Represents a single lookup table for a `Paper IV`_ distortion
transformation.
"""

equinox = """
``double``

The equinox associated with dynamical equatorial or ecliptic
coordinate systems, ``EQUINOXa`` (or ``EPOCH`` in older headers).  Not
applicable to ICRS equatorial or ecliptic coordinates.

An undefined value is represented by NaN.
"""

fix = """
fix(*translate_units=''*, *naxis=0*)

Applies all of the corrections handled separately by
`~pywcs.Wcsprm.datfix`, `~pywcs.Wcsprm.unitfix`,
`~pywcs.Wcsprm.celfix`, `~pywcs.Wcsprm.spcfix` and
`~pywcs.Wcsprm.cylfix`.

- *translate_units*: string. Do potentially unsafe translations of
  non-standard unit strings.

  Although ``"S"`` is commonly used to represent seconds, its
  translation to ``"s"`` is potentially unsafe since the standard
  recognizes ``"S"`` formally as Siemens, however rarely that may be
  used.  The same applies to ``"H"`` for hours (Henry), and ``"D"``
  for days (Debye).

  This string controls what to do in such cases, and is
  case-insensitive.

  - If the string contains ``"s"``, translate ``"S"`` to ``"s"``.

  - If the string contains ``"h"``, translate ``"H"`` to ``"h"``.

  - If the string contains ``"d"``, translate ``"D"`` to ``"d"``.

    Thus ``''`` doesn't do any unsafe translations, whereas ``'shd'``
    does all of them.

- *naxis*: int array[naxis].  Image axis lengths.  If this array is
  set to zero or ``None``, then `~pywcs.Wcsprm.cylfix` will not be
  invoked.

Returns a dictionary containing the following keys, each referring to
a status string for each of the sub-fix functions that were
called:

- `~pywcs.Wcsprm.datfix`

- `~pywcs.Wcsprm.unitfix`

- `~pywcs.Wcsprm.celfix`

- `~pywcs.Wcsprm.spcfix`

- `~pywcs.Wcsprm.cylfix`
"""

get_offset = """
get_offset(*x, y*) -> (*x, y*)

Returns the offset from the distortion table for pixel point (*x, y*).
"""

get_ps = """
get_ps() -> list of tuples

Returns ``PSi_ma`` keywords for each *i* and *m*.  Returned as a list
of tuples of the form (*i*, *m*, *value*):

    - *i*: int.  Axis number, as in ``PSi_ma``, (i.e. 1-relative)

    - *m*: int.  Parameter number, as in ``PSi_ma``, (i.e. 0-relative)

    - *value*: string.  Parameter value.

.. seealso::

   `~pywcs.Wcsprm.set_ps`
"""

get_pv = """
get_pv() -> list of tuples

Returns ``PVi_ma`` keywords for each *i* and *m*.  Returned as a list
of tuples of the form (*i*, *m*, *value*):

    - *i*: int.  Axis number, as in ``PVi_ma``, (i.e. 1-relative)

    - *m*: int.  Parameter number, as in ``PVi_ma``, (i.e. 0-relative)

    - *value*: string. Parameter value.

Note that, if they were not given, `~pywcs.Wcsprm.set` resets the
entries for ``PVi_1a``, ``PVi_2a``, ``PVi_3a``, and ``PVi_4a`` for
longitude axis *i* to match (``phi_0``, ``theta_0``), the native
longitude and latitude of the reference point given by ``LONPOLEa``
and ``LATPOLEa``.

.. seealso::

   `~pywcs.Wcsprm.set_pv`
"""

has_cdi_ja = """
has_cdi_ja() -> bool

Returns ``True`` if ``CDi_ja`` is present.  ``CDi_ja`` is an alternate
specification of the linear transformation matrix, maintained for
historical compatibility.

Matrix elements in the IRAF convention are equivalent to the product
``CDi_ja = CDELTia * PCi_ja``, but the defaults differ from that of
the ``PCi_ja`` matrix.  If one or more ``CDi_ja`` keywords are present
then all unspecified ``CDi_ja`` default to zero.  If no ``CDi_ja`` (or
``CROTAia``) keywords are present, then the header is assumed to be in
``PCi_ja`` form whether or not any ``PCi_ja`` keywords are present
since this results in an interpretation of ``CDELTia`` consistent with
the original FITS specification.

While ``CDi_ja`` may not formally co-exist with ``PCi_ja``, it may
co-exist with ``CDELTia`` and ``CROTAia`` which are to be ignored.

.. seealso::

   `cd`
"""

has_crotaia = """
has_crotaia() -> bool

Returns ``True`` if ``CROTAia`` is present.  ``CROTAia`` is an
alternate specification of the linear transformation matrix,
maintained for historical compatibility.

In the AIPS convention, ``CROTAia`` may only be associated with the
latitude axis of a celestial axis pair.  It specifies a rotation in
the image plane that is applied *after* the ``CDELTia``; any other
``CROTAia`` keywords are ignored.

``CROTAia`` may not formally co-exist with ``PCi_ja``.  ``CROTAia`` and
``CDELTia`` may formally co-exist with ``CDi_ja`` but if so are to be
ignored.

.. seealso::

   `cd`
"""

has_pci_ja = """
has_pci_ja() -> bool

Returns ``True`` if ``PCi_ja`` is present.  ``PCi_ja`` is the
recommended way to specify the linear transformation matrix.

.. seealso::

   `cd`
"""

imgpix_matrix = """
``double array[2][2]`` (read-only)

Inverse of the matrix containing the product of the ``CDELTia``
diagonal matrix and the ``PCi_ja`` matrix.

.. warning::

   This value may not be correct until after `~pywcs.Wcsprm.set` is
   called.
"""

is_unity = """
is_unity() -> bool

Returns ``True`` if the linear transformation matrix
(`~pywcs.Wcsprm.cd`) is unity.

.. warning::

   This value may not be correct until after `~pywcs.Wcsprm.set` is
   called.
"""

lat = """
``int`` (read-only)

The index into the sky coordinate array containing latitude values.
"""

latpole = """
``double``

The native latitude of the celestial pole, ``LATPOLEa`` (deg).
"""

lattyp = """
``string`` (read-only)

Celestial axis type for longitude, e.g. "RA", "DEC", "GLON", "GLAT",
etc. extracted from 'RA--', 'DEC-', 'GLON', 'GLAT', etc. in the first
four characters of ``CTYPEia`` but with trailing dashes removed.

.. warning::

   This value may not be correct until after `~pywcs.Wcsprm.set` is
   called.
"""

lng = """
``int`` (read-only)

The index into the sky coordinate array containing longitude values.
"""

lngtyp = """
``string`` (read-only)

Celestial axis type for longitude, e.g. "RA", "DEC", "GLON", "GLAT",
etc. extracted from 'RA--', 'DEC-', 'GLON', 'GLAT', etc. in the first
four characters of ``CTYPEia`` but with trailing dashes removed.

.. warning::

   This value may not be correct until after `~pywcs.Wcsprm.set` is
   called.
"""

lonpole = """
``double``

The native longitude of the celestial pole, ``LONPOLEa`` (deg).
"""

mix = """
mix(*mixpix, mixcel, vspan, vstep, viter, world, pixcrd, origin*) -> dict

Given either the celestial longitude or latitude plus an element of
the pixel coordinate, solves for the remaining elements by iterating
on the unknown celestial coordinate element using `~pywcs.Wcsprm.s2p`.

- *mixpix*: int.  Which element on the pixel coordinate is given.

- *mixcel*: int.  Which element of the celestial coordinate is
  given. If mixcel* = ``1``, celestial longitude is given in
  ``world[self.lng]``, latitude returned in ``world[self.lat]``.  If
  *mixcel* = ``2``, celestial latitude is given in
  ``world[self.lat]``, longitude returned in ``world[self.lng]``.

- *vspan*: pair of floats.  Solution interval for the celestial
  coordinate, in degrees.  The ordering of the two limits is
  irrelevant.  Longitude ranges may be specified with any convenient
  normalization, for example ``(-120,+120)`` is the same as
  ``(240,480)``, except that the solution will be returned with the
  same normalization, i.e. lie within the interval specified.

- *vstep*: float.  Step size for solution search, in degrees.  If
  ``0``, a sensible, although perhaps non-optimal default will be
  used.

- *viter*: int.  If a solution is not found then the step size will be
  halved and the search recommenced.  *viter* controls how many times
  the step size is halved.  The allowed range is 5 - 10.

- *world*: double array[naxis].  World coordinate elements.
  ``world[self.lng]`` and ``world[self.lat]`` are the celestial
  longitude and latitude, in degrees.  Which is given and which
  returned depends on the value of *mixcel*.  All other elements are
  given.  The results will be written to this array in-place.

- *pixcrd*: double array[naxis].  Pixel coordinate.  The element
  indicated by *mixpix* is given and the remaining elements will be
  written in-place.

%s

Returns dictionary with the following keys:

- *phi* (double array[naxis])

- *theta* (double array[naxis])

  - Longitude and latitude in the native coordinate system of the
    projection, in degrees.

- *imgcrd* (double array[naxis])

  - Image coordinate elements.  ``imgcrd[self.lng]`` and
    ``imgcrd[self.lat]`` are the projected *x*- and *y*-coordinates,
    in decimal degrees.

- *world* (double array[naxis])

  - Another reference to the *world* argument passed in.

**Exceptions:**

- `MemoryError` Memory allocation failed.

- `SingularMatrixError`: Linear transformation matrix is singular.

- `InconsistentAxisTypesError`: Inconsistent or unrecognized
  coordinate axis types.

- `ValueError`: Invalid parameter value.

- `InvalidTransformError`: Invalid coordinate transformation
  parameters.

- `InvalidTransformError` Ill-conditioned coordinate transformation
  parameters.

- `InvalidCoordinateError`: Invalid world coordinate.

- `NoSolutionError`: No solution found in the specified interval.

.. seealso::

   `~pywcs.Wcsprm.lat`, `~pywcs.Wcsprm.lng`

.. note::

  Initially, the specified solution interval is checked to see if it's
  a "crossing" interval.  If it isn't, a search is made for a crossing
  solution by iterating on the unknown celestial coordinate starting
  at the upper limit of the solution interval and decrementing by the
  specified step size.  A crossing is indicated if the trial value of
  the pixel coordinate steps through the value specified.  If a
  crossing interval is found then the solution is determined by a
  modified form of "regula falsi" division of the crossing interval.
  If no crossing interval was found within the specified solution
  interval then a search is made for a "non-crossing" solution as may
  arise from a point of tangency.  The process is complicated by
  having to make allowance for the discontinuities that occur in all
  map projections.

  Once one solution has been determined others may be found by
  subsequent invokations of `~pywcs.Wcsprm.mix` with suitably
  restricted solution intervals.

  Note the circumstance that arises when the solution point lies at a
  native pole of a projection in which the pole is represented as a
  finite curve, for example the zenithals and conics.  In such cases
  two or more valid solutions may exist but `~pywcs.Wcsprm.mix` only
  ever returns one.

  Because of its generality, `~pywcs.Wcsprm.mix` is very
  compute-intensive.  For compute-limited applications, more efficient
  special-case solvers could be written for simple projections, for
  example non-oblique cylindrical projections.
""" % __.ORIGIN()

mjdavg = """
``double``

Modified Julian Date ``(MJD = JD - 2400000.5)``, ``MJD-AVG``,
corresponding to ``DATE-AVG``.

An undefined value is represented by NaN.

.. seealso::

   `~pywcs.Wcsprm.mjdobs`
"""

mjdobs = """
``double``

Modified Julian Date ``(MJD = JD - 2400000.5)``, ``MJD-OBS``,
corresponding to ``DATE-OBS``.

An undefined value is represented by NaN.

.. seealso::

   `~pywcs.Wcsprm.mjdavg`
"""

name = """
``string``

The name given to the coordinate representation ``WCSNAMEa``.
"""

naxis = """
``int`` (read-only)

The number of axes (pixel and coordinate), given by the ``NAXIS`` or
``WCSAXESa`` keyvalues.

The number of coordinate axes is determined at parsing time, and can
not be subsequently changed.

It is determined from the highest of the following:

  1. ``NAXIS``

  2. ``WCSAXESa``

  3. The highest axis number in any parameterized WCS keyword.  The
     keyvalue, as well as the keyword, must be syntactically valid
     otherwise it will not be considered.

If none of these keyword types is present, i.e. if the header only
contains auxiliary WCS keywords for a particular coordinate
representation, then no coordinate description is constructed for it.

This value may differ for different coordinate representations of the
same image.
"""

obsgeo = """
``double array[3]``

Location of the observer in a standard terrestrial reference frame,
``OBSGEO-X``, ``OBSGEO-Y``, ``OBSGEO-Z`` (in meters).

An undefined value is represented by NaN.
"""

p2s = """
p2s(*pixcrd, origin*) -> dict

Converts pixel to sky coordinates.

- *pixcrd*: double array[ncoord][nelem].  Array of pixel coordinates.

%s

Returns a dictionary with the following keys:

- *imgcrd*: double array[ncoord][nelem]

  - Array of intermediate sky coordinates.  For celestial axes,
    ``imgcrd[][self.lng]`` and ``imgcrd[][self.lat]`` are the
    projected *x*-, and *y*-coordinates, in pseudo degrees.  For
    spectral axes, ``imgcrd[][self.spec]`` is the intermediate
    spectral coordinate, in SI units.

- *phi*: double array[ncoord]

- *theta*: double array[ncoord]

  - Longitude and latitude in the native coordinate system of the
    projection, in degrees.

- *world*: double array[ncoord][nelem]

  - Array of sky coordinates.  For celestial axes,
    ``world[][self.lng]`` and ``world[][self.lat]`` are the celestial
    longitude and latitude, in degrees.  For spectral axes,
    ``world[][self.spec]`` is the intermediate spectral coordinate, in
    SI units.

- *stat*: int array[ncoord]

  - Status return value for each coordinate. ``0`` for success,
    ``1+`` for invalid pixel coordinate.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `SingularMatrixError`: Linear transformation matrix is singular.

- `InconsistentAxisTypesError`: Inconsistent or unrecognized
  coordinate axis types.

- `ValueError`: Invalid parameter value.

- `ValueError`: *x*- and *y*-coordinate arrays are not the same size.

- `InvalidTransformError`: Invalid coordinate transformation
  parameters.

- `InvalidTransformError`: Ill-conditioned coordinate transformation
  parameters.

.. seealso::

   `~pywcs.Wcsprm.lat`, `~pywcs.Wcsprm.lng`
""" % __.ORIGIN()

p4_pix2foc = """
p4_pix2foc(*pixcrd, origin*) -> double array[ncoord][nelem]

Convert pixel coordinates to focal plane coordinates using `Paper IV`_
lookup-table distortion correction.

- *pixcrd*: double array[ncoord][nelem].  Array of pixel coordinates.

%s

Returns an array of focal plane coordinates.

- `MemoryError`: Memory allocation failed.

- `ValueError`: Invalid coordinate transformation parameters.
""" % __.ORIGIN()

pc = """
``double array[2][2]``

The ``PCi_ja`` (pixel coordinate) transformation matrix.  The order is::

  [[PC1_1, PC1_2],
   [PC2_1, PC2_2]]

For historical compatibility, two alternate specifications of the
``CDi_ja`` and ``CROTAia`` keywords are supported.

`~pywcs.Wcsprm.has_pci_ja`, `~pywcs.Wcsprm.has_cdi_ja` and
`~pywcs.Wcsprm.has_crotaia` can be used to determine which of these
alternatives are present in the header.
"""

phi0 = """
``double``

The native latitude of the fiducial point, i.e. the point whose
celestial coordinates are given in ``ref[1:2]``.  If undefined (NaN)
the initialization routine, `~pywcs.Wcsprm.set`, will set this to a
projection-specific default.

.. seealso::

   `~pywcs.Wcsprm.theta0`
"""

pix2foc = """
pix2foc(*pixcrd, origin*) -> double array[ncoord][nelem]

Perform both `SIP`_ polynomial and `Paper IV`_ lookup-table distortion
correction in parallel.

- *pixcrd*: double array[ncoord][nelem].  Array of pixel coordinates.

%s

Returns an array of focal plane coordinates.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `ValueError`: Invalid coordinate transformation parameters.
""" % __.ORIGIN()

piximg_matrix = """
``double array[2][2]`` (read-only)

Matrix containing the product of the ``CDELTia`` diagonal matrix and
the ``PCi_ja`` matrix.

.. warning::

   This value may not be correct until after `~pywcs.Wcsprm.set` is
   called.
"""

print_contents = """
print_contents()

Print the contents of the WCS object to stdout.  Probably only useful
for debugging purposes, and may be removed in the future.
"""

radesys = """
``string``

The equatorial or ecliptic coordinate system type, ``RADESYSa``.
"""

restfrq = """
``double``

Rest frequency (Hz) from ``RESTFRQa``.

An undefined value is represented by NaN.
"""

restwav = """
``double``

Rest wavelength (m) from ``RESTWAVa``.

An undefined value is represented by NaN.
"""

s2p = """
s2p(*sky, origin*) -> dict

Transforms sky coordinates to pixel coordinates.

- *sky*: double array[ncoord][nelem].  Array of sky coordinates, in
  decimal degrees.

%s

Returns a dictionary with the following keys:

- *phi*: double array[ncoord]

- *theta*: double array[ncoord]

  - Longitude and latitude in the native coordinate system of the
    projection, in degrees.

- *imgcrd*: double array[ncoord][nelem]

  - Array of intermediate sky coordinates.  For celestial axes,
    ``imgcrd[][self.lng]`` and ``imgcrd[][self.lat]`` are the
    projected *x*-, and *y*-coordinates, in pseudo "degrees".  For
    quadcube projections with a ``CUBEFACE`` axis, the face number is
    also returned in ``imgcrd[][self.cubeface]``.  For spectral axes,
    ``imgcrd[][self.spec]`` is the intermediate spectral coordinate,
    in SI units.

- *pixcrd*: double array[ncoord][nelem]

  - Array of pixel coordinates.  Pixel coordinates are
    zero-based.

- *stat*: int array[ncoord]

  - Status return value for each coordinate. ``0`` for
    success, ``1+`` for invalid pixel coordinate.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `SingularMatrixError`: Linear transformation matrix is singular.

- `InconsistentAxisTypesError` Inconsistent or unrecognized coordinate
  axis types.

- `ValueError`: Invalid parameter value.

- `InvalidTransformError`: Invalid coordinate transformation
  parameters.

- `InvalidTransformError`: Ill-conditioned coordinate transformation
  parameters.

.. seealso::

   `~pywcs.Wcsprm.lat`, `~pywcs.Wcsprm.lng`
""" % (__.ORIGIN())

set = """
set()

Sets up a WCS object for use according to information supplied within
it.

Note that this routine need not be called directly; it will be invoked
by `~pywcs.Wcsprm.p2s` and `~pywcs.Wcsprm.s2p` if necessary.

Some attributes that are based on other attributes (such as
`~pywcs.Wcsprm.lattyp` on `~pywcs.Wcsprm.ctype`) may not be correct
until after `~pywcs.Wcsprm.set` is called.

`~pywcs.Wcsprm.set` strips off trailing blanks in all string members.

`~pywcs.Wcsprm.set` recognizes the ``NCP`` projection and converts it
to the equivalent ``SIN`` projection and it also recognizes ``GLS`` as
a synonym for ``SFL``.  It does alias translation for the AIPS
spectral types (``FREQ-LSR``, ``FELO-HEL``, etc.) but without changing
the input header keywords.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `SingularMatrixError`: Linear transformation matrix is singular.

- `InconsistentAxisTypesError`: Inconsistent or unrecognized
  coordinate axis types.

- `ValueError`: Invalid parameter value.

- `InvalidTransformError`: Invalid coordinate transformation
  parameters.

- `InvalidTransformError`: Ill-conditioned coordinate transformation
  parameters.
"""

set_ps = """
set_ps(list)

Sets `PSi_ma` keywords for each *i* and *m*.  The input must be a
sequence of tuples of the form (*i*, *m*, *value*):

    - *i*: int.  Axis number, as in ``PSi_ma``, (i.e. 1-relative)

    - *m*: int.  Parameter number, as in ``PSi_ma``, (i.e. 0-relative)

    - *value*: string.  Parameter value.

.. seealso::

   `~pywcs.Wcsprm.get_ps`
"""

set_pv = """
set_pv(list)

Sets `PVi_ma` keywords for each *i* and *m*.  The input must be a
sequence of tuples of the form (*i*, *m*, *value*):

    - *i*: int.  Axis number, as in ``PVi_ma``, (i.e. 1-relative)

    - *m*: int.  Parameter number, as in ``PVi_ma``, (i.e. 0-relative)

    - *value*: string.  Parameter value.

.. seealso::

   `~pywcs.Wcsprm.get_pv`
"""

sip = """
Get/set the `~pywcs.Sip` object for performing `SIP`_ distortion
correction.
"""

Sip = """
Sip(*a, b, ap, bp, crpix*)

The `~pywcs.Sip` class performs polynomial distortion correction using
the `SIP`_ convention in both directions.

   Shupe, D. L., M. Moshir, J. Li, D. Makovoz and R. Narron.  2005.
   "The SIP Convention for Representing Distortion in FITS Image
   Headers."  ADASS XIV.

- *a*: double array[m+1][m+1].  The ``A_i_j`` polynomial for pixel to
  focal plane transformation.  Its size must be (*m* + 1, *m* + 1)
  where *m* = ``A_ORDER``.

- *b*: double array[m+1][m+1].  The ``B_i_j`` polynomial for pixel to
  focal plane transformation.  Its size must be (*m* + 1, *m* + 1)
  where *m* = ``B_ORDER``.

- *ap*: double array[m+1][m+1].  The ``AP_i_j`` polynomial for pixel
  to focal plane transformation.  Its size must be (*m* + 1, *m* + 1)
  where *m* = ``AP_ORDER``.

- *bp*: double array[m+1][m+1].  The ``BP_i_j`` polynomial for pixel to
  focal plane transformation.  Its size must be (*m* + 1, *m* + 1) where
  *m* = ``BP_ORDER``.

- *crpix*: double array[2].  The reference pixel.
"""

sip_foc2pix = """
sip_foc2pix(*foccrd, origin*) -> double array[ncoord][nelem]

Convert focal plane coordinates to pixel coordinates using the `SIP`_
polynomial distortion convention.

- *foccrd*: double array[ncoord][nelem].  Array of focal plane
  coordinates.

%s

Returns an array of pixel coordinates.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `ValueError`: Invalid coordinate transformation parameters.
""" % __.ORIGIN()

sip_pix2foc = """
sip_pix2foc(*pixcrd, origin*) -> double array[ncoord][nelem]

Convert pixel coordinates to focal plane coordinates using the `SIP`_
polynomial distortion convention.

- *pixcrd*: double array[ncoord][nelem].  Array of pixel coordinates.

%s

Returns an array of focal plane coordinates.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `ValueError`: Invalid coordinate transformation parameters.
""" % __.ORIGIN()

spcfix = """
spcfix() -> int

Translates AIPS-convention spectral coordinate types.  {``FREQ``,
``VELO``, ``FELO``}-{``OBS``, ``HEL``, ``LSR``} (e.g. ``FREQ-LSR``,
``VELO-OBS``, ``FELO-HEL``)

Returns ``0`` for success; ``-1`` if no change required.
"""

spec = """
``int`` (read-only)

The index containing the spectral axis values.
"""

specsys = """
``string``

Spectral reference frame (standard of rest), ``SPECSYSa``.

.. seealso::

   `~pywcs.Wcsprm.ssysobs`, `~pywcs.Wcsprm.velosys`.
"""

sptr = """
sptr(*ctype, i=-1*)

Translates the spectral axis in a WCS object.  For example, a ``FREQ``
axis may be translated into ``ZOPT-F2W`` and vice versa.

- *ctype*: string.  Required spectral ``CTYPEia``, maximum of 8
  characters.  The first four characters are required to be given and
  are never modified.  The remaining four, the algorithm code, are
  completely determined by, and must be consistent with, the first
  four characters.  Wildcarding may be used, i.e.  if the final three
  characters are specified as ``"???"``, or if just the eighth
  character is specified as ``"?"``, the correct algorithm code will
  be substituted and returned.

- *i*: int.  Index of the spectral axis (0-relative).  If ``i < 0`` (or not
  provided), it will be set to the first spectral axis identified
  from the ``CTYPE`` keyvalues in the FITS header.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `SingularMatrixError`: Linear transformation matrix is singular.

- `InconsistentAxisTypesError`: Inconsistent or unrecognized
  coordinate axis types.

- `ValueError`: Invalid parameter value.

- `InvalidTransformError`: Invalid coordinate transformation
  parameters.

- `InvalidTransformError`: Ill-conditioned coordinate transformation
  parameters.

- `InvalidSubimageSpecificationError`: Invalid subimage specification
  (no spectral axis).
"""

ssysobs = """
``string``

The actual spectral reference frame in which there is no differential
variation in the spectral coordinate across the field-of-view,
``SSYSOBSa``.

.. seealso::

   `~pywcs.Wcsprm.specsys`, `~pywcs.Wcsprm.velosys`
"""

ssyssrc = """
``string``

The spectral reference frame (standard of rest) in which the redshift
was measured, ``SSYSSRCa``.
"""

sub = """
sub(*axes*) -> `~pywcs.WCS` object

Extracts the coordinate description for a subimage from a `~pywcs.WCS`
object.

The world coordinate system of the subimage must be separable in the
sense that the world coordinates at any point in the subimage must
depend only on the pixel coordinates of the axes extracted.  In
practice, this means that the ``PCi_ja`` matrix of the original image
must not contain non-zero off-diagonal terms that associate any of the
subimage axes with any of the non-subimage axes.

- *axes*: int or a sequence.

  - If an int, include the first *N* axes in their original
    order.

  - If a sequence, may contain a combination of image axis numbers
    (1-relative) or special axis identifiers (see below).  Order is
    significant; ``axes[0]`` is the axis number of the input image that
    corresponds to the first axis in the subimage, etc.

  - If ``0``, ``[]`` or ``None``, do a deep copy.

Coordinate axes types may be specified using either strings or
special integer constants.  The available types are:

  - ``'longitude'`` / ``WCSSUB_LONGITUDE``: Celestial longitude

  - ``'latitude'`` / ``WCSSUB_LATITUDE``: Celestial latitude

  - ``'cubeface'`` / ``WCSSUB_CUBEFACE``: Quadcube ``CUBEFACE`` axis

  - ``'spectral'`` / ``WCSSUB_SPECTRAL``: Spectral axis

  - ``'stokes'`` / ``WCSSUB_STOKES``: Stokes axis

Returns a `~pywcs.WCS` object, which is a deep copy of the original
object.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `InvalidSubimageSpecificationError`: Invalid subimage specification
  (no spectral axis).

- `NonseparableSubimageCoordinateSystem`: Non-separable subimage
  coordinate system.

.. note::

  Combinations of subimage axes of particular types may be extracted
  in the same order as they occur in the input image by combining the
  integer constants.  For example::

    wcs.sub([WCSSUB_LONGITUDE | WCSSUB_LATITUDE | WCSSUB_SPECTRAL])

  would extract the longitude, latitude, and spectral axes in the same
  order as the input image.  If one of each were present, the
  resulting object would have three dimensions.

  For convenience, ``WCSSUB_CELESTIAL`` is defined as the combination
  ``WCSSUB_LONGITUDE | WCSSUB_LATITUDE | WCSSUB_CUBEFACE``.

  The codes may also be negated to extract all but the types
  specified, for example::

    wcs.sub([
      WCSSUB_LONGITUDE,
      WCSSUB_LATITUDE,
      WCSSUB_CUBEFACE,
      -(WCSSUB_SPECTRAL | WCSSUB_STOKES)])

  The last of these specifies all axis types other than spectral or
  Stokes.  Extraction is done in the order specified by `axes`, i.e. a
  longitude axis (if present) would be extracted first (via
  ``axes[0]``) and not subsequently (via ``axes[3]``).  Likewise for
  the latitude and cubeface axes in this example.

  The number of dimensions in the returned object may be less than or
  greater than the length of `axes`.  However, it will never exceed
  the number of axes in the input image.
"""

theta0 = """
``double``

The native longitude of the fiducial point, i.e. the point whose
celestial coordinates are given in ``ref[1:2]``.  If undefined (NaN)
the initialization routine, `~pywcs.Wcsprm.set`, will set this to a
projection-specific default.

.. seealso::

   `~pywcs.Wcsprm.phi0`
"""

to_header = """
to_header(*relax=False*) -> string

`to_header` translates a WCS object into a FITS header.

    - If the `~pywcs.Wcsprm.colnum` member is non-zero then a binary
      table image array header will be produced.

    - Otherwise, if the `~pywcs.Wcsprm.colax` member is set non-zero
      then a pixel list header will be produced.

    - Otherwise, a primary image or image extension header will be
      produced.

The output header will almost certainly differ from the input in a
number of respects:

    1. The output header only contains WCS-related keywords.  In
       particular, it does not contain syntactically-required keywords
       such as ``SIMPLE``, ``NAXIS``, ``BITPIX``, or ``END``.

    2. Deprecated (e.g. ``CROTAn``) or non-standard usage will be
       translated to standard (this is partially dependent on whether
       `fix` was applied).

    3. Quantities will be converted to the units used internally,
       basically SI with the addition of degrees.

    4. Floating-point quantities may be given to a different decimal
       precision.

    5. Elements of the ``PCi_j`` matrix will be written if and only if
       they differ from the unit matrix.  Thus, if the matrix is unity
       then no elements will be written.

    6. Additional keywords such as ``WCSAXES``, ``CUNITia``,
       ``LONPOLEa`` and ``LATPOLEa`` may appear.

    7. The original keycomments will be lost, although
       `~pywcs.Wcsprm.to_header` tries hard to write meaningful
       comments.

    8. Keyword order may be changed.

Keywords can be translated between the image array, binary table, and
pixel lists forms by manipulating the `~pywcs.Wcsprm.colnum` or
`~pywcs.Wcsprm.colax` members of the `~pywcs.Wcsprm.WCS` object.

- *relax*: Degree of permissiveness:

    - ``False``: Recognize only FITS keywords defined by the published
      WCS standard.

    - ``True``: Admit all recognized informal extensions of the WCS
      standard.

Returns a raw FITS header as a string.
"""

unitfix = """
unitfix(*translate_units=''*) -> int

Translates non-standard ``CUNITia`` keyvalues.  For example, ``DEG`` ->
``deg``, also stripping off unnecessary whitespace.

- *translate_units*: string.  Do potentially unsafe translations of
  non-standard unit strings.

  Although ``"S"`` is commonly used to represent seconds, its
  recognizes ``"S"`` formally as Siemens, however rarely that may be
  translation to ``"s"`` is potentially unsafe since the standard
  used.  The same applies to ``"H"`` for hours (Henry), and ``"D"``
  for days (Debye).

  This string controls what to do in such cases, and is
  case-insensitive.

  - If the string contains ``"s"``, translate ``"S"`` to ``"s"``.

  - If the string contains ``"h"``, translate ``"H"`` to ``"h"``.

  - If the string contains ``"d"``, translate ``"D"`` to ``"d"``.

  Thus ``''`` doesn't do any unsafe translations, whereas ``'shd'``
  does all of them.

Returns ``0`` for success; ``-1`` if no change required.
"""

velangl = """
``double``

The angle in degrees that should be used to decompose an observed
velocity into radial and transverse components.

An undefined value is represented by NaN.
"""

velosys = """
``double``

The relative radial velocity (m/s) between the observer and the
selected standard of rest in the direction of the celestial reference
coordinate, ``VELOSYSa``.

An undefined value is represented by NaN.

.. seealso::

   `specsys`, `ssysobs`
"""

wcs = """
A `~pywcs.Wcsprm` object to perform the basic `wcslib`_ WCS
tranformation.
"""

Wcs = """
Wcs(*sip, cpdis, wcsprm, det2im*)

Wcs objects amalgamate basic WCS (as provided by `wcslib`_), with
`SIP`_ and `Paper IV`_ distortion operations.

To perform all distortion corrections and WCS tranformation, use
`all_pix2sky`.

- *sip*: A `~pywcs.Sip` object or ``None``

- *cpdis*: A pair of `~pywcs.DistortionLookupTable` objects, or
  ``(None, None)``.

- *wcsprm*: A `~pywcs.Wcsprm` object

- *det2im*: A pair of `~pywcs.DistortionLookupTable` objects, or
   ``(None, None)``.
"""

Wcsprm = """
Wcsprm(*header=None, key=' ', relax=False, naxis=2, keysel=0,
       colsel=None*)

`~pywcs.Wcsprm` is a direct wrapper around `wcslib`_, and provides
access to the core WCS transformations that it supports.

The FITS header parsing enforces correct FITS "keyword = value" syntax
with regard to the equals sign occurring in columns 9 and 10.
However, it does recognize free-format character (NOST 100-2.0,
Sect. 5.2.1), integer (Sect. 5.2.3), and floating-point values
(Sect. 5.2.4) for all keywords.

- *header*: A PyFITS header object or a string containing the raw FITS
  header data or ``None``.  If ``None``, the object will be
  initialized to default values.

- *key*: The key referring to a particular WCS transform in the
  header.  This may be either ``' '`` or ``'A'``-``'Z'`` and
  corresponds to the ``"a"`` part of ``"CTYPEia"``.  (*key*
  may only be provided if *header* is also provided.)

- *relax*: Degree of permissiveness:

    - ``False``: Recognize only FITS keywords defined by the published
      WCS standard.

    - ``True``: Admit all recognized informal extensions of the WCS
      standard.

- *naxis*: The number of sky coordinates axes for the object.
  (*naxis* may only be provided if *header* is ``None``.)

- *keysel*: Vector of flag bits that may be used to restrict the
  keyword types considered:

     - ``WCSHDR_IMGHEAD``: Image header keywords.

     - ``WCSHDR_BIMGARR``: Binary table image array.

     - ``WCSHDR_PIXLIST``: Pixel list keywords.

   If zero, there is no restriction.

- *colsel*: A sequence of table column numbers used to restrict the
  keywords considered.  ``None`` indicates no restriction.

**Exceptions:**

- `MemoryError`: Memory allocation failed.

- `ValueError`: Invalid key.

- `KeyError`: Key not found in FITS header.
"""

zsource = """
``double``

The redshift, ``ZSOURCEa``, of the source.

An undefined value is represented by NaN.
"""
