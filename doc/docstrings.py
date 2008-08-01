# It gets to be really tedious to type long docstrings in ANSI C
# syntax (since multi-line strings literals are not valid).
# Therefore, the docstrings are written here in doc/docstrings.py,
# which are then converted by setup.py into docstrings.h, which is
# included by pywcs.c

alt = """
Character code for alternate coordinate descriptions.  For example,
the C{"a"} in keyword names such as C{CTYPEia}).  This is blank for
the primary coordinate description, or one of the 26 upper-case
letters, A-Z.

@type: string
"""

cd = """
C{CDi_ja} linear transformation matrix.

For historical compatibility, two alternate specifications of the
C{CDi_ja} and C{CROTAia} keywords are supported.  Although these may
not formally co-exist with C{PCi_ja}, the approach here is simply to
ignore them if given in conjunction with C{PCi_ja}.

L{has_pci_ja}, L{has_cdi_ja} and L{has_crotaia} can be used to
determine which of these alternatives are present in the header.

C{CDi_ja} and C{CROTAia} keywords, if found, are to be stored in the
L{cd} and L{crota} arrays which are dimensioned similarly to L{pc} and
L{cdelt}.

These alternate specifications of the linear transformation matrix are
translated immediately to C{PCi_ja} by L{set} and are nowhere visible
to the lower-level routines.  In particular, L{set} resets L{cdelt} to
unity if C{CDi_ja} is present (and no C{PCi_ja}).  If no C{CROTAia} is
associated with the latitude axis, L{set} reverts to a unity C{PCi_ja}
matrix.

@type: array[2][2] of double
"""

cdelt = """
Coordinate increments (C{CDELTia}) for each coord axis.

An undefined value is represented by NaN.

@type: array[naxis] of double
"""

cel_offset = """
If true, an offset will be applied to C{(x, y)} to force C{(x,y) = (0,0)}
at the fiducial point.

@type: boolean
"""

celfix = """
celfix() -> int

Translates AIPS-convention celestial projection types, C{-NCP} and
C{-GLS}.

@return: C{0} for success; C{-1} if no change required.
"""

cname = """
A list of the coordinate axis names, C{CNAMEia}.

@type: list of strings
"""

colax = """
An array recording the column numbers for each axis in a pixel
list.

@type: array[naxis] of int
"""

colnum = """
Where the coordinate representation is associated with an image-array
column in a FITS binary table, this variable may be used to record the
relevant column number.

It should be set to zero for an image header or pixel list.
"""

copy = """
copy()

Creates a deep copy of the WCS object.
"""

cpdis = """
The pre-linear transformation distortion lookup tables, C{CPDIS}.

This is a 2-tuple of the form (x, y)
"""

cqdis = """
The post-linear transformation distortion lookup tables, C{CQDIS}.

This is a 2-tuple of the form (x, y)
"""

crder = """
The random error in each coordinate axis, C{CRDERia}.

An undefined value is represented by NaN.

@see: L{csyer}
@type: array[naxis] of double
"""

crota = """
C{CROTAia} keyvalues for each coordinate axis.

C{CROTAia} is an alternate specification of the linear transformation
matrix, maintained for historical compatibility.

@see: L{cd} for more information.

@type: array[naxis] of double
"""

_crpix_generic = """
Coordinate reference pixels (C{CRPIXja}) for each pixel axis.

B{Note that this coordinate is 1-based, and will be set/returned
exactly as in the header.}

%s

@type: array[naxis] of double
"""

crpix = _crpix_generic % """
This means that L{p2s}(crpix) will not be correct, and one should use
L{p2s_fits}(crpix) instead."""

crval = """
Coordinate reference values (C{CRVALia}) for each coordinate axis.

@type: array[naxis] of double
"""

csyer = """
The systematic error in the coordinate value axes, C{CSYERia}.

An undefined value is represented by NaN.

@see: L{crder}
@type: array[naxis] of double
"""

ctype = """
List of C{CTYPEia} keyvalues.

The L{ctype} keyword values must be in upper case and there must
be zero or one pair of matched celestial axis types, and zero or one
spectral axis.

@type: list of strings
"""

cubeface = """
Index into the C{pixcrd} (pixel coordinate) array for the C{CUBEFACE}
axis.  This is used for quadcube projections where the cube faces are
stored on a separate axis.

The quadcube projections (C{TSC}, C{CSC}, C{QSC}) may be represented
in FITS in either of two ways:

    - The six faces may be laid out in one plane and numbered as
      follows::


                                       0

                              4  3  2  1  4  3  2

                                       5

      Faces 2, 3 and 4 may appear on one side or the other (or both).
      The world-to-pixel routines map faces 2, 3 and 4 to the left but
      the pixel-to-world routines accept them on either side.

    - The C{"COBE"} convention in which the six faces are stored in a
      three-dimensional structure using a C{"CUBEFACE"} axis indexed
      from 0 to 5 as above.

These routines support both methods; L{set} determines which is being
used by the presence or absence of a C{CUBEFACE} axis in L{ctype}.
L{p2s} and L{s2p} translate the C{CUBEFACE} axis representation to the
single plane representation understood by the lower-level projection
routines.

@type: int
"""

cunit = """
List of C{CUNITia} keyvalues which define the units of measurement of
the C{CRVALia}, C{CDELTia} and C{CDi_ja} keywords.

As C{CUNITia} is an optional header keyword, L{cunit} may be left
blank but otherwise is expected to contain a standard units
specification as defined by WCS Paper I.  Utility function
C{wcsutrn()}, (not currently wrapped for Python) is available to
translate commonly used non-standard units specifications but this
must be done as a separate step before invoking L{set}.

For celestial axes, if L{cunit} is not blank, L{set} uses C{wcsunits}
to parse it and scale L{cdelt}, L{crval}, and L{cd} to decimal
degrees.  It then resets L{cunit} to "deg".

For spectral axes, if L{cunit} is not blank, L{set} uses C{wcsunits}
to parse it and scale L{cdelt}, L{crval}, and L{cd} to SI units.  It
then resets L{cunit} accordingly.

L{set} ignores L{cunit} for other coordinate types; L{cunit} may be
used to label coordinate values.

@type: list of strings
"""

cylfix = """
cylfix() -> int

Fixes WCS keyvalues for malformed cylindrical projections.

@return: C{0} for success; C{-1} if no change required.
"""

data = """
The array data for the distortion lookup table.  Should be a 2D numpy
double array.
"""

dateavg = """
Representative mid-point of the date of observation in ISO format,
C{yyyy-mm-ddThh:mm:ss}.

@see: L{dateobs}
@type: string
"""

dateobs = """
Start of the date of observation in ISO format,
C{yyyy-mm-ddThh:mm:ss}.

@see: L{dateavg}
@type: string
"""

datfix = """
datfix() -> int

Translates the old C{DATE-OBS} date format to year-2000 standard form
C{(yyyy-mm-ddThh:mm:ss)} and derives C{MJD-OBS} from it if not already
set.  Alternatively, if C{mjdobs} is set and C{dateobs} isn't, then
L{datfix} derives C{dateobs} from it.  If both are set but disagree by
more than half a day then C{ValueError} is raised.

@return: C{0} for success; C{-1} if no change required.
"""

Distortion = """
A class to perform the transformations and distortions outlined in
Paper IV.

Currently, only one direction -- from pixel to world coordinates -- is
supported.

Additionally, only lookup table distortions ("-TAB") are supported,
not polynomial or spline distortions.
"""

DistortionLookupTable = """
Represents a single lookup table for a distortion (Paper IV)
transformation.

These objects are used for setting the cpdis and cqdis lookup tables
in a Distortion object.
"""

distortion_cd = """
C{CDi_ja} linear transformation matrix.

For historical compatibility, either a single C{CDi_ja} matrix is
supported or a C{PCi_ja} matrix and C{CDELTia} keywords.

L{has_pc} can be used to determine which of these alternatives are
present in the header.

Setting C{cd} resets L{cdelt} to unity.

@type: array[2][2] of double
"""

distortion_crpix = _crpix_generic % ""

distortion_has_pc = """
has_pc() -> bool

Returns True if the distortion is using PC/CDELT, otherwise it is using CD.
"""

distortion_p2s = """
p2s(pixcrd) -> world

Converts pixel to world coordinates using the transformations and
distortions outlined in Paper IV.

B{The pixel coordinates given are 0-based (like array indices in C and
Python).  If your pixel coordinates are 1-based (like array indices in
Fortran), use L{p2s_fits} instead.}

@param pixcrd: Array of pixel coordinates.

@type pixcrd: numpy array[ncoord][nelem] of double

@return: array[ncoord][nelem] of double

Array of world coordinates.

@raises MemoryError: Memory allocation failed.
@raises SingularMatrixError: Linear transformation matrix is singular.
@raises InconsistentAxisTypesError: Inconsistent or unrecognized
    coordinate axis types.
@raises ValueError: Invalid parameter value.
@raises ValueError: x- and y-coordinate arrays are not the same size.
@raises InvalidTransformError: Invalid coordinate transformation
    parameters.
@raises InvalidTransformError: Ill-conditioned coordinate transformation
    parameters.
"""

distortion_p2s_fits = """
p2s_fits(pixcrd) -> world

Identical to L{p2s}, except pixel coordinates are 1-based (like array
indices in Fortran), instead of 0-based (like array indices C and
Python).
"""

distortion_pixel2world = """
An alias for L{p2s}.
"""

equinox = """
The equinox associated with dynamical equatorial or ecliptic
coordinate systems, C{EQUINOXa} (or C{EPOCH} in older headers).  Not
applicable to ICRS equatorial or ecliptic coordinates.

An undefined value is represented by NaN.

@type: float
"""

fix = """
fix(translate_units='', naxis=0) -> dict

Applies all of the corrections handled separately by L{datfix},
L{unitfix}, L{celfix}, L{spcfix} and L{cylfix}.

@param translate_units: Do potentially unsafe translations of
    non-standard unit strings.

    Although C{"S"} is commonly used to represent seconds, its
    translation to C{"s"} is potentially unsafe since the standard
    recognizes C{"S"} formally as Siemens, however rarely that may be
    used.  The same applies to C{"H"} for hours (Henry), and C{"D"}
    for days (Debye).

    This string controls what to do in such cases, and is
    case-insensitive.

        - If the string contains C{"s"}, translate C{"S"} to C{"s"}.
        - If the string contains C{"h"}, translate C{"H"} to C{"h"}.
        - If the string contains C{"d"}, translate C{"D"} to C{"d"}.

    Thus C{''} doesn't do any unsafe translations, whereas C{'shd'}
    does all of them.

@type translate_units: string

@param naxis: Image axis lengths.  If this array pointer is set to
    zero, then L{cylfix} will not be invoked.

@type naxis: array[naxis] of int

@return: A dictionary containing the following keys, each referring to
    a status string for each of the sub-fix functions that were
    called: L{datfix}, L{unitfix}, L{celfix}, L{spcfix}, L{cylfix}.
"""

get_offset = """
get_offset(x, y) -> offset

Returns the offset from the distortion table for pixel point (x, y).
"""

get_ps = """
get_ps() -> list of tuples

Returns C{PSi_ma} keywords for each I{i} and I{m}.  Returned as a list
of tuples of the form (I{i}, I{m}, I{value}):

    - I{i}: axis number, as in C{PSi_ma}, (i.e. 1-relative)
    - I{m}: parameter number, as in C{PSi_ma}, (i.e. 0-relative)
    - I{value}: parameter value (as a string)
"""

get_pv = """
get_pv() -> list of tuples

Returns C{PVi_ma} keywords for each I{i} and I{m}.  Returned as a list
of tuples of the form (I{i}, I{m}, I{value}):

    - I{i}: axis number, as in C{PVi_ma}, (i.e. 1-relative)
    - I{m}: parameter number, as in C{PVi_ma}, (i.e. 0-relative)
    - I{value}: parameter value (as a string)
"""

has_cdi_ja = """
has_cdi_ja() -> bool

Returns C{True} if C{CDi_ja} is present.  C{CDi_ja} is an alternate
specification of the linear transformation matrix, maintained for
historical compatibility.

Matrix elements in the IRAF convention are equivalent to the product
C{CDi_ja = CDELTia * PCi_ja}, but the defaults differ from that of the
C{PCi_ja} matrix.  If one or more C{CDi_ja} keywords are present then
all unspecified C{CDi_ja} default to zero.  If no C{CDi_ja} (or
C{CROTAia}) keywords are present, then the header is assumed to be in
C{PCi_ja} form whether or not any C{PCi_ja} keywords are present since
this results in an interpretation of C{CDELTia} consistent with the
original FITS specification.

While C{CDi_ja} may not formally co-exist with C{PCi_ja}, it may
co-exist with C{CDELTia} and C{CROTAia} which are to be ignored.

@see: L{cd} for more information.
"""

has_crotaia = """
has_crotaia() -> bool

Returns True if C{CROTAia} is present.  C{CROTAia} is an alternate
specification of the linear transformation matrix, maintained for
historical compatibility.

In the AIPS convention, C{CROTAia} may only be associated with the
latitude axis of a celestial axis pair.  It specifies a rotation in
the image plane that is applied AFTER the C{CDELTia}; any other
C{CROTAia} keywords are ignored.

C{CROTAia} may not formally co-exist with C{PCi_ja}.  C{CROTAia} and
C{CDELTia} may formally co-exist with C{CDi_ja} but if so are to be
ignored.

@see: L{cd} for more information.
"""

has_pci_ja = """
has_pci_ja() -> bool

Returns True if C{PCi_ja} is present.  C{PCi_ja} is the recommended way
to specify the linear transformation matrix.

@see: L{cd} for more information.
"""

imgpix_matrix = """
Inverse of the matrix containing the product of the C{CDELTia}
diagonal matrix and the C{PCi_ja} matrix.

I{This value may not be correct until after L{set} is called.}

@type: array[2][2] of double
"""

is_unity = """
is_unity() -> bool

Returns True if the linear transformation matrix (L{cd}) is unity.

I{This value may not be correct until after L{set} is called.}

@see: L{cd}
"""

lat = """
The index into the world coordinate array containing latitude values.
B{[Read only]}.

@see: L{lng}
@type: int
"""

latpole = """
The native latitude of the celestial pole, C{LATPOLEa} (deg).

@see: L{lonpole}
@type: float
"""

lattyp = """
Celestial axis type for latitude, e.g. RA.  B{[Read only]}.

I{This value may not be correct until after L{set} is called.}

@see: L{ctype}
@type: string
"""

lng = """
The index into the world coordinate array containing longitude values.
B{[Read only]}.

@see: L{lat}
@type: int
"""

lngtyp = """
Celestial axis type for longitude, e.g. DEC.  B{[Read only]}.

I{This value may not be correct until after L{set} is called.}

@see: L{ctype}
@type: string
"""

lonpole = """
The native longitude of the celestial pole, C{LONPOLEa} (deg).

@see: L{latpole}
@type: float
"""

mix = """
mix(mixpix, mixcel, vspan, vstep, viter, world, pixcrd) -> dict

Given either the celestial longitude or latitude plus an element of
the pixel coordinate, solves for the remaining elements by iterating
on the unknown celestial coordinate element using L{s2p}.

B{The pixel coordinates given and returned are 0-based (like array
indices in C and Python).  If your pixel coordinates are 1-based (like
array indices in Fortran), use L{mix_fits} instead.}

@param mixpix: Which element on the pixel coordinate is given.
@type mixpix: int

@param mixcel: Which element of the celestial coordinate is given. If
    C{1}, celestial longitude is given in C{world[self.L{lng}]},
    latitude returned in C{world[self.L{lat}]}.  If C{2}, celestial
    latitude is given in C{world[self.L{lat}]}, longitude returned in
    C{world[self.L{lng}]}

@type mixcel: int

@param vspan: Solution interval for the celestial coordinate, in
    degrees.  The ordering of the two limits is irrelevant.  Longitude
    ranges may be specified with any convenient normalization, for
    example C{(-120,+120)} is the same as C{(240,480)}, except that
    the solution will be returned with the same normalization,
    i.e. lie within the interval specified.

@type vspan: sequence of two floats

@param vstep: Step size for solution search, in degrees.  If C{0}, a
    sensible, although perhaps non-optimal default will be used.

@type vstep: float

@param viter: If a solution is not found then the step size will be
    halved and the search recommenced.  C{viter} controls how many
    times the step size is halved.  The allowed range is 5 - 10.

@type viter: int

@param world: World coordinate elements.  C{world[self.L{lng}]} and
    C{world[self.L{lat}]} are the celestial longitude and latitude, in
    degrees.  Which is given and which returned depends on the value
    of C{mixcel}.  All other elements are given.  The results will be
    written to this array in-place.

@type world: array[naxis] of double

@param pixcrd: Pixel coordinate.  The element indicated by C{mixpix}
    is given and the remaining elements will be written in-place.
    B{The pixel coordinates given and returned are zero-based.  If
    your pixel coordinates are one-based, use L{mix_fits} instead.}

@type pixcrd: array[naxis] of double

@return: A dictionary with the following keys:

    - C{phi} (array[naxis] of double)
    - C{theta} (array[naxis] of double)

        - Longitude and latitude in the native coordinate system of the
          projection, in degrees.

    - C{imgcrd} (array[naxis] of double)

        - Image coordinate elements.  C{imgcrd[self.L{lng}]} and
          C{imgcrd[self.L{lat}]} are the projected I{x}- and
          I{y}-coordinates, in decimal degrees.

@raise MemoryError: Memory allocation failed.
@raise SingularMatrixError: Linear transformation matrix is singular.
@raise InconsistentAxisTypesError: Inconsistent or unrecognized coordinate
    axis types.
@raise ValueError: Invalid parameter value.
@raise InvalidTransformError: Invalid coordinate transformation parameters.
@raise InvalidTransformError: Ill-conditioned coordinate
transformation parameters.
@raise InvalidCoordinateError: Invalid world coordinate.
@raise NoSolutionError: No solution found in the specified interval.
"""

mix_fits = """
mix_fits(mixpix, mixcel, vspan, vstep, viter, world, pixcrd) -> dict

Identical to L{mix}, except pixel coordinates are 1-based (like array
indices in Fortran), instead of 0-based (like array indices C and
Python).
"""

mjdavg = """
Modified Julian Date (MJD = JD - 2400000.5), C{MJD-AVG}, corresponding
to C{DATE-AVG}.

An undefined value is represented by NaN.

@see: L{mjdobs}
@type: float
"""

mjdobs = """
Modified Julian Date (MJD = JD - 2400000.5), C{MJD-OBS}, corresponding
to C{DATE-OBS}.

An undefined value is represented by NaN.

@see: L{mjdavg}
@type: float
"""

name = """
The name given to the coordinate representation C{WCSNAMEa}.

@type: string
"""

naxis = """
The number of axes (pixel and coordinate), given by the C{NAXIS} or
C{WCSAXESa} keyvalues.  B{[Read only]}.

The number of coordinate axes is determined at parsing time, and can
not be subsequently changed.

It is determined from the highest of the following:

  1. C{NAXIS}

  2. C{WCSAXESa}

  3. The highest axis number in any parameterized WCS keyword.  The
     keyvalue, as well as the keyword, must be syntactically valid
     otherwise it will not be considered.

If none of these keyword types is present, i.e. if the header only
contains auxiliary WCS keywords for a particular coordinate
representation, then no coordinate description is constructed for it.

This value may differ for different coordinate representations of the
same image.

@type: int
"""

obsgeo = """
Location of the observer in a standard terrestrial reference frame,
C{OBSGEO-X}, C{OBSGEO-Y}, C{OBSGEO-Z} (in metres).

An undefined value is represented by NaN.

@type: array[3] of double
"""

p2s = """
p2s(pixcrd) -> dict

Converts pixel to world coordinates.

B{The pixel coordinates given are 0-based (like array indices in C and
Python).  If your pixel coordinates are 1-based (like array indices in
Fortran), use L{p2s_fits} instead.}

@param pixcrd: Array of pixel coordinates.

@type pixcrd: numpy array[ncoord][nelem] of double

@return: A dictionary with the following keys:

        - C{imgcrd} (array[ncoord][nelem] of double)

            - Array of intermediate world coordinates.  For celestial
              axes, C{imgcrd[][self.L{lng}]} and
              C{imgcrd[][self.L{lat}]} are the projected I{x}-, and
              I{y}-coordinates, in decimal degrees.  For spectral
              axes, C{imgcrd[][self.L{spec}]} is the intermediate
              spectral coordinate, in SI units.

        - C{phi} (array[ncoord] of double)

        - C{theta} (array[ncoord] of double)

            - Longitude and latitude in the native coordinate system
              of the projection, in degrees.

        - C{world} (array[ncoord][nelem] of double)

            - Array of world coordinates.  For celestial axes,
              C{world[][self.L{lng}]} and C{world[][self.L{lat}]} are
              the celestial longitude and latitude, in decimal
              degrees.  For spectral axes, C{world[][self.L{spec}]} is
              the intermediate spectral coordinate, in SI units.

        - C{stat} (array[ncoord] of int)

            - Status return value for each coordinate. C{0} for success, C{1}
              for invalid pixel coordinate.

@raises MemoryError: Memory allocation failed.
@raises SingularMatrixError: Linear transformation matrix is singular.
@raises InconsistentAxisTypesError: Inconsistent or unrecognized
    coordinate axis types.
@raises ValueError: Invalid parameter value.
@raises ValueError: x- and y-coordinate arrays are not the same size.
@raises InvalidTransformError: Invalid coordinate transformation
    parameters.
@raises InvalidTransformError: Ill-conditioned coordinate transformation
    parameters.
"""

p2s_fits = """
p2s_fits(pixcrd) -> dict

Identical to L{p2s}, except pixel coordinates are 1-based (like array
indices in Fortran), instead of 0-based (like array indices C and
Python).
"""

pc = """
The C{PCi_ja} (pixel coordinate) transformation matrix.  The order is::

  [[PC1_1, PC1_2],
   [PC2_1, PC2_2]]

@type: array[2][2] of double
"""

phi0 = """
The native latitude of the fiducial point, i.e. the point whose
celestial coordinates are given in ref[1:2].  If undefined (NaN) the
initialization routine, L{set}, will set this to a projection-specific
default.

@see: L{theta0}
@type: float
"""

piximg_matrix = """
Matrix containing the product of the C{CDELTia} diagonal matrix and
the C{PCi_ja} matrix.

I{This value may not be correct until after L{set} is called.}

@type: array[2][2] of double
"""

print_contents = """
print_contents()

Print the contents of the WCS object to stdout.  Probably only useful
for debugging purposes, and may be removed in the future.
"""

pywcs = """
The routines in this module implement the FITS World Coordinate System
(WCS) standard which defines methods to be used for computing world
coordinates from image pixel coordinates, and vice versa.

The basic workflow is as follows:

    1. C{import pywcs}

    2. Call the C{pywcs.WCS} constructor with a PyFITS header object.

    3. Optionally, if the FITS file uses any deprecated or
    non-standard features, you may need to call one of the C{fix}
    methods on the object.

    4. Convert coordinates using the C{pixel2world()} or
    C{world2pixel()} methods.

Short example::

    import numpy
    import pywcs
    import pyfits

    hdulist = pyfits.open("test.fits")

    # Parse the WCS keywords in the primary HDU
    wcs = pywcs.WCS(hdulist[0].header)

    # Print out the "name" of the WCS, as defined in the FITS header
    print wcs.name

    # Some interesting pixel coordinates
    pixcrd = numpy.array([[0,0],[24,38],[45,98]], numpy.float_)

    # Convert pixel coordinates to world coordinates
    world = wcs.pixel2world(pixcrd)
    print world
"""

radesys = """
The equatorial or ecliptic coordinate system type, C{RADESYSa}.

@type: string
"""

restfrq = """
Rest frequency (Hz) from C{RESTFRQa}.

An undefined value is represented by NaN.

@see: L{restwav}
@type: float
"""

restwav = """
Rest wavelength (m) from C{RESTWAVa}.

An undefined value is represented by NaN.

@see: L{restfrq}
@type: float
"""

s2p = """
s2p(world) -> dict

Transforms world coordinates to pixel coordinates.

B{The pixel coordinates returned are 0-based (like array indices in C
and Python).  If you require pixel coordinates to be 1-based (like
array indices in Fortran), use L{s2p_fits} instead.}

@param world: Array of world coordinates, in decimal degrees.

@type world: array[ncoord][nelem] of double

@return: A dictionary with the following keys:
        - C{phi} (array[ncoord] of double)
        - C{theta} (array[ncoord] of double)

            - Longitude and latitude in the native coordinate system
              of the projection, in degrees.

        - C{imgcrd} (array[ncoord][nelem] of double)

            - Array of intermediate world coordinates.  For celestial
              axes, C{imgcrd[][self.L{lng}]} and
              C{imgcrd[][self.L{lat}]} are the projected I{x}-, and
              I{y}-coordinates, in "degrees".  For quadcube
              projections with a C{CUBEFACE} axis, the face number is
              also returned in C{imgcrd[][self.L{cubeface}]}.  For
              spectral axes, C{imgcrd[][self.L{spec}]} is the
              intermediate spectral coordinate, in SI units.

        - C{pixcrd} (array[ncoord][nelem] of double)

            - Array of pixel coordinates.  B{Pixel coordinates are
              zero-based.}

        - C{stat} (array[ncoord] of int)

            - Status return value for each coordinate. C{0} for
              success, C{1} for invalid pixel coordinate.

@raises MemoryError: Memory allocation failed.
@raises SingularMatrixError: Linear transformation matrix is singular.
@raises InconsistentAxisTypesError: Inconsistent or unrecognized
    coordinate axis types.
@raises ValueError: Invalid parameter value.
@raises InvalidTransformError: Invalid coordinate transformation
    parameters.
@raises InvalidTransformError: Ill-conditioned coordinate
    transformation parameters.
"""

s2p_fits = """
s2p_fits(pixcrd) -> dict

Identical to L{s2p}, except pixel coordinates are 1-based (like array
indices in Fortran), instead of 0-based (like array indices C and
Python).
"""

set = """
set()

Sets up a WCS object for use according to information supplied within
it.

Note that this routine need not be called directly; it will be invoked
by L{p2s} and L{s2p} if necessary.

Some attributes that are based on other attributes (such as L{lattyp} on
L{ctype}) may not be correct until after L{set} is called.

C{set} strips off trailing blanks in all string members.

Among all the obvious details, C{set} recognizes the C{NCP} projection
and converts it to the equivalent C{SIN} projection and it also
recognizes C{GLS} as a synonym for C{SFL}.  It does alias translation
for the AIPS spectral types (C{FREQ-LSR}, C{FELO-HEL}, etc.) but
without changing the input header keywords.

@raises MemoryError: Memory allocation failed.
@raises SingularMatrixError: Linear transformation matrix is singular.
@raises InconsistentAxisTypesError: Inconsistent or unrecognized
    coordinate axis types.
@raises ValueError: Invalid parameter value.
@raises InvalidTransformError: Invalid coordinate transformation
    parameters.
@raises InvalidTransformError: Ill-conditioned coordinate transformation
    parameters.
"""

set_ps = """
set_ps(list)

Sets C{PSi_ma} keywords for each I{i} and I{m}.  The input must be a
sequence of tuples of the form (I{i}, I{m}, I{value}):

    - I{i}: axis number, as in C{PSi_ma}, (i.e. 1-relative)
    - I{m}: parameter number, as in C{PSi_ma}, (i.e. 0-relative)
    - I{value}: parameter value (as a string)
"""

set_pv = """
set_pv(list)

Sets C{PVi_ma} keywords for each I{i} and I{m}.  The input must be a
sequence of tuples of the form (I{i}, I{m}, I{value}):

    - I{i}: axis number, as in C{PVi_ma}, (i.e. 1-relative)
    - I{m}: parameter number, as in C{PVi_ma}, (i.e. 0-relative)
    - I{value}: parameter value (as a string)
"""

spcfix = """
spcfix() -> int

Translates AIPS-convention spectral coordinate types.
{C{FREQ},C{VELO},C{FELO}}-{C{OBS},C{HEL},C{LSR}} (e.g. C{FREQ-LSR},
C{VELO-OBS}, C{FELO-HEL})

@return: C{0} for success; C{-1} if no change required.
"""

spec = """
The index containing the spectral axis values.  B{[Read only]}.

@type: int
"""

specsys = """
Spectral reference frame (standard of rest), C{SPECSYSa}

@see: L{ssysobs}, L{velosys}.
@type: string
"""

sptr = """
sptr(ctype, i=-1)

Translates the spectral axis in a WCS object.  For example, a C{FREQ}
axis may be translated into C{ZOPT-F2W} and vice versa.

@param ctype: Required spectral C{CTYPEia}.  Wildcarding may be used,
    i.e.  if the final three characters are specified as C{"???"}, or
    if just the eighth character is specified as C{"?"}, the correct
    algorithm code will be substituted and returned.
@type ctype: string

@param i: Index of the spectral axis (0-rel).  If C{i < 0} (or not
    provided), it will be set to the first spectral axis identified
    from the C{CTYPE} keyvalues in the FITS header.
@type i: int

@raises MemoryError: Memory allocation failed.
@raises SingularMatrixError: Linear transformation matrix is singular.
@raises InconsistentAxisTypesError: Inconsistent or unrecognized
    coordinate axis types.
@raises ValueError: Invalid parameter value.
@raises InvalidTransformError: Invalid coordinate transformation
    parameters.
@raises InvalidTransformError: Ill-conditioned coordinate
    transformation parameters.
@raises InvalidSubimageSpecificationError: Invalid subimage
    specification (no spectral axis).
"""

ssysobs = """
The actual spectral reference frame in which there is no differential
variation in the spectral coordinate across the field-of-view,
C{SSYSOBSa}.

@see: L{specsys}, L{velosys}
@type: string
"""

ssyssrc = """
The spectral reference frame (standard of rest) in which the redshift
was measured, C{SSYSSRCa}.

@see: L{zsource}
@type: string
"""

theta0 = """
The native longitude of the fiducial point, i.e. the point whose
celestial coordinates are given in ref[1:2].  If undefined (NaN) the
initialization routine, L{set}, will set this to a projection-specific
default.

@see: L{phi0}
@type: float
"""

to_header = """
to_header(relax=False) -> string

L{to_header} translates a WCS object into a FITS header.

    - If the L{colnum} member is non-zero then a binary table image
      array header will be produced.

    - Otherwise, if the L{colax} member is set non-zero then a pixel
      list header will be produced.

    - Otherwise, a primary image or image extension header will be
      produced.

The output header will almost certainly differ from the input in a
number of respects:

    1. The output header only contains WCS-related keywords.  In
       particular, it does not contain syntactically-required keywords
       such as C{SIMPLE}, C{NAXIS}, C{BITPIX}, or C{END}.

    2. Deprecated (e.g. C{CROTAn}) or non-standard usage will be
       translated to standard (this is partially dependent on whether
       L{fix} was applied).

    3. Quantities will be converted to the units used internally,
       basically SI with the addition of degrees.

    4. Floating-point quantities may be given to a different decimal
       precision.

    5. Elements of the C{PCi_j} matrix will be written if and only if
       they differ from the unit matrix.  Thus, if the matrix is unity
       then no elements will be written.

    6. Additional keywords such as C{WCSAXES}, C{CUNITia}, C{LONPOLEa}
       and C{LATPOLEa} may appear.

    7. The original keycomments will be lost, although L{to_header}
       tries hard to write meaningful comments.

    8. Keyword order may be changed.

Keywords can be translated between the image array, binary table, and
pixel lists forms by manipulating the L{colnum} or L{colax} members of
the WCS object.

@param relax: Degree of permissiveness:
    - C{False}: Recognize only FITS keywords defined by the
      published WCS standard.
    - C{True}: Admit all recognized informal extensions of the
      WCS standard.
@type relax: bool

@return: A raw FITS header as a string.
"""

unitfix = """
unitfix(translate_units='') -> int

Translates non-standard C{CUNITia} keyvalues.  For example, C{DEG} ->
C{deg}, also stripping off unnecessary whitespace.

@param translate_units: Do potentially unsafe translations of
    non-standard unit strings.

    Although C{"S"} is commonly used to represent seconds, its
    recognizes C{"S"} formally as Siemens, however rarely that may be
    translation to C{"s"} is potentially unsafe since the standard
    used.  The same applies to C{"H"} for hours (Henry), and C{"D"}
    for days (Debye).

    This string controls what to do in such cases, and is
    case-insensitive.

        - If the string contains C{"s"}, translate C{"S"} to C{"s"}.
        - If the string contains C{"h"}, translate C{"H"} to C{"h"}.
        - If the string contains C{"d"}, translate C{"D"} to C{"d"}.

    Thus C{''} doesn't do any unsafe translations, whereas C{'shd'}
    does all of them.

@return: C{0} for success; C{-1} if no change required.
"""

velangl = """
The angle in degrees that should be used to decompose an observed
velocity into radial and transverse components.

An undefined value is represented by NaN.

@type: float
"""

velosys = """
The relative radial velocity (m/s) between the observer and the
selected standard of rest in the direction of the celestial reference
coordinate, C{VELOSYSa}.

An undefined value is represented by NaN.

@see: L{specsys}, L{ssysobs}
@type: float
"""

WCS = """
WCS(header=None, key=' ', relax=False, naxes=2)

WCS objects convert between pixel and world coordinates, based on
the WCS settings in a FITS file.

The FITS header parsing enforces correct FITS "keyword = value" syntax
with regard to C{"= "} occurring in columns 9 and 10.  However, it
does recognize free-format character (NOST 100-2.0, Sect. 5.2.1),
integer (Sect. 5.2.3), and floating-point values (Sect. 5.2.4) for all
keywords.

@param header: A PyFITS header object or a string containing the raw
    FITS header data.  If header is not provided, the object will be
    initialized to default values.
@type header: PyFITS header object or string

@param key: The key referring to a particular WCS transform in the
    header.  This may be either C{' '} or C{'A'}-C{'Z'} and
    corresponds to the C{"a"} part of C{"CTYPEia"}.  (C{key}
    may only be provided if C{header} is also provided.)
@type key: string

@param relax: Degree of permissiveness:
    - C{False}: Recognize only FITS keywords defined by the
      published WCS standard.
    - C{True}: Admit all recognized informal extensions of the
      WCS standard.
    (C{relax} may be provided only if C{header} is also provided.)
@type relax: bool

@param naxes: The number of world coordinates axes for the object.
    (C{naxes} may only be provided if C{header} is C{None}.)
@type naxes: int

@raises MemoryError: Memory allocation failed.
@raises ValueError: Invalid key.
@raises KeyError: Key not found in FITS header.
"""

zsource = """
The redshift, C{ZSOURCEa}, of the source.

An undefined value is represented by NaN.

@see: L{ssyssrc}
@type: float
"""
