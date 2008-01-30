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

linear transformation matrix are supported, those associated with the
For historical compatibility, two alternate specifications of the
C{CDi_ja} and C{CROTAia} keywords.  Although these may not formally
co-exist with C{PCi_ja}, the approach here is simply to ignore them if
given in conjunction with C{PCi_ja}.

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

Any undefined elements are returned as NaN.

@type: array[naxis] of double
"""

celfix = """
celfix() -> int

Translates AIPS-convention celestial projection types, C{-NCP}
and C{-GLS}.

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

Creates a deep copy of the Wcs object.
"""

crder = """
The random error in each coordinate axis, C{CRDERia}.

Any undefined elements are returned as NaN.

@see: L{csyer}
@type: array[naxis] of double
"""

crota = """
C{CROTAia} keyvalues for each coord axis.

C{CROTAia} is an alternate
specification of the linear transformation matrix, maintained for
historical compatibility.

@see: L{cd} for more information.

@type: array[2][2] of double
"""

crpix = """
Coordinate reference pixels (C{CRPIXja}) for each pixel axis.

@type: array[naxis] of double
"""

crval = """
Coordinate reference values (C{CRVALia}) for each coordinate axis.

@type: array[naxis] of double
"""

csyer = """
The systematic error in the coordinate value axes, C{CSYERia}.

Any undefined elements are returned as NaN.

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
to parse it and scale L{cdelt}, L{crval}, and L{cd} to degrees.  It
then resets L{cunit} to "deg".

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

equinox = """
The equinox associated with dynamical equatorial or ecliptic
coordinate systems, C{EQUINOXa} (or C{EPOCH} in older headers).  Not
applicable to ICRS equatorial or ecliptic coordinates.

If the equinox is undefined, NaN is returned.

@type: float
"""

fix = """fix(translate_units='', naxis=0) -> dict

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

lat = """
The index into the world coordinate array containing latitude values.

@see: L{lng}
@type: int
"""

latpole = """
The native latitude of the celestial pole, C{LATPOLEa} (deg).

@see: L{lonpole}
@type: float
"""

lng = """
The index into the world coordinate array containing longitude values.

@see: L{lat}
@type: int
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

@param mixpix: Which element on the pixel coordinate is given.
@type mixpix: int

@param mixcel: Which element of the celestial coordinate is given. If
    C{1}, celestial longitude is given in C{world[self.L{lng}]},
    latitude returned in C{world[self.L{lat}]}.  If C{2}, celestial
    latitude is given in world[self.lat], longitude returned in
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
    of mixcel.  All other elements are given.  The results will be
    written to this array in-place.

@type world: array[naxis] of double

@param pixcrd: Pixel coordinate.  The element indicated by mixpix is
    given and the remaining elements will be written in-place.

@type pixcrd: array[naxis] of double

@return: A dictionary with the following keys:

    - C{phi} (array[naxis] of double)
    - C{theta} (array[naxis] of double)

        - Longitude and latitude in the native coordinate system of the
          projection, in degrees.

    - C{imgcrd} (array[naxis] of double)

        - Image coordinate elements.  C{imgcrd[self.L{lng}]} and
          C{imgcrd[self.L{lat}]} are the projected I{x}- and
          I{y}-coordinates, in "degrees".

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

mjdavg = """
Modified Julian Date (MJD = JD - 2400000.5), C{MJD-AVG}, corresponding
to C{DATE-AVG}.

If mjdavg is undefined, NaN is returned.

@see: L{mjdobs}
@type: float
"""

mjdobs = """
Modified Julian Date (MJD = JD - 2400000.5), C{MJD-OBS}, corresponding
to C{DATE-OBS}.

If mjdobs is undefined, NaN is returned.

@see: L{mjdavg}
@type: float
"""

name = """
The name given to the coordinate representation C{WCSNAMEa}.

@type: string
"""

naxis = """
The number of axes (pixel and coordinate), given by the C{NAXIS} or
C{WCSAXESa} keyvalues.

@type: int
"""

obsgeo = """
Location of the observer in a standard terrestrial reference frame,
C{OBSGEO-X}, C{OBSGEO-Y}, C{OBSGEO-Z} (in metres).

Any undefined elements are returned as NaN.

@type: array[3] of double
"""

p2s = """
p2s(pixcrd) -> dict

Converts pixel to world coordinates.

@param pixcrd: Array of pixel coordinates.

@type pixcrd: numpy array[ncoord][nelem] of double

@return: A dictionary with the following keys:

        - C{imgcrd} (array[ncoord][nelem] of double)

            - Array of intermediate world coordinates.  For celestial
              axes, C{imgcrd[][self.L{lng}]} and
              C{imgcrd[][self.L{lat}]} are the projected I{x}-, and
              I{y}-coordinates, in "degrees".  For spectral axes,
              C{imgcrd[][self.L{spec}]} is the intermediate spectral
              coordinate, in SI units.

        - C{phi} (array[ncoord] of double)

        - C{theta} (array[ncoord] of double)

            - Longitude and latitude in the native coordinate system
              of the projection, in degrees.

        - C{world} (array[ncoord][nelem] of double)

            - Array of world coordinates.  For celestial axes,
              C{world[][self.L{lng}]} and C{world[][self.L{lat}]} are
              the celestial longitude and latitude, in degrees.  For
              spectral axes, C{world[][self.L{spec}]} is the
              intermediate spectral coordinate, in SI units.

        - C{stat} (array[ncoord] of int)

            - Status return value for each coordinate. C{0} for success, C{1}
              for invalid pixel coordinate.

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

parse_image_header = """
parse_image_header(header, relax=False) -> dict of L{Wcs} objects

Parses a FITS image header, either that of a primary HDU or of an
image extension.  All WCS keywords defined in Papers I, II, and III
are recognized, and also those used by the AIPS convention and certain
other keywords that existed in early drafts of the WCS papers.

Given a string containing a FITS image header, C{parse_image_header()}
identifies and reads all WCS keywords for the primary coordinate
representation and up to 26 alternate representations.  It returns
this information as a list of C{Wcs} objects.

C{parse_image_header} fills in information associated with coordinate
lookup tables.

C{parse_image_header} determines the number of coordinate axes
independently for each alternate coordinate representation (denoted by
the C{"a"} value in keywords like C{CTYPEia}) from the higher of

    - C{NAXIS}

    - C{WCSAXES}

    - The highest axis number in any parameterized WCS keyword.  The
      keyvalue, as well as the keyword, must be syntactically valid
      otherwise it will not be considered.

If none of these keyword types is present, i.e. if the header only
contains auxiliary WCS keywords for a particular coordinate
representation, then no coordinate description is constructed for it.

C{parse_image_header} enforces correct FITS "keyword = value" syntax
with regard to C{"= "} occurring in columns 9 and 10.  However, it
does recognize free-format character (NOST 100-2.0, Sect. 5.2.1),
integer (Sect. 5.2.3), and floating-point values (Sect. 5.2.4) for all
keywords.

@param header: String containing the (entire) FITS image header from
    which to identify and construct the coordinate representations.
@type header: string
@param relax: Degree of permissiveness:
    - C{False}: Recognize only FITS keywords defined by the
      published WCS standard.
    - C{True}: Admit all recognized informal extensions of the
      WCS standard.
@type relax: bool

@returns: A dictionary where the keys are single-character strings
   that identify a particular WCS transformation.  For example, the
   C{"a"} in keyword names such as C{CTYPEia}).  This is blank (C{' '})
   for the primary coordinate description, or one of the 26 upper-case
   letters, A-Z.

@raises MemoryError: Memory allocation failed.
"""
#    /* TODO: Deal with this next paragraph in a Pythonic way
#    For negative values of ctrl (see below), header[] is modified so
#    "that WCS keyrecords processed by wcspih() are removed from it.\n" */

pc = """
The C{PCi_ja} (pixel coordinate) transformation matrix.  The order is::

  [[PC1_1, PC1_2],
   [PC2_1, PC2_2]]

@type: array[2][2] of double
"""

print_contents = """
print_contents()

Print the contents of the Wcs object to stdout.  Probably only useful
for debugging purposes, and may be removed in the future.
"""

ps = """
C{PSi_ma} keywords for each I{i} and I{m}.  Returned as a list of
tuples of the form (I{i}, I{m}, I{value}):

    - I{i}: axis number, as in C{PSi_ma}, (i.e. 1-relative)
    - I{m}: parameter number, as in C{PSi_ma}, (i.e. 0-relative)
    - I{value}: parameter value (as a string)

@type: list of tuples
"""

pv = """
C{PVi_ma} keywords for each I{i} and I{m}.  Returned as a list of
tuples of the form (I{i}, I{m}, I{value}):

    - I{i}: axis number, as in C{PVi_ma}, (i.e. 1-relative)
    - I{m}: parameter number, as in C{PVi_ma}, (i.e. 0-relative)
    - I{value}: parameter value (as a string)

@type: list of tuples
"""

pywcs = """
The routines in this module implement the FITS World Coordinate System
(WCS) standard which defines methods to be used for computing world
coordinates from image pixel coordinates, and vice versa.

The basic workflow is as follows:

    1. C{import pywcs}

    2. Call C{pywcs.parse_image_header()} with a raw FITS header or
    C{pywcs.parse_hdu()} with a PyFITS HDU object.  This will return a
    dictionary containing all of the WCS transformations specified in
    the FITS file.

    3. Optionally, if the FITS file uses any deprecated or
    non-standard features, you may need to call one of the C{fix}
    methods on the object.

    4. Select the correct C{Wcs} transform object, and convert
    coordinates using the C{p2s()} or C{s2p()} methods.

Short example::

    import numpy
    import pywcs
    import pyfits

    hdulist = pyfits.open("test.fits")

    # Parse the WCS keywords in the primary HDU
    wcsdict = pywcs.parse_hdu(hdulist[0])

    # Get the "blank" or primary Wcs
    wcs = wcsdict[' ']

    # Print out the "name" of the Wcs, as defined in the FITS header
    print wcs.name

    # Some interesting pixel coordinates
    pixcrd = numpy.array([[0,0],[24,38],[45,98]], numpy.float_)

    # Convert pixel coordinates to world coordinates
    result = wcs.p2s(pixcrd)
    world = result['world']
    print world
"""

radesys = """
The equatorial or ecliptic coordinate system type, C{RADESYSa}.

@type: string
"""

restfrq = """
Rest frequency (Hz) from C{RESTFRQa}.

If restfrq is undefined, NaN is returned.

@see: L{restwav}
@type: float
"""

restwav = """
Rest wavelength (m) from C{RESTWAVa}.

If restwav is undefined, NaN is returned.

@see: L{restfrq}
@type: float
"""

s2p = """
s2p(world) -> dict

Transforms world coordinates to pixel coordinates.


@param world: Array of world coordinates.
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

            - Array of pixel coordinates.

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

set = """
set()

Sets up a Wcs object for use according to information supplied within
it.

C{set} recognizes the C{NCP} projection and converts it to the
equivalent C{SIN} projection and it also recognizes C{GLS} as a
synonym for C{SFL}.  It does alias translation for the AIPS spectral
types (C{FREQ-LSR}, C{FELO-HEL}, etc.) but without changing the input
header keywords.

Note that this routine need not be called directly; it will be invoked
by L{p2s} and L{s2p} if necessary.

C{set} strips off trailing blanks in all string members.

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

spcfix = """
spcfix() -> int

Translates AIPS-convention spectral coordinate types.
{C{FREQ},C{VELO},C{FELO}}-{C{OBS},C{HEL},C{LSR}} (e.g. C{FREQ-LSR},
C{VELO-OBS}, C{FELO-HEL})

@return: C{0} for success; C{-1} if no change required.
"""

spec = """
The index containing the spectral axis values.

@type: int
"""

specsys = """
Spectral reference frame (standard of rest), SPECSYSa

@see: L{ssysobs}, L{velosys}.
@type: string
"""

sptr = """
sptr(ctype, i=-1)

Translates the spectral axis in a Wcs object.  For example, a C{FREQ}
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

If velangl is undefined, NaN is returned.

@type: float
"""

velosys = """
The relative radial velocity (m/s) between the observer and the
selected standard of rest in the direction of the celestial reference
coordinate, C{VELOSYSa}.

If velosys is undefined, NaN is returned.

@see: L{specsys}, L{ssysobs}
@type: float
"""

Wcs = """
Wcs objects can convert between pixel and world coordinates, based on
the WCS settings in a FITS file.

To create Wcs objects, one would normally use
pywcs.parse_image_header.
"""

zsource = """
The redshift, C{ZSOURCEa}, of the source.

@see: L{ssyssrc}
@type: float
"""
