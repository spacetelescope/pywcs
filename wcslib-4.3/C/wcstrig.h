/*============================================================================

    WCSLIB 4.3 - an implementation of the FITS WCS standard.
    Copyright (C) 1995-2007, Mark Calabretta

    This file is part of WCSLIB.

    WCSLIB is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
    more details.

    You should have received a copy of the GNU Lesser General Public License
    along with WCSLIB.  If not, see <http://www.gnu.org/licenses/>.

    Correspondence concerning WCSLIB may be directed to:
       Internet email: mcalabre@atnf.csiro.au
       Postal address: Dr. Mark Calabretta
                       Australia Telescope National Facility, CSIRO
                       PO Box 76
                       Epping NSW 1710
                       AUSTRALIA

    Author: Mark Calabretta, Australia Telescope National Facility
    http://www.atnf.csiro.au/~mcalabre/index.html
    $Id: wcstrig.h,v 4.3 2007/12/27 05:41:36 cal103 Exp $
*===========================================================================*/

#ifndef WCSLIB_WCSTRIG
#define WCSLIB_WCSTRIG

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef sunos4
#define HAS_SINCOS
#endif

#ifdef VXWORKS
#ifndef POWERPC
#define HAS_SINCOS
#endif
#endif

#if defined (__GLIBC__) && defined(__USE_GNU)
#define HAS_SINCOS
#endif

#ifdef WCSTRIG_MACRO

/* Macro implementation of the trigd functions. */
#include "wcsmath.h"

#define cosd(X) cos((X)*D2R)
#define sind(X) sin((X)*D2R)
#define tand(X) tan((X)*D2R)
#define acosd(X) acos(X)*R2D
#define asind(X) asin(X)*R2D
#define atand(X) atan(X)*R2D
#define atan2d(Y,X) atan2(Y,X)*R2D
#ifdef HAS_SINCOS
    #define sincosd(X,S,C) sincos((X)*D2R,(S),(C))
#else
    #define sincosd(X,S,C) *(S) = sin((X)*D2R); *(C) = cos((X)*D2R);
#endif

#else

/* Use WCSLIB wrappers or native trigd functions. */

double cosd(double);
double sind(double);
double tand(double);
double acosd(double);
double asind(double);
double atand(double);
double atan2d(double, double);
double sincosd(double, double*, double*);

/* Domain tolerance for asin and acos functions. */
#define WCSTRIG_TOL 1e-10

#endif /* WCSTRIG_MACRO */


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_WCSTRIG */
