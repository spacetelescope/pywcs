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
    $Id: wcstrig.c,v 4.3 2007/12/27 05:41:36 cal103 Exp $
*===========================================================================*/

#include <math.h>
#include "wcsmath.h"
#include "wcstrig.h"

double cosd(angle)

double angle;

{
   /* MGD: This appears to be an attempt at a speed improvement, but
      profiling reveals it is much faster to skip the short-circuiting
      and just use the system cos.
    */

   /* double resid; */

   /* if (fmod(angle,90.0) == 0.0) { */
   /*    resid = fabs(fmod(angle,360.0)); */
   /*    if (resid == 0.0) { */
   /*       return 1.0; */
   /*    } else if (resid == 90.0) { */
   /*       return 0.0; */
   /*    } else if (resid == 180.0) { */
   /*       return -1.0; */
   /*    } else if (resid == 270.0) { */
   /*       return 0.0; */
   /*    } */
   /* } */

   return cos(angle*D2R);
}

/*--------------------------------------------------------------------------*/

double sind(angle)

double angle;

{
   /* MGD: This appears to be an attempt at a speed improvement, but
      profiling reveals it is much faster to skip the short-circuiting
      and just use the system sin.
    */

   /* double resid; */

   /* if (fmod(angle,90.0) == 0.0) { */
   /*    resid = fmod(angle-90.0,360.0); */
   /*    if (resid == 0.0) { */
   /*       return 1.0; */
   /*    } else if (resid == 90.0) { */
   /*       return 0.0; */
   /*    } else if (resid == 180.0) { */
   /*       return -1.0; */
   /*    } else if (resid == 270.0) { */
   /*       return 0.0; */
   /*    } */
   /* } */

   return sin(angle*D2R);
}

/*--------------------------------------------------------------------------*/

double tand(angle)

double angle;

{
   /* MGD: This appears to be an attempt at a speed improvement, but
      profiling reveals it is much faster to skip the short-circuiting
      and just use the system cos.
    */

   /* double resid; */

   /* resid = fmod(angle,360.0); */
   /* if (resid == 0.0 || fabs(resid) == 180.0) { */
   /*    return 0.0; */
   /* } else if (resid == 45.0 || resid == 225.0) { */
   /*    return 1.0; */
   /* } else if (resid == -135.0 || resid == -315.0) { */
   /*    return -1.0; */
   /* } */

   return tan(angle*D2R);
}

/*--------------------------------------------------------------------------*/

double acosd(v)

double v;

{
   if (v >= 1.0) {
      if (v-1.0 <  WCSTRIG_TOL) return 0.0;
   } else if (v == 0.0) {
      return 90.0;
   } else if (v <= -1.0) {
      if (v+1.0 > -WCSTRIG_TOL) return 180.0;
   }

   return acos(v)*R2D;
}

/*--------------------------------------------------------------------------*/

double asind(v)

double v;

{
   if (v <= -1.0) {
      if (v+1.0 > -WCSTRIG_TOL) return -90.0;
   } else if (v == 0.0) {
      return 0.0;
   } else if (v >= 1.0) {
      if (v-1.0 <  WCSTRIG_TOL) return 90.0;
   }

   return asin(v)*R2D;
}

/*--------------------------------------------------------------------------*/

double atand(v)

double v;

{
   if (v == -1.0) {
      return -45.0;
   } else if (v == 0.0) {
      return 0.0;
   } else if (v == 1.0) {
      return 45.0;
   }

   return atan(v)*R2D;
}

/*--------------------------------------------------------------------------*/

double atan2d(y, x)

double x, y;

{
   if (y == 0.0) {
      if (x >= 0.0) {
         return 0.0;
      } else if (x < 0.0) {
         return 180.0;
      }
   } else if (x == 0.0) {
      if (y > 0.0) {
         return 90.0;
      } else if (y < 0.0) {
         return -90.0;
      }
   }

   return atan2(y,x)*R2D;
}
