/*============================================================================

  WCSLIB 4.7 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2011, Mark Calabretta

  This file is part of WCSLIB.

  WCSLIB is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

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
  Module author: Michael Droettboom
  http://www.atnf.csiro.au/~mcalabre/index.html
  $Id$
*===========================================================================*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsprintf.h"

/*--------------------------------------------------------------------------*/

int wcserr_ini(
  struct wcserr *err)

{
  if (err == 0x0) {
    return 1;
  }

  memset(err, 0, sizeof(struct wcserr));

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcserr_set(
  struct wcserr **err,
  int status,
  const char *function,
  const char *file,
  int line_no,
  const char *format,
  ...)

{
  va_list argp;

  if (err == 0x0) {
    return status;
  }

  if (*err == 0x0) {
    *err = calloc(1, sizeof(struct wcserr));
  }

  (*err)->status   = status;
  (*err)->function = function;
  (*err)->file     = file;
  (*err)->line_no  = line_no;

  va_start(argp, format);
  vsnprintf((*err)->msg, WCSERR_MSG_LENGTH, format, argp);
  va_end(argp);

  return status;
}

/*--------------------------------------------------------------------------*/

int wcserr_copy(
  const struct wcserr *src,
  struct wcserr **dst)

{
  if (src == 0x0) {
    return 1;
  }

  if (dst == 0x0) {
    return src->status;
  }

  if (*dst == 0x0) {
    *dst = calloc(1, sizeof(struct wcserr));
  }

  memcpy(*dst, src, sizeof(struct wcserr));
  return (*dst)->status;
}

/*--------------------------------------------------------------------------*/

int wcserr_prt(
  const struct wcserr *err,
  const char *prefix)

{
  if (err == 0x0) {
    return 0;
  }

  if (prefix == 0x0) prefix = "";

  wcsprintf("%sERROR %d in %s() at line %d of file %s:\n%s  %s.\n", prefix,
    err->status, err->function, err->line_no, err->file, prefix, err->msg);

  return 0;
}
