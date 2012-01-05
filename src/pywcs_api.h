/*
Copyright (C) 2008-2012 Association of Universities for Research in Astronomy (AURA)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    3. The name of AURA and its representatives may not be used to
      endorse or promote products derived from this software without
      specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef PYWCS_API_H
#define PYWCS_API_H

#include "pyutil.h"
#include "distortion.h"
#include "pipeline.h"
#include "sip.h"
#include "wcs.h"
#include "wcsprintf.h"

/*
HOW TO UPDATE THE PUBLIC API

This code uses a table of function pointers to dynamically expose the public API
to other code that wants to use pywcs from C.

Each function should be:

  1) Declared, as usual for C, in a .h file

  2) Defined in a .c file that is compiled as part of the _pywcs.so file

  3) Have a macro that maps the function name to a position in the
     function table.  That macro should go in this file (pywcs_api.h)

  4) An entry in the function table, which lives in pywcs_api.c

Every time the function signatures change, or functions are added or
removed from the table, the value of REVISION should be incremented.
This allows for a rudimentary version check upon dynamic linking to
the pywcs module.
 */

#define REVISION 3

#ifdef PYWCS_BUILD

int _setup_api(PyObject* m);

#else

#if defined(NO_IMPORT_PYWCS_API)
extern void** PyWcs_API;
#else
void** PyWcs_API;
#endif /* defined(NO_IMPORT_PYWCS_API) */

/* Function macros that delegate to a function pointer in the PyWCS_API table */
#define PyWcs_GetCVersion (*(int (*)(void)) PyWcs_API[0])
#define wcsprm_python2c (*(void (*)(struct wcsprm*)) PyWcs_API[1])
#define wcsprm_c2python (*(void (*)(struct wcsprm*)) PyWcs_API[2])
#define distortion_lookup_t_init (*(int (*)(distortion_lookup_t* lookup)) PyWcs_API[3])
#define distortion_lookup_t_free (*(void (*)(distortion_lookup_t* lookup)) PyWcs_API[4])
#define get_distortion_offset (*(double (*)(const distortion_lookup_t*, const double* const)) PyWcs_API[5])
#define p4_pix2foc (*(int (*)(const unsigned int, const distortion_lookup_t**, const unsigned int, const double *, double *)) PyWcs_API[6])
#define p4_pix2deltas (*(int (*)(const unsigned int, const distortion_lookup_t**, const unsigned int, const double *, double *)) PyWcs_API[7])
#define sip_clear (*(void (*)(sip_t*) PyWcs_API[8]))
#define sip_init (*(int (*)(sip_t*, unsigned int, double*, unsigned int, double*, unsigned int, double*, unsigned int, double*, double*)) PyWcs_API[9])
#define sip_free (*(void (*)(sip_t*) PyWcs_API[10]))
#define sip_pix2foc (*(int (*)(sip_t*, unsigned int, unsigned int, double*, double*)) PyWcs_API[11])
#define sip_pix2deltas (*(int (*)(sip_t*, unsigned int, unsigned int, double*, double*)) PyWcs_API[12])
#define sip_foc2pix (*(int (*)(sip_t*, unsigned int, unsigned int, double*, double*)) PyWcs_API[13])
#define sip_foc2deltas (*(int (*)(sip_t*, unsigned int, unsigned int, double*, double*)) PyWcs_API[14])
#define pipeline_clear (*(void (*)(pipeline_t*)) PyWcs_API[15])
#define pipeline_init (*(void (*)(pipeline_t*, sip_t*, distortion_lookup_t**, struct wcsprm*)) PyWcs_API[16])
#define pipeline_free (*(void (*)(pipeline_t*)) PyWcs_API[17])
#define pipeline_all_pixel2world (*(int (*)(pipeline_t*, unsigned int, unsigned int, double*, double*)) PyWcs_API[18])
#define pipeline_pix2foc (*(int (*)(pipeline_t*, unsigned int, unsigned int, double*, double*)) PyWcs_API[19])
#define wcsp2s (*(int (*)(struct wcsprm *, int, int, const double[], double[], double[], double[], double[], int[])) PyWcs_API[20])
#define wcss2p (*(int (*)(struct wcsprm *, int, int, const double[], double[], double[], double[], double[], int[])) PyWcs_API[21])
#define wcsprt (*(int (*)(struct wcsprm *)) PyWcs_API[22])
#define wcslib_get_error_message (*(const char* (*)(int)) PyWcs_API[23])
#define wcsprintf_buf (*(const char * (*)()) PyWCS_API[24])

#ifndef NO_IMPORT_PYWCS_API
int
import_pywcs(void) {
  PyObject *pywcs_module = NULL;
  PyObject *c_api        = NULL;
  int       status       = -1;

  #if PY_VERSION_HEX >= 0x03020000
    PyWcs_API = (void **)PyCapsule_Import("pywcs._pywcs._WCS_API", 0);
    if (PyWcs_API == NULL) goto exit;
  #else
    pywcs_module = PyImport_ImportModule("pywcs._pywcs");
    if (pywcs_module == NULL) goto exit;

    c_api = PyObject_GetAttrString(pywcs_module, "_PYWCS_API");
    if (c_api == NULL) goto exit;

    if (PyCObject_Check(c_api)) {
      PyWcs_API = (void **)PyCObject_AsVoidPtr(c_api);
    } else {
      goto exit;
    }
  #endif

  /* Perform runtime check of C API version */
  if (REVISION != PyWcs_GetCVersion()) {
    PyErr_Format(
                 PyExc_ImportError, "module compiled against "        \
                 "ABI version '%x' but this version of pywcs is '%x'", \
                 (int)REVISION, (int)PyWcs_GetCVersion());
    return -1;
  }

 exit:
  Py_XDECREF(pywcs_module);
  Py_XDECREF(c_api);

  return status;
}

#endif /* !defined(NO_IMPORT_PYWCS_API) */

#endif /* PYWCS_BUILD */

#endif /* PYWCS_API_H */
