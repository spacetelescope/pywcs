/*
Copyright (C) 2008 Association of Universities for Research in Astronomy (AURA)

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

/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __PYUTIL_H__
#define __PYUTIL_H__

#include "util.h"

#define PY_ARRAY_UNIQUE_SYMBOL pywcs_numpy_api

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include "str_list_proxy.h"

/* Py_ssize_t for old Pythons */
/* This code is as recommended by: */
/* http://www.python.org/dev/peps/pep-0353/#conversion-guidelines */
#if PY_VERSION_HEX < 0x02050000
   #if !defined(PY_SSIZE_T_MIN)
      typedef int Py_ssize_t;
      # define PY_SSIZE_T_MAX INT_MAX
      # define PY_SSIZE_T_MIN INT_MIN
   #endif
#define lenfunc inquiry
#define ssizeargfunc intargfunc
#define ssizeobjargproc intobjargproc
#endif

PyObject*
PyArrayProxy_New(
    PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data);

static inline void
offset_c_array(
    double* value,
    npy_intp size,
    double offset) {
  double* end = value + size;

  for ( ; value != end; ++value) {
    *value += offset;
  }
}

static inline
void nan2undefined(
    double* value,
    unsigned int nvalues) {

  double* end = value + nvalues;

  for ( ; value != end; ++value) {
    if (isnan64(*value)) {
      *value = UNDEFINED;
    }
  }
}

static inline
void undefined2nan(
    double* value,
    unsigned int nvalues) {

  double* end = value + nvalues;

  for ( ; value != end; ++value) {
    if (*value == UNDEFINED) {
      *value = (double)NPY_NAN;
    }
  }
}

void
preoffset_array(
    PyArrayObject* array,
    int value);

void
unoffset_array(
    PyArrayObject* array,
    int value);

void
copy_array_to_c_double(
    PyArrayObject* array,
    double* dest);

void
copy_array_to_c_int(
    PyArrayObject* array,
    int* dest);

/**
 Returns TRUE if pointer is NULL, and sets Python exception
*/
int
is_null(/*@null@*/ void *);

typedef void (*value_fixer_t)(double*, unsigned int);

void
wcsprm_c2python(
    /*@null@*/ struct wcsprm* x);

void
wcsprm_python2c(
    /*@null@*/ struct wcsprm* x);

/***************************************************************************
 * Exceptions                                                              *
 ***************************************************************************/

extern PyObject* WcsExc_SingularMatrix;
extern PyObject* WcsExc_InconsistentAxisTypes;
extern PyObject* WcsExc_InvalidTransform;
extern PyObject* WcsExc_InvalidCoordinate;
extern PyObject* WcsExc_NoSolution;
extern PyObject* WcsExc_InvalidSubimageSpecification;
extern PyObject* WcsExc_NonseparableSubimageCoordinateSystem;
extern PyObject* WcsExc_NoWcsKeywordsFound;
extern PyObject* WcsExc_InvalidTabularParameters;

/* This is an array mapping the wcs status codes to Python exception
 * types.  The exception string is stored as part of wcslib itself in
 * wcs_errmsg.
 */
extern PyObject** wcs_errexc[14];
#define WCS_ERRMSG_MAX 14
#define WCSFIX_ERRMSG_MAX 11

int
_define_exceptions(PyObject* m);

const char*
wcslib_get_error_message(int stat);

/***************************************************************************
  Property helpers
 ***************************************************************************/
static inline int
check_delete(
    const char* propname,
    PyObject* value) {

  PyObject* ignored;

  if (value == NULL) {
    ignored = PyErr_Format(PyExc_TypeError, "'%s' can not be deleted", propname);
    return -1;
  }

  return 0;
}

static inline PyObject*
get_string(
    /*@unused@*/ const char* propname,
    const char* value) {

  return PyString_FromString(value);
}

int
set_string(
    const char* propname,
    PyObject* value,
    char* dest,
    Py_ssize_t maxlen);

static inline PyObject*
get_bool(
    /*@unused@*/ const char* propname,
    long value) {

  return PyBool_FromLong(value);
}

int
set_bool(
    const char* propname,
    PyObject* value,
    int* dest);

static inline PyObject*
get_int(
    /*@unused@*/ const char* propname,
    long value) {

  return PyInt_FromLong(value);
}

int
set_int(
    const char* propname,
    PyObject* value,
    int* dest);

static inline PyObject*
get_double(
    const char* propname,
    double value) {

  return PyFloat_FromDouble(value);
}

int
set_double(
    const char* propname,
    PyObject* value,
    double* dest);

/*@null@*/ static inline PyObject*
get_double_array(
    /*@unused@*/ const char* propname,
    double* value,
    npy_intp ndims,
    const npy_intp* dims,
    /*@shared@*/ PyObject* owner) {

  return PyArrayProxy_New(owner, ndims, dims, PyArray_DOUBLE, value);
}

int
set_double_array(
    const char* propname,
    PyObject* value,
    npy_int ndims,
    const npy_intp* dims,
    double* dest);

/*@null@*/ static inline PyObject*
get_int_array(
    /*@unused@*/ const char* propname,
    int* value,
    npy_int ndims,
    const npy_intp* dims,
    /*@shared@*/ PyObject* owner) {

  return PyArrayProxy_New(owner, ndims, dims, PyArray_INT, value);
}

int
set_int_array(
    const char* propname,
    PyObject* value,
    npy_int ndims,
    const npy_intp* dims,
    int* dest);

/* Defined in str_list_proxy.h */
PyObject *
PyStrListProxy_New(
    PyObject* owner,
    Py_ssize_t size,
    char (*array)[72]);

static inline PyObject*
get_str_list(
    /*@unused@*/ const char* propname,
    char (*array)[72],
    Py_ssize_t len,
    PyObject* owner) {

  return PyStrListProxy_New(owner, len, array);
}

int
set_str_list(
    const char* propname,
    PyObject* value,
    Py_ssize_t len,
    Py_ssize_t maxlen,
    char (*dest)[72]);

PyObject*
get_pscards(
    const char* propname,
    struct pscard* ps,
    int nps);

int
set_pscards(
    const char* propname,
    PyObject* value,
    struct pscard** ps,
    int *nps,
    int *npsmax);

PyObject*
get_pvcards(
    const char* propname,
    struct pvcard* pv,
    int npv);

int
set_pvcards(
    const char* propname,
    PyObject* value,
    struct pvcard** pv,
    int *npv,
    int *npvmax);

PyObject*
get_deepcopy(
    PyObject* obj,
    PyObject* memo);

/***************************************************************************
  Miscellaneous helper functions
 ***************************************************************************/

int
parse_unsafe_unit_conversion_spec(
    const char* arg, int* ctrl);

#endif /* __PYUTIL_H__ */
