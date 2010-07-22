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

#define NO_IMPORT_ARRAY

/* util.h must be imported first */
#include "pyutil.h"

/*@null@*/ PyObject*
PyArrayProxy_New(
    /*@shared@*/ PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data) {

  PyArray_Descr* type_descr = NULL;
  PyObject*      result     = NULL;

  type_descr = (PyArray_Descr*)PyArray_DescrFromType(typenum);
  if (type_descr == NULL) {
    return NULL;
  }

  result = (PyObject*)PyArray_NewFromDescr(
      &PyArray_Type,
      type_descr,
      nd, (npy_intp*)dims,
      NULL,
      (void*)data,
      NPY_CONTIGUOUS,
      NULL);

  if (result == NULL) {
    return NULL;
  }
  Py_INCREF(self);
  PyArray_BASE(result) = (PyObject*)self;
  return result;
}

void
preoffset_array(
    PyArrayObject* array,
    int value) {

  npy_intp  size;
  double   *data;

  if (value == 1) {
    return;
  }

  size = PyArray_Size((PyObject*)array);
  data = (double*)PyArray_DATA(array);
  offset_c_array(data, size, (double)(1 - value));
}

void
unoffset_array(
    PyArrayObject* array,
    int value) {

  npy_intp  size;
  double   *data;

  if (value == 1) {
    return;
  }

  size = PyArray_Size((PyObject*)array);
  data = (double*)PyArray_DATA(array);
  offset_c_array(data, size, (double)-(1 - value));
}

void
copy_array_to_c_double(
    PyArrayObject* array,
    double* dest) {

  npy_intp size = 1;
  double*  data = NULL;

  size = PyArray_Size((PyObject*)array);
  data = (double*)PyArray_DATA(array);

  memcpy(dest, data, size * sizeof(double));
}

void
copy_array_to_c_int(
    PyArrayObject* array,
    int* dest) {

  npy_intp size = 1;
  int*     data = NULL;

  size = PyArray_Size((PyObject*)array);
  data = (int*)PyArray_DATA(array);

  memcpy(dest, data, size * sizeof(int));
}

int
is_null(
    /*@null@*/ void *p) {

  if (p == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return 1;
  }
  return 0;
}

/* wcslib represents undefined values using its own special constant,
   UNDEFINED.  To be consistent with the Pythonic way of doing things,
   it's nicer to represent undefined values using NaN.  Unfortunately,
   in order to get nice mutable arrays in Python, Python must be able
   to edit the wcsprm values directly.  The solution is to store NaNs
   in the struct "canonically", but convert those NaNs to/from
   UNDEFINED around every call into a wcslib function.  It's not as
   computationally expensive as it sounds, as all these arrays are
   quite small.
*/

static inline void
wcsprm_fix_values(
    struct wcsprm* x,
    value_fixer_t value_fixer) {

  unsigned int naxis = (unsigned int)x->naxis;

  value_fixer(x->cd, 4);
  value_fixer(x->cdelt, naxis);
  value_fixer(x->crder, naxis);
  value_fixer(x->crota, naxis);
  value_fixer(x->crpix, naxis);
  value_fixer(x->crval, naxis);
  value_fixer(x->csyer, naxis);
  value_fixer(&x->equinox, 1);
  value_fixer(&x->mjdavg, 1);
  value_fixer(&x->mjdobs, 1);
  value_fixer(x->obsgeo, 3);
  value_fixer(&x->restfrq, 1);
  value_fixer(&x->restwav, 1);
  value_fixer(&x->velangl, 1);
  value_fixer(&x->velosys, 1);
  value_fixer(&x->zsource, 1);
}

void
wcsprm_c2python(
    /*@null@*/ struct wcsprm* x) {

  if (x != NULL) {
    wcsprm_fix_values(x, &undefined2nan);
  }
}

void
wcsprm_python2c(
    /*@null@*/ struct wcsprm* x) {

  if (x != NULL) {
    wcsprm_fix_values(x, &nan2undefined);
  }
}

/***************************************************************************
 * Exceptions                                                              *
 ***************************************************************************/

PyObject* WcsExc_SingularMatrix;
PyObject* WcsExc_InconsistentAxisTypes;
PyObject* WcsExc_InvalidTransform;
PyObject* WcsExc_InvalidCoordinate;
PyObject* WcsExc_NoSolution;
PyObject* WcsExc_InvalidSubimageSpecification;
PyObject* WcsExc_NonseparableSubimageCoordinateSystem;
PyObject* WcsExc_NoWcsKeywordsFound;
PyObject* WcsExc_InvalidTabularParameters;

/* This is an array mapping the wcs status codes to Python exception
 * types.  The exception string is stored as part of wcslib itself in
 * wcs_errmsg.
 */
PyObject** wcs_errexc[14];

#define DEFINE_EXCEPTION(exc) \
  WcsExc_##exc = PyErr_NewException("_pywcs." #exc "Error", PyExc_ValueError, NULL); \
  if (WcsExc_##exc == NULL) \
    return 1; \
  PyModule_AddObject(m, #exc "Error", WcsExc_##exc); \

int
_define_exceptions(
    PyObject* m) {

  DEFINE_EXCEPTION(SingularMatrix);
  DEFINE_EXCEPTION(InconsistentAxisTypes);
  DEFINE_EXCEPTION(InvalidTransform);
  DEFINE_EXCEPTION(InvalidCoordinate);
  DEFINE_EXCEPTION(NoSolution);
  DEFINE_EXCEPTION(InvalidSubimageSpecification);
  DEFINE_EXCEPTION(NonseparableSubimageCoordinateSystem);
  DEFINE_EXCEPTION(NoWcsKeywordsFound);
  DEFINE_EXCEPTION(InvalidTabularParameters);
  return 0;
}

const char*
wcslib_get_error_message(int status) {
    return wcsp2s_errmsg[status];
}

/***************************************************************************
  Property helpers
 ***************************************************************************/

#define SHAPE_STR_LEN 128

/* Helper function to display the desired shape of an array as a
   string, eg. 2x2 */
static void
shape_to_string(
    npy_int ndims,
    const npy_intp* dims,
    char* str /* [SHAPE_STR_LEN] */) {

  int i;
  char value[32]; /* More than large enough to hold string rep of a
                     64-bit integer (way overkill) */

  if (ndims > 3) {
    strncpy(str, "ERROR", 6);
    return;
  }

  str[0] = 0;
  for (i = 0; i < ndims; ++i) {
      snprintf(value, 32, "%d", (int)dims[i]);
    strncat(str, value, 32);
    if (i != ndims - 1) {
      strncat(str, "x", 2);
    }
  }
}

/* get_string is inlined */

int
set_string(
    const char* propname,
    PyObject* value,
    char* dest,
    Py_ssize_t maxlen) {

  char*      buffer;
  Py_ssize_t len;

  if (check_delete(propname, value)) {
    return -1;
  }

  if (PyString_AsStringAndSize(value, &buffer, &len) == -1) {
    return -1;
  }

  if (len > maxlen) {
    PyErr_Format(
        PyExc_ValueError,
        "'%s' must be less than %u characters",
        propname,
        (unsigned int)maxlen);
    return -1;
  }

  strncpy(dest, buffer, (size_t)maxlen);

  return 0;
}

/* get_bool is inlined */

int
set_bool(
    const char* propname,
    PyObject* value,
    int* dest) {
  long value_int;

  if (check_delete(propname, value)) {
    return -1;
  }

  value_int = PyInt_AsLong(value);
  if (value_int == -1 && PyErr_Occurred()) {
    return -1;
  }

  *dest = value_int ? 1 : 0;

  return 0;
}

/* get_int is inlined */

int
set_int(
    const char* propname,
    PyObject* value,
    int* dest) {
  long value_int;

  if (check_delete(propname, value)) {
    return -1;
  }

  value_int = PyInt_AsLong(value);
  if (value_int == -1 && PyErr_Occurred()) {
    return -1;
  }

  if ((unsigned long)value_int > 0x7fffffff) {
    return -1;
  }

  *dest = (int)value_int;

  return 0;
}

/* get_double is inlined */

int
set_double(
    const char* propname,
    PyObject* value,
    double* dest) {

  if (check_delete(propname, value)) {
    return -1;
  }

  *dest = PyFloat_AsDouble(value);

  if (PyErr_Occurred()) {
    return -1;
  } else {
    return 0;
  }
}

/* get_double_array is inlined */

int
set_double_array(
    const char* propname,
    PyObject* value,
    npy_int ndims,
    const npy_intp* dims,
    double* dest) {

  PyArrayObject* value_array = NULL;
  npy_int        i           = 0;
  PyObject*      ignored     = NULL;
  char           shape_str[SHAPE_STR_LEN];

  if (check_delete(propname, value)) {
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_DOUBLE,
                                                          ndims, ndims);
  if (value_array == NULL) {
    return -1;
  }

  if (dims != NULL) {
    for (i = 0; i < ndims; ++i) {
      if (PyArray_DIM(value_array, i) != dims[i]) {
        shape_to_string(ndims, dims, shape_str);
        ignored = PyErr_Format(
            PyExc_ValueError,
            "'%s' array is the wrong shape, must be %s",
            propname, shape_str);
        Py_DECREF(value_array);
        return -1;
      }
    }
  }

  copy_array_to_c_double(value_array, dest);

  Py_DECREF(value_array);

  return 0;
}

int
set_int_array(
    const char* propname,
    PyObject* value,
    npy_int ndims,
    const npy_intp* dims,
    int* dest) {
  PyArrayObject* value_array = NULL;
  npy_int        i           = 0;
  PyObject*      ignored     = NULL;
  char           shape_str[SHAPE_STR_LEN];

  if (check_delete(propname, value)) {
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_INT,
                                                          ndims, ndims);
  if (value_array == NULL) {
    return -1;
  }

  if (dims != NULL) {
    for (i = 0; i < ndims; ++i) {
      if (PyArray_DIM(value_array, i) != dims[i]) {
        shape_to_string(ndims, dims, shape_str);
        ignored = PyErr_Format(
            PyExc_ValueError,
            "'%s' array is the wrong shape, must be %s",
            propname, shape_str);
        Py_DECREF(value_array);
        return -1;
      }
    }
  }

  copy_array_to_c_int(value_array, dest);

  Py_DECREF(value_array);

  return 0;
}

/* get_str_list is inlined */

int
set_str_list(
    const char* propname,
    PyObject* value,
    Py_ssize_t len,
    Py_ssize_t maxlen,
    char (*dest)[72]) {

  PyObject*  str      = NULL;
  char*      str_char = NULL;
  Py_ssize_t str_len  = 0;
  Py_ssize_t i        = 0;
  PyObject*  ignored  = NULL;

  if (check_delete(propname, value)) {
    return -1;
  }

  if (maxlen == 0) {
    maxlen = 68;
  }

  if (!PySequence_Check(value)) {
    ignored = PyErr_Format(
        PyExc_TypeError,
        "'%s' must be a sequence of strings",
        propname);
    return -1;
  }

  if (PySequence_Size(value) != len) {
    ignored = PyErr_Format(
        PyExc_ValueError,
        "len(%s) must be %u",
        propname,
        (unsigned int)len);
    return -1;
  }

  /* We go through the list twice, once to verify that the list is
     in the correct format, and then again to do the data copy.  This
     way, we won't partially copy the contents and then throw an
     exception. */
  for (i = 0; i < len; ++i) {
    str = PySequence_GetItem(value, i);
    if (str == NULL) {
      return -1;
    }

    if (!PyString_CheckExact(str)) {
      ignored = PyErr_Format(
          PyExc_TypeError,
          "'%s' must be a sequence of strings",
          propname);
      Py_DECREF(str);
      return -1;
    }
    if (PyString_Size(str) > maxlen) {
      ignored = PyErr_Format(
          PyExc_TypeError,
          "Each string in '%s' must be less than %u characters",
          propname, (unsigned int)maxlen);
      Py_DECREF(str);
      return -1;
    }

    Py_DECREF(str);
  }

  for (i = 0; i < len; ++i) {
    str = PySequence_GetItem(value, i);
    if (str == NULL) {
      /* Theoretically, something has gone really wrong here, since
         we've already verified the list. */
      ignored = PyErr_Format(
          PyExc_RuntimeError,
          "Input values have changed underneath us.  Something is seriously wrong.");
      return -1;
    }

    /* We already know its a string of the correct length */
    if (PyString_AsStringAndSize(str, &str_char, &str_len)) {
      /* Theoretically, something has gone really wrong here, since
         we've already verified the list. */
      ignored = PyErr_Format(
          PyExc_RuntimeError,
          "Input values have changed underneath us.  Something is seriously wrong.");
      Py_DECREF(str);
      return -1;
    }

    strncpy(dest[i], str_char, (size_t)maxlen);

    Py_DECREF(str);
  }

  return 0;
}

/*@null@*/ PyObject*
get_pscards(
    /*@unused@*/ const char* propname,
    struct pscard* ps,
    int nps) {

  PyObject*  result    = NULL;
  PyObject*  subresult = NULL;
  Py_ssize_t i         = 0;

  if (nps < 0) {
    return NULL;
  }

  result = PyList_New((Py_ssize_t)nps);
  if (result == NULL) {
    return NULL;
  }

  for (i = 0; i < (Py_ssize_t)nps; ++i) {
    subresult = Py_BuildValue("iis", ps[i].i, ps[i].m, ps[i].value);
    if (subresult == NULL) {
      Py_DECREF(result);
      return NULL;
    }

    if (PyList_SetItem(result, i, subresult)) {
      Py_DECREF(subresult);
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

int
set_pscards(
    /*@unused@*/ const char* propname,
    PyObject* value,
    struct pscard** ps,
    int *nps,
    int *npsmax) {

  PyObject*  subvalue  = NULL;
  Py_ssize_t i         = 0;
  Py_ssize_t size      = 0;
  int        ival      = 0;
  int        mval      = 0;
  char*      strvalue  = 0;

  if (!PySequence_Check(value))
    return -1;
  size = PySequence_Size(value);
  if (size > 0x7fffffff) {
    return -1;
  }

  /* Verify the entire list for correct types first, so we don't have
     to undo anything copied into the canonical array. */
  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (subvalue == NULL) {
      return -1;
    }
    if (!PyArg_ParseTuple(subvalue, "iis", &ival, &mval, &strvalue)) {
      Py_DECREF(subvalue);
      return -1;
    }
    Py_DECREF(subvalue);
  }

  if (size > (Py_ssize_t)*npsmax) {
    free(*ps);
    *ps = malloc(sizeof(struct pvcard) * size);
    if (*ps == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
      return -1;
    }
    *npsmax = (int)size;
  }

  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (subvalue == NULL) {
      return -1;
    }
    if (!PyArg_ParseTuple(subvalue, "iis", &ival, &mval, &strvalue)) {
      Py_DECREF(subvalue);
      return -1;
    }
    Py_DECREF(subvalue);

    (*ps)[i].i = ival;
    (*ps)[i].m = mval;
    strncpy((*ps)[i].value, strvalue, 72);
    (*ps)[i].value[71] = '\0';
    (*nps) = i + 1;
  }

  return 0;
}

/*@null@*/ PyObject*
get_pvcards(
    /*@unused@*/ const char* propname,
    struct pvcard* pv,
    int npv) {

  PyObject*  result    = NULL;
  PyObject*  subresult = NULL;
  Py_ssize_t i         = 0;

  result = PyList_New((Py_ssize_t)npv);
  if (result == NULL) {
    return NULL;
  }

  for (i = 0; i < (Py_ssize_t)npv; ++i) {
    subresult = Py_BuildValue("iid", pv[i].i, pv[i].m, pv[i].value);
    if (subresult == NULL) {
      Py_DECREF(result);
      return NULL;
    }

    if (PyList_SetItem(result, i, subresult)) {
      Py_DECREF(subresult);
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

int
set_pvcards(
    /*@propname@*/ const char* propname,
    PyObject* value,
    struct pvcard** pv,
    int *npv,
    int *npvmax) {

  PyObject*  subvalue  = NULL;
  int        i         = 0;
  Py_ssize_t size      = 0;
  int        ival      = 0;
  int        mval      = 0;
  double     dblvalue  = 0.0;

  if (!PySequence_Check(value)) {
    return -1;
  }
  size = PySequence_Size(value);
  if (size > 0x7fffffff) {
    return -1;
  }

  /* Verify the entire list for correct types first, so we don't have
     to undo anything copied into the canonical array. */
  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (subvalue == NULL) {
      return -1;
    }
    if (!PyArg_ParseTuple(subvalue, "iid", &ival, &mval, &dblvalue)) {
      Py_DECREF(subvalue);
      return -1;
    }
    Py_DECREF(subvalue);
  }

  if (size > (Py_ssize_t)*npvmax) {
    free(*pv);
    *pv = malloc(sizeof(struct pvcard) * size);
    if (*pv == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
      return -1;
    }
    *npvmax = (int)size;
  }

  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (subvalue == NULL) {
      return -1;
    }
    if (!PyArg_ParseTuple(subvalue, "iid", &ival, &mval, &dblvalue)) {
      Py_DECREF(subvalue);
      return -1;
    }
    Py_DECREF(subvalue);

    (*pv)[i].i = ival;
    (*pv)[i].m = mval;
    (*pv)[i].value = dblvalue;
    (*npv) = i + 1;
  }

  return 0;
}

PyObject*
get_deepcopy(
    PyObject* obj,
    PyObject* memo) {

  if (PyObject_HasAttrString(obj, "__deepcopy__")) {
    return PyObject_CallMethod(obj, "__deepcopy__", "O", memo);
  } else {
    return PyObject_CallMethod(obj, "__copy__", "");
  }
}

