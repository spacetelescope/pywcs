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

#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <wcs.h>
#include <wcsfix.h>
#include <wcshdr.h>
#include <wcsmath.h>
#include "isnan.h"

/* Py_ssize_t for old Pythons */
/* This code is as recommended by: */
/* http://www.python.org/dev/peps/pep-0353/#conversion-guidelines */
#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
# define PY_SSIZE_T_MAX INT_MAX
# define PY_SSIZE_T_MIN INT_MIN
#endif

/*
 It gets to be really tedious to type long docstrings in ANSI C syntax (since
 multi-line strings literals are not valid).  Therefore, the docstrings are
 written in doc/docstrings.py, which are then converted by setup.py into
 docstrings.h, which we include here.
*/
#include "docstrings.h"

/* TODO: What does this need to read the tables for, and how does that
 * work?
 */

/***************************************************************************
 * Helper functions                                                        *
 ***************************************************************************/

enum e_altlin {
  has_pc = 1,
  has_cd = 2,
  has_crota = 4
};

static int
parse_unsafe_unit_conversion_spec(const char* arg) {
  int ctrl = 0;
  if (strchr(arg, 's') || strchr(arg, 'S'))
    ctrl |= 1;
  if (strchr(arg, 'h') || strchr(arg, 'H'))
    ctrl |= 2;
  if (strchr(arg, 'd') || strchr(arg, 'D'))
    ctrl |= 4;
  return ctrl;
}

static inline void
offset_c_array(double* value, size_t size, double offset) {
  double* end = value + size;

  for ( ; value != end; ++value)
    *value += offset;
}

static void
offset_array(PyArrayObject* array, double value) {
  int     size = 1;
  int     i    = 0;
  double* data = NULL;

  for (i = 0; i < PyArray_NDIM(array); ++i)
    size *= PyArray_DIM(array, i);

  data = (double*)PyArray_DATA(array);

  offset_c_array(data, size, value);
}

static void
copy_array_to_c_double(PyArrayObject* array, double* dest) {
  int     size = 1;
  int     i    = 0;
  double* data = NULL;

  data = (double*)PyArray_DATA(array);

  for (i = 0; i < PyArray_NDIM(array); ++i)
    size *= PyArray_DIM(array, i);

  for (i = 0; i < size; ++i, ++dest, ++data) {
    if (isnan64(*data))
      *dest = UNDEFINED;
    else
      *dest = *data;
  }
}

static void
copy_array_to_c_int(PyArrayObject* array, int* dest) {
  int  size = 1;
  int  i    = 0;
  int* data = NULL;

  data = (int*)PyArray_DATA(array);

  for (i = 0; i < PyArray_NDIM(array); ++i)
    size *= PyArray_DIM(array, i);

  for (i = 0; i < size; ++i, ++dest, ++data)
    *dest = *data;
}

static int
is_valid_alt_key(const char* key) {
  if (key[1] != 0 ||
      !(key[0] == ' ' ||
        (key[0] >= 'A' && key[0] <= 'Z'))) {
    PyErr_SetString(PyExc_ValueError, "key must be ' ' or 'A'-'Z'");
    return 0;
  }

  return 1;
}

static inline
void nan2undefined(double* value, size_t nvalues) {
  double* end = value + nvalues;

  for ( ; value != end; ++value)
    if (isnan64(*value))
      *value = UNDEFINED;
}

static inline
void undefined2nan(double* value, size_t nvalues) {
  double* end = value + nvalues;
  double  v   = 0;

  for ( ; value != end; ++value) {
    v = *value;
    *value = (v == UNDEFINED) ? NAN : v;
  }
}

/***************************************************************************
 * Exceptions                                                              *
 ***************************************************************************/

static PyObject* WcsExc_SingularMatrix;
static PyObject* WcsExc_InconsistentAxisTypes;
static PyObject* WcsExc_InvalidTransform;
static PyObject* WcsExc_InvalidCoordinate;
static PyObject* WcsExc_NoSolution;
static PyObject* WcsExc_InvalidSubimageSpecification;
static PyObject* WcsExc_NonseparableSubimageCoordinateSystem;

/* This is an array mapping the wcs status codes to Python exception
 * types.  The exception string is stored as part of wcslib itself in
 * wcs_errmsg.
 */
static PyObject** wcs_errexc[] = {
  /* 0 */ NULL,                         /* Success */
  /* 1 */ &PyExc_MemoryError,           /* Null wcsprm pointer passed */
  /* 2 */ &PyExc_MemoryError,           /* Memory allocation failed */
  /* 3 */ &WcsExc_SingularMatrix,       /* Linear transformation matrix is singular */
  /* 4 */ &WcsExc_InconsistentAxisTypes, /* Inconsistent or unrecognized coordinate axis types */
  /* 5 */ &PyExc_ValueError,            /* Invalid parameter value */
  /* 6 */ &WcsExc_InvalidTransform,     /* Invalid coordinate transformation parameters */
  /* 7 */ &WcsExc_InvalidTransform,     /* Ill-conditioned coordinate transformation parameters */
  /* 8 */ &WcsExc_InvalidCoordinate,    /* One or more of the pixel coordinates were invalid, */
                                        /* as indicated by the stat vector */
  /* 9 */ &WcsExc_InvalidCoordinate,    /* One or more of the world coordinates were invalid, */
                                        /* as indicated by the stat vector */
  /*10 */ &WcsExc_InvalidCoordinate,    /* Invalid world coordinate */
  /*11 */ &WcsExc_NoSolution,           /* no solution found in the specified interval */
  /*12 */ &WcsExc_InvalidSubimageSpecification, /* Invalid subimage specification (no spectral axis) */
  /*13 */ &WcsExc_NonseparableSubimageCoordinateSystem /* Non-separable subimage coordinate system */
};
#define WCS_ERRMSG_MAX 14
#define WCSFIX_ERRMSG_MAX 11

/***************************************************************************
 * PyWcsprm object
 ***************************************************************************/

static PyTypeObject PyWcsprmType;

typedef struct {
  PyObject_HEAD
  struct wcsprm* x;
} PyWcsprm;

static PyObject*
PyWcsprmArrayProxy_New(PyWcsprm* self, int nd, const npy_intp* dims,
                       int typenum, const void* data) {
  PyArray_Descr* type_descr = NULL;
  PyObject*      result     = NULL;

  type_descr = PyArray_DescrFromType(typenum);
  if (type_descr == NULL)
    return NULL;
  result = PyArray_NewFromDescr(&PyArray_Type, type_descr,
                                nd, (npy_intp*)dims, NULL, (void*)data,
                                NPY_CONTIGUOUS | NPY_WRITEABLE,
                                NULL);

  if (result == NULL)
    return NULL;
  Py_INCREF(self);
  PyArray_BASE(result) = (PyObject*)self;
  return result;
}

static PyObject *
PyWcsprmListProxy_New(PyWcsprm* owner, Py_ssize_t size, char (*array)[72]);

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
static void
PyWcsprm_Python2C(PyWcsprm* self) {
  struct wcsprm* x     = self->x;
  int            naxis = x->naxis;

  nan2undefined(x->cd, 4);
  nan2undefined(x->cdelt, naxis);
  nan2undefined(x->crder, naxis);
  nan2undefined(x->crota, naxis);
  nan2undefined(x->crpix, naxis);
  nan2undefined(x->crval, naxis);
  nan2undefined(x->csyer, naxis);
  nan2undefined(&x->equinox, 1);
  nan2undefined(&x->mjdavg, 1);
  nan2undefined(&x->mjdobs, 1);
  nan2undefined(x->obsgeo, 3);
  nan2undefined(&x->restfrq, 1);
  nan2undefined(&x->restwav, 1);
  nan2undefined(&x->velangl, 1);
  nan2undefined(&x->velosys, 1);
  nan2undefined(&x->zsource, 1);
  offset_c_array(x->crpix, naxis, 1.0);
}

static void
PyWcsprm_C2Python(PyWcsprm* self) {
  struct wcsprm* x     = self->x;
  int            naxis = x->naxis;

  undefined2nan(x->cd, 4);
  undefined2nan(x->cdelt, naxis);
  undefined2nan(x->crder, naxis);
  undefined2nan(x->crota, naxis);
  undefined2nan(x->crpix, naxis);
  undefined2nan(x->crval, naxis);
  undefined2nan(x->csyer, naxis);
  undefined2nan(&x->equinox, 1);
  undefined2nan(&x->mjdavg, 1);
  undefined2nan(&x->mjdobs, 1);
  undefined2nan(x->obsgeo, 3);
  undefined2nan(&x->restfrq, 1);
  undefined2nan(&x->restwav, 1);
  undefined2nan(&x->velangl, 1);
  undefined2nan(&x->velosys, 1);
  undefined2nan(&x->zsource, 1);
  offset_c_array(x->crpix, naxis, -1.0);
}

/***************************************************************************
 * PyWcsprm methods
 */

static void
PyWcsprm_dealloc(PyWcsprm* self) {
  wcsfree(self->x);
  free(self->x);
  self->ob_type->tp_free((PyObject*)self);
}

static PyWcsprm*
PyWcsprm_cnew(struct wcsprm* prm) {
  PyWcsprm* self;
  self = (PyWcsprm*)(&PyWcsprmType)->tp_alloc(&PyWcsprmType, 0);
  if (self != NULL)
    self->x = prm;
  return self;
}

static PyObject *
PyWcsprm_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
  PyWcsprm* self;
  self = (PyWcsprm*)type->tp_alloc(type, 0);
  if (self != NULL)
    self->x = NULL;
  return (PyObject*)self;
}

static int
PyWcsprm_init(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  struct wcsprm* obj;
  int            status;
  PyObject*      header_obj    = NULL;
  char *         header        = NULL;
  Py_ssize_t     header_length = 0;
  Py_ssize_t     nkeyrec       = 0;
  char *         key           = " ";
  int            relax         = 0;
  int            naxis         = 2;
  int            ctrl          = 0;
  int            nreject       = 0;
  int            nwcs          = 0;
  struct wcsprm* wcs           = NULL;
  int            i             = 0;

  /* If we've already been initialized, free our wcsprm object */
  if (self->x != NULL) {
    wcsfree(self->x);
    free(self->x);
    self->x = NULL;
  }

  const char* keywords[] = {"header", "key", "relax", "naxis", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Osii", (char **)keywords,
                                   &header_obj, &key, &relax, &naxis))
    return -1;

  if (header_obj == NULL || header_obj == Py_None) {
    if (relax != 0 || key[0] != ' ' || key[1] != 0) {
      PyErr_SetString(PyExc_ValueError,
                      "If no header is provided, relax and key may not be "
                      "provided either.");
      return -1;
    }

    if (naxis < 1 || naxis > 15) {
      PyErr_SetString(PyExc_ValueError,
                      "naxis must be in range 1-15");
      return -1;
    }

    obj = malloc(sizeof(struct wcsprm));
    if (obj == NULL) {
      PyErr_SetString(PyExc_MemoryError,
                      "Could not allocate wcsprm object");
      return -1;
    }

    obj->flag = -1;
    status = wcsini(1, naxis, obj);

    if (status != 0) {
      free(obj);
      PyErr_SetString(PyExc_MemoryError,
                      "Could not initialize wcsprm object");
      return -1;
    }

    self->x = obj;
    PyWcsprm_C2Python(self);

    return 0;
  } else { /* header != NULL */
    if (PyString_AsStringAndSize(header_obj, &header, &header_length))
      return -1;

    if (relax)
      relax = WCSHDR_all;

    if (!is_valid_alt_key(key))
      return -1;

    if (naxis != 2) {
      PyErr_SetString(PyExc_ValueError,
                      "naxis may not be provided if a header is provided.");
      return -1;
    }

    nkeyrec = header_length / 80;

    status = wcspih(header,
                    nkeyrec,
                    relax,
                    ctrl,
                    &nreject,
                    &nwcs,
                    &wcs);

    if (status != 0) {
      PyErr_SetString(PyExc_MemoryError,
                      "Memory allocation error.");
      return -1;
    }

    /* Find the desired WCS */
    for (i = 0; i < nwcs; ++i)
      if (wcs[i].alt[0] == key[0])
        break;

    if (i >= nwcs) {
      wcsvfree(&nwcs, &wcs);
      PyErr_Format(PyExc_KeyError,
                   "No WCS with key '%s' was found in the given header",
                   key);
      return -1;
    }

    obj = malloc(sizeof(struct wcsprm));
    if (obj == NULL) {
      wcsvfree(&nwcs, &wcs);
      PyErr_SetString(PyExc_MemoryError,
                      "Could not allocate wcsprm object");
      return -1;
    }

    obj->flag = -1;
    if (wcscopy(1, wcs + i, obj) != 0) {
      wcsfree(obj);
      free(obj);
      wcsvfree(&nwcs, &wcs);
      PyErr_SetString(PyExc_MemoryError,
                      "Could not initialize wcsprm object");
      return -1;
    }

    self->x = obj;
    PyWcsprm_C2Python(self);
    wcsvfree(&nwcs, &wcs);
    return 0;
  }
}

static PyObject*
PyWcsprm_repr(PyWcsprm* self) {
  if (self->x == NULL || self->x->wcsname == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  return PyString_FromFormat("<WCS object '%s'>", self->x->wcsname);
}

static PyObject*
PyWcsprm_copy(PyWcsprm* self) {
  PyWcsprm*      copy      = NULL;
  struct wcsprm* copy_data = NULL;
  int            status;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  copy_data = malloc(sizeof(struct wcsprm));
  if (copy_data == NULL) {
      PyErr_SetString(PyExc_MemoryError,
                      "Memory allocation failed.");
      return NULL;
  }

  copy = PyWcsprm_cnew(copy_data);
  if (copy == NULL) {
    return NULL;
  }

  status = wcscopy(1, self->x, copy->x);

  if (status == 0) {
    return (PyObject*)copy;
  } else if (status > 0 && status < WCS_ERRMSG_MAX) {
    Py_XDECREF(copy);
    PyErr_SetString(*wcs_errexc[status], wcscopy_errmsg[status]);
    return NULL;
  } else {
    Py_XDECREF(copy);
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_celfix(PyWcsprm* self) {
  int status = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  PyWcsprm_Python2C(self);
  status = celfix(self->x);
  PyWcsprm_C2Python(self);

  if (status == -1 || status == 0) {
    return PyInt_FromLong(status);
  } else if (status > 0 && status < WCSFIX_ERRMSG_MAX) {
    PyErr_SetString(PyExc_ValueError, wcsfix_errmsg[status]);
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_cylfix(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  PyObject*      naxis_obj   = NULL;
  PyArrayObject* naxis_array = NULL;
  int*           naxis       = NULL;
  int            status      = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  const char* keywords[] = {"naxis", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", (char **)keywords,
                                   &naxis_obj))
    return NULL;

  if (naxis_obj != NULL) {
    naxis_array = (PyArrayObject*)PyArray_ContiguousFromAny(naxis_obj, 1, 1,
                                                               PyArray_INT);
    if (naxis_array == NULL)
      return NULL;
    if (PyArray_DIM(naxis_array, 0) != self->x->naxis) {
      PyErr_Format(PyExc_ValueError,
                   "naxis must be same length as the number of axes of "
                   "the WCS object (%d).",
                   self->x->naxis);
      Py_DECREF(naxis_array);
      return NULL;
    }
    naxis = (int*)PyArray_DATA(naxis_array);
  }

  PyWcsprm_Python2C(self);
  status = cylfix(naxis, self->x);
  PyWcsprm_C2Python(self);

  Py_XDECREF(naxis_array);

  if (status == -1 || status == 0) {
    return PyInt_FromLong(status);
  } else if (status > 0 && status < WCSFIX_ERRMSG_MAX) {
    PyErr_SetString(PyExc_ValueError, wcsfix_errmsg[status]);
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_datfix(PyWcsprm* self) {
  int status = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  status = datfix(self->x);

  if (status == -1 || status == 0) {
    return PyInt_FromLong(status);
  } else if (status > 0 && status < WCSFIX_ERRMSG_MAX) {
    PyErr_SetString(PyExc_ValueError, wcsfix_errmsg[status]);
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_fix(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  char*          translate_units = NULL;
  int            ctrl            = 0;
  PyObject*      naxis_obj       = NULL;
  PyArrayObject* naxis_array     = NULL;
  int*           naxis           = NULL;
  int            stat[NWCSFIX];
  int            status          = 0;
  PyObject*      subresult;
  PyObject*      result;
  int            i               = 0;
  int            msg_index       = 0;
  struct message_map_entry {
    const char* name;
    const int index;
  };
  const struct message_map_entry message_map[NWCSFIX] = {
    {"datfix", DATFIX},
    {"unitfix", UNITFIX},
    {"celfix", CELFIX},
    {"spcfix", SPCFIX},
    {"cylfix", CYLFIX}
  };

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  const char* keywords[] = {"translate_units", "naxis", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|sO", (char **)keywords,
                                   &translate_units, &naxis_obj))
    return NULL;

  if (translate_units != NULL)
    ctrl = parse_unsafe_unit_conversion_spec(translate_units);

  if (naxis_obj != NULL) {
    naxis_array = (PyArrayObject*)PyArray_ContiguousFromAny(naxis_obj, 1, 1,
                                                               PyArray_INT);
    if (naxis_array == NULL)
      return NULL;
    if (PyArray_DIM(naxis_array, 0) != self->x->naxis) {
      PyErr_Format(PyExc_ValueError,
                   "naxis must be same length as the number of axes of "
                   "the WCS object (%d).",
                   self->x->naxis);
      Py_DECREF(naxis_array);
      return NULL;
    }
    naxis = (int*)PyArray_DATA(naxis_array);
  }

  PyWcsprm_Python2C(self);
  status = wcsfix(ctrl, naxis, self->x, stat);
  PyWcsprm_C2Python(self);

  /* We're done with this already, so deref now so we don't have to remember
     later */
  Py_XDECREF(naxis_array);

  result = PyDict_New();
  if (result == NULL)
    return NULL;

  for (i = 0; i < NWCSFIX; ++i) {
    msg_index = stat[message_map[i].index];
    if (msg_index >= 0 && msg_index < 11) {
      subresult = PyString_FromString(wcsfix_errmsg[msg_index]);
      if (subresult == NULL ||
          PyDict_SetItemString(result, "datfix", subresult)) {
        Py_XDECREF(subresult);
        Py_XDECREF(result);
        return NULL;
      }
      Py_XDECREF(subresult);
    }
  }

  return result;
}

static PyObject*
PyWcsprm_get_ps(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  PyObject* result    = NULL;
  PyObject* subresult = NULL;
  int       i         = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if (self->x->ps == NULL) {
    PyErr_SetString(PyExc_AssertionError, "No PSi_ma records present.");
    return NULL;
  }

  result = PyList_New(self->x->nps);
  if (result == NULL)
    return NULL;

  for (i = 0; i < self->x->nps; ++i) {
    subresult = Py_BuildValue("iis",
                              self->x->ps[i].i,
                              self->x->ps[i].m,
                              self->x->ps[i].value);
    if (subresult == NULL || PyList_SetItem(result, i, subresult)) {
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

static PyObject*
PyWcsprm_get_pv(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  PyObject* result    = NULL;
  PyObject* subresult = NULL;
  int       i         = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if (self->x->pv == NULL) {
    PyErr_SetString(PyExc_AssertionError, "No PVi_ma records present.");
    return NULL;
  }

  result = PyList_New(self->x->npv);
  if (result == NULL)
    return NULL;

  for (i = 0; i < self->x->npv; ++i) {
    subresult = Py_BuildValue("iid",
                              self->x->pv[i].i,
                              self->x->pv[i].m,
                              self->x->pv[i].value);
    if (subresult == NULL || PyList_SetItem(result, i, subresult)) {
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

// TODO: Replace these with some "hasattr" magic.

static PyObject*
PyWcsprm_has_cdi_ja(PyWcsprm* self) {
  int result = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  result = self->x->altlin & has_cd;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_has_crotaia(PyWcsprm* self) {
  int result = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  result = self->x->altlin & has_crota;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_has_pci_ja(PyWcsprm* self) {
  int result = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  result = self->x->altlin & has_pc;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_is_unity(PyWcsprm* self) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyBool_FromLong(self->x->lin.unity);
}

static PyObject*
PyWcsprm_mix(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  int            mixpix     = 0;
  int            mixcel     = 0;
  double         vspan[2]   = {0, 0};
  double         vstep      = 0;
  int            viter      = 0;
  Py_ssize_t     naxis      = 0;
  PyObject*      world_obj  = NULL;
  PyObject*      pixcrd_obj = NULL;
  PyArrayObject* world      = NULL;
  PyArrayObject* phi        = NULL;
  PyArrayObject* theta      = NULL;
  PyArrayObject* imgcrd     = NULL;
  PyArrayObject* pixcrd     = NULL;
  int            status     = -1;
  PyObject*      result     = NULL;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if (!PyArg_ParseTuple(args, "ii(dd)diOO",
                        &mixpix, &mixcel, &vspan[0], &vspan[1],
                        &vstep, &viter, &world_obj, &pixcrd_obj))
    return NULL;

  if (viter < 5 || viter > 10) {
    PyErr_SetString(PyExc_ValueError,
                    "viter must be in the range 5 - 10");
    goto __PyWcsprm_mix_exit;
  }

  world = (PyArrayObject*)PyArray_ContiguousFromAny
    (world_obj, PyArray_DOUBLE, 1, 1);
  if (world == NULL) {
    PyErr_SetString(PyExc_TypeError,
                    "Argument 6 (world) must be a 1-dimensional numpy array");
    goto __PyWcsprm_mix_exit;
  }
  if (PyArray_DIM(world, 0) != self->x->naxis) {
    PyErr_Format(PyExc_TypeError,
                 "Argument 6 (world) must be the same length as the number "
                 "of axes (%d)",
                 self->x->naxis);
    goto __PyWcsprm_mix_exit;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny
    (pixcrd_obj, PyArray_DOUBLE, 1, 1);
  if (pixcrd == NULL) {
    PyErr_SetString(PyExc_TypeError,
                    "Argument 7 (pixcrd) must be a 1-dimensional numpy array");
    goto __PyWcsprm_mix_exit;
  }
  if (PyArray_DIM(pixcrd, 0) != self->x->naxis) {
    PyErr_Format(PyExc_TypeError,
                 "Argument 7 (pixcrd) must be the same length as the "
                 "number of axes (%d)",
                 self->x->naxis);
    goto __PyWcsprm_mix_exit;
  }

  if (mixpix < 1 || mixpix > 2) {
    PyErr_SetString(PyExc_ValueError,
                    "Argument 1 (mixpix) must specify a pixel coordinate "
                    "axis number");
    goto __PyWcsprm_mix_exit;
  }

  if (mixcel < 1 || mixcel > 2) {
    PyErr_SetString(PyExc_ValueError,
                    "Argument 2 (mixcel) must specify a celestial coordinate "
                    "axis number (1 for latitude, 2 for longitude)");
    goto __PyWcsprm_mix_exit;
  }

  /* Now we allocate a bunch of numpy arrays to store the
   * results in.
   */
  naxis = (Py_ssize_t)&self->x->naxis;
  phi = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto __PyWcsprm_mix_exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, PyArray_DOUBLE);
  if (theta == NULL) {
    status = 2;
    goto __PyWcsprm_mix_exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, PyArray_DOUBLE);
  if (imgcrd == NULL) {
    status = 2;
    goto __PyWcsprm_mix_exit;
  }

  /* Convert pixel coordinates to 1-based */
  offset_array(pixcrd, 1.0);
  PyWcsprm_Python2C(self);
  status = wcsmix(self->x,
                  mixpix,
                  mixcel,
                  vspan,
                  vstep,
                  viter,
                  (double*)PyArray_DATA(world),
                  (double*)PyArray_DATA(phi),
                  (double*)PyArray_DATA(theta),
                  (double*)PyArray_DATA(imgcrd),
                  (double*)PyArray_DATA(pixcrd));
  PyWcsprm_C2Python(self);
  /* Convert pixel coordinates back to 0-based) */
  offset_array(pixcrd, -1.0);

  if (status == 0) {
    result = PyDict_New();
    if (result == NULL ||
        PyDict_SetItemString(result, "imgcrd", (PyObject*)imgcrd) ||
        PyDict_SetItemString(result, "phi", (PyObject*)phi) ||
        PyDict_SetItemString(result, "theta", (PyObject*)theta))
      status = 2;
  }

 __PyWcsprm_mix_exit:
  Py_XDECREF(world);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(imgcrd);
  Py_XDECREF(pixcrd);

  if (status == -1) {
    /* The error message has already been set */
    return NULL;
  } else if (status == 0) {
    return result;
  } else if (status > 0 && status < WCS_ERRMSG_MAX) {
    PyErr_SetString(*wcs_errexc[status], wcsmix_errmsg[status]);
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_p2s(PyWcsprm* self, PyObject* arg) {
  PyArrayObject* pixcrd  = NULL;
  PyArrayObject* imgcrd  = NULL;
  PyArrayObject* phi     = NULL;
  PyArrayObject* theta   = NULL;
  PyArrayObject* world   = NULL;
  PyArrayObject* stat    = NULL;
  PyObject*      result  = NULL;
  int            status  = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny
    (arg, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL)
    return NULL;

  /* Now we allocate a bunch of numpy arrays to store the results in.
   */
  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (imgcrd == NULL) {
    status = 2;
    goto __PyWcsprm_p2s_exit;
  }

  phi = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto __PyWcsprm_p2s_exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (theta == NULL) {
    status = 2;
    goto __PyWcsprm_p2s_exit;
  }

  world = (PyArrayObject*)PyArray_SimpleNew
    (2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (world == NULL) {
    status = 2;
    goto __PyWcsprm_p2s_exit;
  }

  stat = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(pixcrd), PyArray_INT);
  if (stat == NULL) {
    status = 2;
    goto __PyWcsprm_p2s_exit;
  }

  /* Adjust pixel coordinates to be 1-based */
  offset_array(pixcrd, 1.0);

  /* Make the call */
  PyWcsprm_Python2C(self);
  status = wcsp2s(self->x,
                  PyArray_DIM(pixcrd, 0),
                  PyArray_DIM(pixcrd, 1),
                  (double*)PyArray_DATA(pixcrd),
                  (double*)PyArray_DATA(imgcrd),
                  (double*)PyArray_DATA(phi),
                  (double*)PyArray_DATA(theta),
                  (double*)PyArray_DATA(world),
                  (int*)PyArray_DATA(stat));
  PyWcsprm_C2Python(self);
  /* Adjust pixel coordinates back to 0-based */
  offset_array(pixcrd, -1.0);

  if (status == 0 || status == 8) {
    result = PyDict_New();
    if (result == NULL ||
        PyDict_SetItemString(result, "imgcrd", (PyObject*)imgcrd) ||
        PyDict_SetItemString(result, "phi", (PyObject*)phi) ||
        PyDict_SetItemString(result, "theta", (PyObject*)theta) ||
        PyDict_SetItemString(result, "world", (PyObject*)world) ||
        PyDict_SetItemString(result, "stat", (PyObject*)stat))
      status = 2;
  }

 __PyWcsprm_p2s_exit:
  Py_XDECREF(pixcrd);
  Py_XDECREF(imgcrd);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(world);
  Py_XDECREF(stat);

  if (status == 0 || status == 8) {
    return result;
  } else if (status > 0 && status < WCS_ERRMSG_MAX) {
    PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
    return NULL;
  } else if (status == -1) {
    PyErr_SetString(PyExc_ValueError,
                    "Invalid argument.");
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_s2p(PyWcsprm* self, PyObject* arg) {
  PyArrayObject* world   = NULL;
  PyArrayObject* phi     = NULL;
  PyArrayObject* theta   = NULL;
  PyArrayObject* imgcrd  = NULL;
  PyArrayObject* pixcrd  = NULL;
  PyArrayObject* stat    = NULL;
  PyObject*      result  = NULL;
  int            status  = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  world = (PyArrayObject*)PyArray_ContiguousFromAny
    (arg, PyArray_DOUBLE, 2, 2);
  if (world == NULL)
    return NULL;

  /* Now we allocate a bunch of numpy arrays to store the
   * results in.
   */
  phi = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(world), PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto __PyWcsprm_s2p_exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(world), PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto __PyWcsprm_s2p_exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (2, PyArray_DIMS(world), PyArray_DOUBLE);
  if (theta == NULL) {
    status = 2;
    goto __PyWcsprm_s2p_exit;
  }

  pixcrd = (PyArrayObject*)PyArray_SimpleNew
    (2, PyArray_DIMS(world), PyArray_DOUBLE);
  if (pixcrd == NULL) {
    status = 2;
    goto __PyWcsprm_s2p_exit;
  }

  stat = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(world), PyArray_INT);
  if (stat == NULL) {
    status = 2;
    goto __PyWcsprm_s2p_exit;
  }

  /* Make the call */
  PyWcsprm_Python2C(self);
  status = wcss2p(self->x,
                  PyArray_DIM(world, 0),
                  PyArray_DIM(world, 1),
                  (double*)PyArray_DATA(world),
                  (double*)PyArray_DATA(phi),
                  (double*)PyArray_DATA(theta),
                  (double*)PyArray_DATA(imgcrd),
                  (double*)PyArray_DATA(pixcrd),
                  (int*)PyArray_DATA(stat));
  PyWcsprm_C2Python(self);

  /* Adjust pixel coordinates to be zero-based */
  offset_array(pixcrd, -1.0);

  if (status == 0 || status == 9) {
    result = PyDict_New();
    if (result == NULL ||
        PyDict_SetItemString(result, "phi", (PyObject*)phi) ||
        PyDict_SetItemString(result, "theta", (PyObject*)theta) ||
        PyDict_SetItemString(result, "imgcrd", (PyObject*)imgcrd) ||
        PyDict_SetItemString(result, "pixcrd", (PyObject*)pixcrd) ||
        PyDict_SetItemString(result, "stat", (PyObject*)stat))
      status = 2;
  }

 __PyWcsprm_s2p_exit:
  Py_XDECREF(pixcrd);
  Py_XDECREF(imgcrd);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(world);
  Py_XDECREF(stat);

  if (status == 0 || status == 9) {
    return result;
  } else if (status > 0 && status < WCS_ERRMSG_MAX) {
    PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static int
PyWcsprm_cset(PyWcsprm* self) {
  int status = 0;
  PyWcsprm_Python2C(self);
  status = wcsset(self->x);
  PyWcsprm_C2Python(self);

  if (status == 0) {
    return 0;
  } else if (status > 0 && status < WCS_ERRMSG_MAX) {
    PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
    return -1;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return -1;
  }
}

static PyObject*
PyWcsprm_set(PyWcsprm* self) {
  if (PyWcsprm_cset(self))
    return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
PyWcsprm_set_ps(PyWcsprm* self, PyObject* arg, PyObject* kwds) {
  PyObject*  subvalue  = NULL;
  int        i         = 0;
  Py_ssize_t size      = 0;
  int        ival      = 0;
  int        mval      = 0;
  char*      value     = 0;

  if (self->x == NULL || self->x->ps == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if (!PySequence_Check(arg))
    return NULL;
  size = PySequence_Size(arg);

  /* Verify the entire list for correct types first, so we don't have
     to undo anything copied into the canonical array. */
  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(arg, i);
    if (subvalue == NULL)
      return NULL;
    if (!PyArg_ParseTuple(subvalue, "iis", &ival, &mval, &value))
      return NULL;
    Py_DECREF(subvalue);
  }

  if (size > self->x->npsmax) {
    free(self->x->ps);
    self->x->ps = malloc(sizeof(struct pscard) * size);
    if (self->x->ps == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
      return NULL;
    }
    self->x->npsmax = size;
  }

  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(arg, i);
    if (subvalue == NULL)
      return NULL;
    if (!PyArg_ParseTuple(subvalue, "iis", &ival, &mval, &value))
      return NULL;
    Py_DECREF(subvalue);

    self->x->ps[i].i = ival;
    self->x->ps[i].m = mval;
    strncpy(self->x->ps[i].value, value, 72);
    self->x->ps[i].value[71] = 0;
    self->x->nps = i + 1;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
PyWcsprm_set_pv(PyWcsprm* self, PyObject* arg, PyObject* kwds) {
  PyObject*  subvalue  = NULL;
  int        i         = 0;
  Py_ssize_t size      = 0;
  int        ival      = 0;
  int        mval      = 0;
  double     value     = 0;

  if (self->x == NULL || self->x->ps == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if (!PySequence_Check(arg))
    return NULL;
  size = PySequence_Size(arg);

  /* Verify the entire list for correct types first, so we don't have
     to undo anything copied into the canonical array. */
  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(arg, i);
    if (subvalue == NULL)
      return NULL;
    if (!PyArg_ParseTuple(subvalue, "iid", &ival, &mval, &value))
      return NULL;
    Py_DECREF(subvalue);
  }

  if (size > self->x->npvmax) {
    free(self->x->pv);
    self->x->pv = malloc(sizeof(struct pscard) * size);
    if (self->x->pv == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
      return NULL;
    }
    self->x->npvmax = size;
  }

  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(arg, i);
    if (subvalue == NULL)
      return NULL;
    if (!PyArg_ParseTuple(subvalue, "iid", &ival, &mval, &value)) {
      Py_DECREF(subvalue);
      return NULL;
    }
    Py_DECREF(subvalue);

    self->x->pv[i].i = ival;
    self->x->pv[i].m = mval;
    self->x->pv[i].value = value;
    self->x->npv = i + 1;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/* TODO: This is convenient for debugging for now -- but it's not very
 * Pythonic.  It should probably be hooked into __str__ or something.
 */
static PyObject*
PyWcsprm_print_contents(PyWcsprm* self) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  if (PyWcsprm_set(self) == NULL)
    return NULL;

  PyWcsprm_Python2C(self);
  wcsprt(self->x);
  PyWcsprm_C2Python(self);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
PyWcsprm_spcfix(PyWcsprm* self) {
  int status = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  PyWcsprm_Python2C(self);
  status = spcfix(self->x);
  PyWcsprm_C2Python(self);

  if (status == -1 || status == 0) {
    return PyInt_FromLong(status);
  } else if (status > 0 && status < WCSFIX_ERRMSG_MAX) {
    PyErr_SetString(PyExc_ValueError, wcsfix_errmsg[status]);
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_sptr(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  int   i        = -1;
  char* py_ctype = 0;
  char  ctype[9];
  int   status   = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  const char* keywords[] = {"i", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", (char **)keywords,
                                   &py_ctype, &i))
    return NULL;

  if (strlen(py_ctype) > 8) {
    PyErr_SetString(PyExc_ValueError,
                    "ctype string has more than 8 characters.");
  }

  strncpy(ctype, py_ctype, 9);

  PyWcsprm_Python2C(self);
  status = wcssptr(self->x, &i, ctype);
  PyWcsprm_C2Python(self);

  if (status == 0) {
    Py_INCREF(Py_None);
    return Py_None;
  } else if (status > 0 && status < WCS_ERRMSG_MAX) {
    PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_to_header(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  int       relax   = 0;
  int       nkeyrec = 0;
  char*     header  = NULL;
  int       status  = 0;
  PyObject* result  = NULL;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  const char* keywords[] = {"relax", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i", (char **)keywords,
                                   &relax))
    return NULL;

  if (relax)
    relax = -1;

  PyWcsprm_Python2C(self);
  status = wcshdo(relax, self->x, &nkeyrec, &header);
  PyWcsprm_C2Python(self);

  if (status != 0) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    free(header);
    return NULL;
  }

  /* Just return the raw header string.  PyFITS on the Python side will help
     to parse and use this information. */
  result = PyString_FromStringAndSize(header, nkeyrec * 80);

  free(header);

  return result;
}

static PyObject*
PyWcsprm_unitfix(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  char* translate_units = NULL;
  int   ctrl            = 0;
  int   status          = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  const char* keywords[] = {"translate_units", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|s", (char **)keywords,
                                   &translate_units))
    return NULL;

  if (translate_units != NULL)
    ctrl = parse_unsafe_unit_conversion_spec(translate_units);

  status = unitfix(ctrl, self->x);

  if (status == -1 || status == 0) {
    return PyInt_FromLong(status);
  } else if (status > 0 && status < WCSFIX_ERRMSG_MAX) {
    PyErr_SetString(PyExc_ValueError, wcsfix_errmsg[status]);
    return NULL;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}


/***************************************************************************
 * Member getters
 */
static PyObject*
PyWcsprm_get_alt(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->alt == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  /* Force a null-termination of this single-character string */
  self->x->alt[1] = 0;
  return PyString_FromString(self->x->alt);
}

static int
PyWcsprm_set_alt(PyWcsprm* self, PyObject* value, void* closure) {
  char* value_string = NULL;

  if (self->x == NULL || self->x->alt == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x->alt[0] = ' ';
    self->x->alt[1] = 0;
    return 0;
  }

  value_string = PyString_AsString(value);
  if (value_string == NULL)
    return -1;

  if (!is_valid_alt_key(value_string))
    return -1;

  self->x->alt[0] = value_string[0];
  self->x->alt[1] = 0;

  return 0;
}

static PyObject*
PyWcsprm_get_cd(PyWcsprm* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x == NULL || self->x->cd == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if ((self->x->altlin & has_cd) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No cd is present.");
    return NULL;
  }

  return PyWcsprmArrayProxy_New(self, 2, dims, PyArray_DOUBLE, self->x->cd);
}

static int
PyWcsprm_set_cd(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->cd == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    self->x->altlin &= ~has_cd;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_DOUBLE,
                                                          2, 2);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2 || PyArray_DIM(value_array, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "cd must be a 2x2 array");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->cd);
  self->x->altlin |= has_cd;

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_cname(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->cname == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyWcsprmListProxy_New(self, self->x->naxis, self->x->cname);
}

static int
PyWcsprm_set_cname(PyWcsprm* self, PyObject* value, void* closure) {
  PyObject*  str      = NULL;
  char*      str_char = NULL;
  Py_ssize_t str_len  = 0;
  int        i        = 0;

  if (self->x == NULL || self->x->cname == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (!PySequence_Check(value)) {
    PyErr_SetString(PyExc_TypeError, "cname must be a sequence of strings");
    return -1;
  }

  if (PySequence_Size(value) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(cname) != naxis");
    return -1;
  }

  for (i = 0; i < self->x->naxis; ++i) {
    str = PySequence_GetItem(value, i);
    if (str == NULL)
      return -1;

    if (PyString_AsStringAndSize(str, &str_char, &str_len))
      return -1;

    if (str_len > 68) {
      PyErr_SetString(PyExc_ValueError,
                      "Each string must be less than 68 characters");
      Py_DECREF(str);
      return -1;
    }

    strncpy(self->x->cname[i], str_char, 72);
  }

  return 0;
}

static PyObject*
PyWcsprm_get_cdelt(PyWcsprm* self, void* closure) {
  Py_ssize_t naxis = 0;

  if (self->x == NULL || self->x->cdelt == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  naxis = self->x->naxis;

  return PyWcsprmArrayProxy_New(self, 1, &naxis, PyArray_DOUBLE,
                                self->x->cdelt);
}

static int
PyWcsprm_set_cdelt(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->cdelt == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the cdelt attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);

  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(cdelt) != naxis");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->cdelt);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_cel_offset(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyBool_FromLong(self->x->cel.offset);
}

static int
PyWcsprm_set_cel_offset(PyWcsprm* self, PyObject* value, void* closure) {
  long value_int;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the cel_offset attribute");
    return -1;
  }

  value_int = PyInt_AsLong(value);
  if (value_int == -1 && PyErr_Occurred())
    return -1;

  self->x->cel.offset = value_int ? 1 : 0;

  return 0;
}

static PyObject*
PyWcsprm_get_colnum(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyInt_FromLong(self->x->colnum);
}

static int
PyWcsprm_set_colnum(PyWcsprm* self, PyObject* value, void* closure) {
  long value_int;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the colnum attribute");
    return -1;
  }

  value_int = PyInt_AsLong(value);
  if (value_int == -1 && PyErr_Occurred())
    return -1;

  self->x->colnum = value_int;

  return 0;
}

static PyObject*
PyWcsprm_get_colax(PyWcsprm* self, void* closure) {
  Py_ssize_t naxis = 0;

  if (self->x == NULL || self->x->colax == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  naxis = self->x->naxis;

  return PyWcsprmArrayProxy_New(self, 1, &naxis, PyArray_INT,
                                self->x->colax);
}

static int
PyWcsprm_set_colax(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->colax == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the colax attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(colax) != naxis");
    return -1;
  }

  copy_array_to_c_int(value_array, self->x->colax);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_crder(PyWcsprm* self, void* closure) {
  Py_ssize_t naxis = 0;

  if (self->x == NULL || self->x->crder == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  naxis = self->x->naxis;

  return PyWcsprmArrayProxy_New(self, 1, &naxis, PyArray_DOUBLE,
                                self->x->crder);
}

static int
PyWcsprm_set_crder(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->crder == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the crder attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(crder) != naxis");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->crder);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_crota(PyWcsprm* self, void* closure) {
  Py_ssize_t naxis = 0;

  if (self->x == NULL || self->x->crota == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if ((self->x->altlin & has_crota) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No crota is present.");
    return NULL;
  }

  naxis = self->x->naxis;

  return PyWcsprmArrayProxy_New(self, 1, &naxis, PyArray_DOUBLE,
                                self->x->crota);
}

static int
PyWcsprm_set_crota(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->crota == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* Deletion */
    self->x->altlin &= ~has_crota;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(crota) != naxis");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->crota);
  self->x->altlin |= has_crota;

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_crpix(PyWcsprm* self, void* closure) {
  Py_ssize_t naxis = 0;

  if (self->x == NULL || self->x->crpix == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  naxis = self->x->naxis;

  return PyWcsprmArrayProxy_New(self, 1, &naxis, PyArray_DOUBLE,
                                self->x->crpix);
}

static int
PyWcsprm_set_crpix(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->crpix == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the crpix attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(crpix) != naxis");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->crpix);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_crval(PyWcsprm* self, void* closure) {
  Py_ssize_t naxis = 0;

  if (self->x == NULL || self->x->crval == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  naxis = self->x->naxis;

  return PyWcsprmArrayProxy_New(self, 1, &naxis, PyArray_DOUBLE,
                                self->x->crval);
}

static int
PyWcsprm_set_crval(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->crval == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the crval attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(crval) != naxis");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->crval);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_csyer(PyWcsprm* self, void* closure) {
  Py_ssize_t naxis;

  if (self->x == NULL || self->x->csyer == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  naxis = self->x->naxis;

  return PyWcsprmArrayProxy_New(self, 1, &naxis, PyArray_DOUBLE,
                                self->x->csyer);
}

static int
PyWcsprm_set_csyer(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->csyer == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the csyer attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(csyer) != naxis");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->csyer);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_ctype(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->ctype == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyWcsprmListProxy_New(self, self->x->naxis, self->x->ctype);
}

static int
PyWcsprm_set_ctype(PyWcsprm* self, PyObject* value, void* closure) {
  PyObject*  str      = NULL;
  char*      str_char = NULL;
  Py_ssize_t str_len  = 0;
  int        i        = 0;

  if (self->x == NULL || self->x->ctype == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (!PySequence_Check(value)) {
    PyErr_SetString(PyExc_TypeError, "ctype must be a sequence of strings");
    return -1;
  }

  if (PySequence_Size(value) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(ctype) != naxis");
    return -1;
  }

  for (i = 0; i < self->x->naxis; ++i) {
    str = PySequence_GetItem(value, i);
    if (str == NULL)
      return -1;

    if (PyString_AsStringAndSize(str, &str_char, &str_len))
      return -1;

    if (str_len > 68) {
      PyErr_SetString(PyExc_ValueError,
                      "Each string must be less than 68 characters");
      Py_DECREF(str);
      return -1;
    }

    strncpy(self->x->ctype[i], str_char, 72);

    Py_DECREF(str);
  }

  return 0;
}

static PyObject*
PyWcsprm_get_cubeface(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyBool_FromLong(self->x->cubeface);
}

static int
PyWcsprm_set_cubeface(PyWcsprm* self, PyObject* value, void* closure) {
  long value_bool;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the cubeface attribute");
    return -1;
  }

  value_bool = PyInt_AsLong(value);
  if (value_bool == -1 && PyErr_Occurred())
    return -1;

  self->x->cubeface = value_bool ? 1 : 0;

  return 0;
}

static PyObject*
PyWcsprm_get_cunit(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->cunit == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyWcsprmListProxy_New(self, self->x->naxis, self->x->cunit);
}

static int
PyWcsprm_set_cunit(PyWcsprm* self, PyObject* value, void* closure) {
  PyObject*  str      = NULL;
  char*      str_char = NULL;
  Py_ssize_t str_len  = 0;
  int        i        = 0;

  if (self->x == NULL || self->x->cunit == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (!PySequence_Check(value)) {
    PyErr_SetString(PyExc_TypeError, "cunit must be a sequence of strings");
    return -1;
  }

  if (PySequence_Size(value) != self->x->naxis) {
    PyErr_SetString(PyExc_ValueError, "len(cunit) != naxis");
    return -1;
  }

  for (i = 0; i < self->x->naxis; ++i) {
    str = PySequence_GetItem(value, i);
    if (str == NULL)
      return -1;

    if (PyString_AsStringAndSize(str, &str_char, &str_len))
      return -1;

    if (str_len > 68) {
      PyErr_SetString(PyExc_ValueError,
                      "Each string must be less than 68 characters");
      Py_DECREF(str);
      return -1;
    }

    strncpy(self->x->cunit[i], str_char, 72);

    Py_DECREF(str);
  }

  return 0;
}

static PyObject*
PyWcsprm_get_dateavg(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->dateavg == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyString_FromString(self->x->dateavg);
}

static int
PyWcsprm_set_dateavg(PyWcsprm* self, PyObject* value, void* closure) {
  char*      value_string = NULL;
  Py_ssize_t value_length = 0;

  if (self->x == NULL || self->x->dateavg == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Can not delete the dateavg attribute.");
    return -1;
  }

  if (PyString_AsStringAndSize(value, &value_string, &value_length))
    return -1;

  if (value_length >= 72) {
    PyErr_SetString(PyExc_ValueError, "len(dateavg) must be less than 72");
    return -1;
  }

  strncpy(self->x->dateavg, value_string, 72);

  return 0;
}

static PyObject*
PyWcsprm_get_dateobs(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->dateobs == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyString_FromString(self->x->dateobs);
}

static int
PyWcsprm_set_dateobs(PyWcsprm* self, PyObject* value, void* closure) {
  char*      value_string = NULL;
  Py_ssize_t value_length = 0;

  if (self->x == NULL || self->x->dateobs == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Can not delete the dateobs attribute.");
    return -1;
  }

  if (PyString_AsStringAndSize(value, &value_string, &value_length))
    return -1;

  if (value_length >= 72) {
    PyErr_SetString(PyExc_ValueError, "len(dateobs) must be less than 72");
    return -1;
  }

  strncpy(self->x->dateobs, value_string, 72);

  return 0;
}

static PyObject*
PyWcsprm_get_equinox(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->equinox);
}

static int
PyWcsprm_set_equinox(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x->equinox = UNDEFINED;
    return 0;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->equinox = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_imgpix_matrix(PyWcsprm* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x == NULL || self->x->lin.imgpix == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyWcsprmArrayProxy_New(self, 2, dims, PyArray_DOUBLE,
                                self->x->lin.imgpix);
}

static PyObject*
PyWcsprm_get_lat(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyInt_FromLong(self->x->lat);
}

static PyObject*
PyWcsprm_get_latpole(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyFloat_FromDouble(self->x->latpole);
}

static int
PyWcsprm_set_latpole(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    PyErr_SetString(PyExc_TypeError, "Can not delete the latpole attribute");
    return -1;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->latpole = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_lattyp(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->lattyp == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyString_FromString(self->x->lattyp);
}

static PyObject*
PyWcsprm_get_lng(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyInt_FromLong(self->x->lng);
}

static PyObject*
PyWcsprm_get_lngtyp(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->lngtyp == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyString_FromString(self->x->lngtyp);
}

static PyObject*
PyWcsprm_get_lonpole(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyFloat_FromDouble(self->x->lonpole);
}

static int
PyWcsprm_set_lonpole(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    PyErr_SetString(PyExc_TypeError, "Can not delete the lonpole attribute");
    return -1;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->lonpole = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_mjdavg(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->mjdavg);
}

static int
PyWcsprm_set_mjdavg(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    PyErr_SetString(PyExc_TypeError, "Can not delete the mjdavg attribute");
    return -1;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->mjdavg = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_mjdobs(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->mjdobs);
}

static int
PyWcsprm_set_mjdobs(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    PyErr_SetString(PyExc_TypeError, "Can not delete the mjdobs attribute");
    return -1;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->mjdobs = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_name(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->wcsname == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyString_FromString(self->x->wcsname);
}

static int
PyWcsprm_set_name(PyWcsprm* self, PyObject* value, void* closure) {
  char*      value_string = NULL;
  Py_ssize_t value_length = 0;

  if (self->x == NULL || self->x->wcsname == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Can not delete the name attribute");
    return -1;
  }

  if (PyString_AsStringAndSize(value, &value_string, &value_length))
    return -1;

  if (value_length >= 72) {
    PyErr_SetString(PyExc_ValueError, "len(name) must be less than 72");
    return -1;
  }

  strncpy(self->x->wcsname, value_string, 72);

  return 0;
}

static PyObject*
PyWcsprm_get_naxis(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyInt_FromLong(self->x->naxis);
}

static PyObject*
PyWcsprm_get_obsgeo(PyWcsprm* self, void* closure) {
  Py_ssize_t size = 3;

  if (self->x == NULL || self->x->obsgeo == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyWcsprmArrayProxy_New(self, 1, &size, PyArray_DOUBLE,
                                self->x->obsgeo);
}

static int
PyWcsprm_set_obsgeo(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->obsgeo == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the obsgeo attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 3) {
    PyErr_SetString(PyExc_ValueError, "len(obsgeo) != 3");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->obsgeo);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_pc(PyWcsprm* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x == NULL || self->x->pc == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if ((self->x->altlin & has_pc) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No pc is present.");
    return NULL;
  }

  return PyWcsprmArrayProxy_New(self, 2, dims, PyArray_DOUBLE, self->x->pc);
}

static int
PyWcsprm_set_pc(PyWcsprm* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x == NULL || self->x->pc == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x->altlin &= ~has_pc;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          2, 2);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2 || PyArray_DIM(value_array, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "pc must be a 2x2 array");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x->pc);
  self->x->altlin |= has_pc;

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyWcsprm_get_phi0(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->cel.phi0);
}

static int
PyWcsprm_set_phi0(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    PyErr_SetString(PyExc_TypeError, "Can not delete the phi0 attribute");
    return -1;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->cel.phi0 = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_piximg_matrix(PyWcsprm* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x == NULL || self->x->lin.piximg == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyWcsprmArrayProxy_New(self, 2, dims, PyArray_DOUBLE,
                                self->x->lin.piximg);
}

static PyObject*
PyWcsprm_get_radesys(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->radesys == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyString_FromString(self->x->radesys);
}

static int
PyWcsprm_set_radesys(PyWcsprm* self, PyObject* value, void* closure) {
  char*      value_string = NULL;
  Py_ssize_t value_length = 0;

  if (self->x == NULL || self->x->radesys == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Can not delete the radesys attribute");
    return -1;
  }

  if (PyString_AsStringAndSize(value, &value_string, &value_length))
    return -1;

  if (value_length >= 72) {
    PyErr_SetString(PyExc_ValueError, "len(name) must be less than 72");
    return -1;
  }

  strncpy(self->x->radesys, value_string, 72);

  return 0;
}

static PyObject*
PyWcsprm_get_restfrq(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->restfrq);
}

static int
PyWcsprm_set_restfrq(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x->restfrq = UNDEFINED;
    return 0;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->restfrq = PyFloat_AsDouble(value);;

  return 0;
}

static PyObject*
PyWcsprm_get_restwav(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->restwav);
}

static int
PyWcsprm_set_restwav(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x->restwav = UNDEFINED;
    return 0;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->restwav = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_spec(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyInt_FromLong(self->x->spec);
}

static PyObject*
PyWcsprm_get_specsys(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->specsys == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyString_FromString(self->x->specsys);
}

static int
PyWcsprm_set_specsys(PyWcsprm* self, PyObject* value, void* closure) {
  char*      value_string = NULL;
  Py_ssize_t value_length = 0;

  if (self->x == NULL || self->x->specsys == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Can not delete the specsys attribute");
    return -1;
  }

  if (PyString_AsStringAndSize(value, &value_string, &value_length))
    return -1;

  if (value_length >= 72) {
    PyErr_SetString(PyExc_ValueError, "len(name) must be less than 72");
    return -1;
  }

  strncpy(self->x->specsys, value_string, 72);

  return 0;
}

static PyObject*
PyWcsprm_get_ssysobs(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->ssysobs == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyString_FromString(self->x->ssysobs);
}

static int
PyWcsprm_set_ssysobs(PyWcsprm* self, PyObject* value, void* closure) {
  char*      value_string = NULL;
  Py_ssize_t value_length = 0;

  if (self->x == NULL || self->x->ssysobs == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Can not delete the ssysobs attribute");
    return -1;
  }

  if (PyString_AsStringAndSize(value, &value_string, &value_length))
    return -1;

  if (value_length >= 72) {
    PyErr_SetString(PyExc_ValueError, "len(name) must be less than 72");
    return -1;
  }

  strncpy(self->x->ssysobs, value_string, 72);

  return 0;
}

static PyObject*
PyWcsprm_get_ssyssrc(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->ssyssrc == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyString_FromString(self->x->ssyssrc);
}

static int
PyWcsprm_set_ssyssrc(PyWcsprm* self, PyObject* value, void* closure) {
  char*      value_string = NULL;
  Py_ssize_t value_length = 0;

  if (self->x == NULL || self->x->ssyssrc == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Can not delete the ssyssrc attribute");
    return -1;
  }

  if (PyString_AsStringAndSize(value, &value_string, &value_length))
    return -1;

  if (value_length >= 72) {
    PyErr_SetString(PyExc_ValueError, "len(name) must be less than 72");
    return -1;
  }

  strncpy(self->x->ssyssrc, value_string, 72);

  return 0;
}

static PyObject*
PyWcsprm_get_theta0(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->cel.theta0);
}

static int
PyWcsprm_set_theta0(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    PyErr_SetString(PyExc_TypeError, "Can not delete the theta0 attribute");
    return -1;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->cel.theta0 = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_velangl(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->velangl);
}

static int
PyWcsprm_set_velangl(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x->velangl = UNDEFINED;
    return 0;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->velangl = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_velosys(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->velosys);
}

static int
PyWcsprm_set_velosys(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x->velosys = UNDEFINED;
    return 0;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->velosys = PyFloat_AsDouble(value);

  return 0;
}

static PyObject*
PyWcsprm_get_zsource(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyFloat_FromDouble(self->x->zsource);
}

static int
PyWcsprm_set_zsource(PyWcsprm* self, PyObject* value, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x->zsource = UNDEFINED;
    return 0;
  }

  if (!PyFloat_Check(value))
    return -1;

  self->x->zsource = PyFloat_AsDouble(value);

  return 0;
}

/***************************************************************************
 * PyWcsprm definition structures
 */

static PyMemberDef PyWcsprm_members[] = {
  {NULL}
};

static PyGetSetDef PyWcsprm_getset[] = {
  {"alt", (getter)PyWcsprm_get_alt, (setter)PyWcsprm_set_alt, (char *)doc_alt},
  {"cd", (getter)PyWcsprm_get_cd, (setter)PyWcsprm_set_cd, (char *)doc_cd},
  {"cdelt", (getter)PyWcsprm_get_cdelt, (setter)PyWcsprm_set_cdelt, (char *)doc_cdelt},
  {"cel_offset", (getter)PyWcsprm_get_cel_offset, (setter)PyWcsprm_set_cel_offset, (char *)doc_cel_offset},
  {"cname", (getter)PyWcsprm_get_cname, (setter)PyWcsprm_set_cname, (char *)doc_cname},
  {"colax", (getter)PyWcsprm_get_colax, (setter)PyWcsprm_set_colax, (char *)doc_colax},
  {"colnum", (getter)PyWcsprm_get_colnum, (setter)PyWcsprm_set_colnum, (char *)doc_colnum},
  {"crder", (getter)PyWcsprm_get_crder, (setter)PyWcsprm_set_crder, (char *)doc_crder},
  {"crota", (getter)PyWcsprm_get_crota, (setter)PyWcsprm_set_crota, (char *)doc_crota},
  {"crpix", (getter)PyWcsprm_get_crpix, (setter)PyWcsprm_set_crpix, (char *)doc_crpix},
  {"crval", (getter)PyWcsprm_get_crval, (setter)PyWcsprm_set_crval, (char *)doc_crval},
  {"csyer", (getter)PyWcsprm_get_csyer, (setter)PyWcsprm_set_csyer, (char *)doc_csyer},
  {"ctype", (getter)PyWcsprm_get_ctype, (setter)PyWcsprm_set_ctype, (char *)doc_ctype},
  {"cubeface", (getter)PyWcsprm_get_cubeface, (setter)PyWcsprm_set_cubeface, (char *)doc_cubeface},
  {"cunit", (getter)PyWcsprm_get_cunit, (setter)PyWcsprm_set_cunit, (char *)doc_cunit},
  {"dateavg", (getter)PyWcsprm_get_dateavg, (setter)PyWcsprm_set_dateavg, (char *)doc_dateavg},
  {"dateobs", (getter)PyWcsprm_get_dateobs, (setter)PyWcsprm_set_dateobs, (char *)doc_dateobs},
  {"equinox", (getter)PyWcsprm_get_equinox, (setter)PyWcsprm_set_equinox, (char *)doc_equinox},
  {"imgpix_matrix", (getter)PyWcsprm_get_imgpix_matrix, NULL, (char *)doc_imgpix_matrix},
  {"lat", (getter)PyWcsprm_get_lat, NULL, (char *)doc_lat},
  {"latpole", (getter)PyWcsprm_get_latpole, (setter)PyWcsprm_set_latpole, (char *)doc_latpole},
  {"lattyp", (getter)PyWcsprm_get_lattyp, NULL, (char *)doc_lattyp},
  {"lng", (getter)PyWcsprm_get_lng, NULL, (char *)doc_lng},
  {"lngtyp", (getter)PyWcsprm_get_lngtyp, NULL, (char *)doc_lngtyp},
  {"lonpole", (getter)PyWcsprm_get_lonpole, (setter)PyWcsprm_set_lonpole, (char *)doc_lonpole},
  {"mjdavg", (getter)PyWcsprm_get_mjdavg, (setter)PyWcsprm_set_mjdavg, (char *)doc_mjdavg},
  {"mjdobs", (getter)PyWcsprm_get_mjdobs, (setter)PyWcsprm_set_mjdobs, (char *)doc_mjdobs},
  {"name", (getter)PyWcsprm_get_name, (setter)PyWcsprm_set_name, (char *)doc_name},
  {"naxis", (getter)PyWcsprm_get_naxis, NULL, (char *)doc_naxis},
  {"obsgeo", (getter)PyWcsprm_get_obsgeo, (setter)PyWcsprm_set_obsgeo, (char *)doc_obsgeo},
  {"pc", (getter)PyWcsprm_get_pc, (setter)PyWcsprm_set_pc, (char *)doc_pc},
  {"phi0", (getter)PyWcsprm_get_phi0, (setter)PyWcsprm_set_phi0, (char *)doc_phi0},
  {"piximg_matrix", (getter)PyWcsprm_get_piximg_matrix, NULL, (char *)doc_piximg_matrix},
  {"radesys", (getter)PyWcsprm_get_radesys, (setter)PyWcsprm_set_radesys, (char *)doc_radesys},
  {"restfrq", (getter)PyWcsprm_get_restfrq, (setter)PyWcsprm_set_restfrq, (char *)doc_restfrq},
  {"restwav", (getter)PyWcsprm_get_restwav, (setter)PyWcsprm_set_restwav, (char *)doc_restwav},
  {"spec", (getter)PyWcsprm_get_spec, NULL, (char *)doc_spec},
  {"specsys", (getter)PyWcsprm_get_specsys, (setter)PyWcsprm_set_specsys, (char *)doc_specsys},
  {"ssysobs", (getter)PyWcsprm_get_ssysobs, (setter)PyWcsprm_set_ssysobs, (char *)doc_ssysobs},
  {"ssyssrc", (getter)PyWcsprm_get_ssyssrc, (setter)PyWcsprm_set_ssyssrc, (char *)doc_ssyssrc},
  {"theta0", (getter)PyWcsprm_get_theta0, (setter)PyWcsprm_set_theta0, (char *)doc_theta0},
  {"velangl", (getter)PyWcsprm_get_velangl, (setter)PyWcsprm_set_velangl, (char *)doc_velangl},
  {"velosys", (getter)PyWcsprm_get_velosys, (setter)PyWcsprm_set_velosys, (char *)doc_velosys},
  {"zsource", (getter)PyWcsprm_get_zsource, (setter)PyWcsprm_set_zsource, (char *)doc_zsource},
  {NULL}
};

static PyMethodDef PyWcsprm_methods[] = {
  {"celfix", (PyCFunction)PyWcsprm_celfix, METH_NOARGS, doc_celfix},
  {"copy", (PyCFunction)PyWcsprm_copy, METH_NOARGS, doc_copy},
  {"__copy__", (PyCFunction)PyWcsprm_copy, METH_NOARGS, doc_copy},
  {"cylfix", (PyCFunction)PyWcsprm_cylfix, METH_VARARGS, doc_cylfix},
  {"datfix", (PyCFunction)PyWcsprm_datfix, METH_NOARGS, doc_datfix},
  {"fix", (PyCFunction)PyWcsprm_fix, METH_VARARGS, doc_fix},
  {"get_ps", (PyCFunction)PyWcsprm_get_ps, METH_NOARGS, doc_get_ps},
  {"get_pv", (PyCFunction)PyWcsprm_get_pv, METH_NOARGS, doc_get_pv},
  {"has_cdi_ja", (PyCFunction)PyWcsprm_has_cdi_ja, METH_NOARGS, doc_has_cdi_ja},
  {"has_crotaia", (PyCFunction)PyWcsprm_has_crotaia, METH_NOARGS, doc_has_crotaia},
  {"has_pci_ja", (PyCFunction)PyWcsprm_has_pci_ja, METH_NOARGS, doc_has_pci_ja},
  {"is_unity", (PyCFunction)PyWcsprm_is_unity, METH_NOARGS, doc_is_unity},
  {"mix", (PyCFunction)PyWcsprm_mix, METH_VARARGS, doc_mix},
  {"p2s", (PyCFunction)PyWcsprm_p2s, METH_O, doc_p2s},
  {"print_contents", (PyCFunction)PyWcsprm_print_contents, METH_NOARGS, doc_print_contents},
  {"s2p", (PyCFunction)PyWcsprm_s2p, METH_O, doc_s2p},
  {"set", (PyCFunction)PyWcsprm_set, METH_NOARGS, doc_set},
  {"set_ps", (PyCFunction)PyWcsprm_set_ps, METH_O, doc_set_ps},
  {"set_pv", (PyCFunction)PyWcsprm_set_pv, METH_O, doc_set_pv},
  {"spcfix", (PyCFunction)PyWcsprm_spcfix, METH_NOARGS, doc_spcfix},
  {"sptr", (PyCFunction)PyWcsprm_sptr, METH_NOARGS, doc_sptr},
  {"to_header", (PyCFunction)PyWcsprm_to_header, METH_VARARGS, doc_to_header},
  {"unitfix", (PyCFunction)PyWcsprm_unitfix, METH_VARARGS, doc_unitfix},
  {NULL}
};

static PyTypeObject PyWcsprmType = {
  PyObject_HEAD_INIT(NULL)
  0,                          /*ob_size*/
  "pywcs._WCS",               /*tp_name*/
  sizeof(PyWcsprm),           /*tp_basicsize*/
  0,                          /*tp_itemsize*/
  (destructor)PyWcsprm_dealloc, /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  (reprfunc)PyWcsprm_repr,    /*tp_repr*/
  0,                          /*tp_as_number*/
  0,                          /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash */
  0,                          /*tp_call*/
  (reprfunc)PyWcsprm_repr,    /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  doc_WCS,                    /* tp_doc */
  0,                          /* tp_traverse */
  0,                          /* tp_clear */
  0,                          /* tp_richcompare */
  0,                          /* tp_weaklistoffset */
  0,                          /* tp_iter */
  0,                          /* tp_iternext */
  PyWcsprm_methods,           /* tp_methods */
  PyWcsprm_members,           /* tp_members */
  PyWcsprm_getset,            /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tp_dictoffset */
  (initproc)PyWcsprm_init,    /* tp_init */
  0,                          /* tp_alloc */
  PyWcsprm_new,               /* tp_new */
};

/***************************************************************************
 * List-of-strings proxy object
 ***************************************************************************/

static PyTypeObject PyWcsprmListProxyType;

typedef struct {
  PyObject_HEAD
  PyObject* pywcsprm;
  Py_ssize_t size;
  char (*array)[72];
} PyWcsprmListProxy;

static void
PyWcsprmListProxy_dealloc(PyWcsprmListProxy* self) {
  Py_XDECREF(self->pywcsprm);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
PyWcsprmListProxy_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
  PyWcsprmListProxy* self = NULL;

  self = (PyWcsprmListProxy*)type->tp_alloc(type, 0);
  if (self != NULL)
    self->pywcsprm = NULL;
  return (PyObject*)self;
}

static int
PyWcsprmListProxy_traverse(PyWcsprmListProxy* self, visitproc visit, void *arg)
{
  int vret;

  if (self->pywcsprm) {
    vret = visit(self->pywcsprm, arg);
    if (vret != 0)
      return vret;
  }

  return 0;
}

static int
PyWcsprmListProxy_clear(PyWcsprmListProxy *self)
{
  PyObject *tmp;

  tmp = self->pywcsprm;
  self->pywcsprm = NULL;
  Py_XDECREF(tmp);

  return 0;
}

static PyObject *
PyWcsprmListProxy_New(PyWcsprm* owner, Py_ssize_t size, char (*array)[72]) {
  PyWcsprmListProxy* self = NULL;

  self = (PyWcsprmListProxy*)PyWcsprmListProxyType.tp_alloc(&PyWcsprmListProxyType, 0);
  if (self == NULL)
    return NULL;

  Py_INCREF(owner);
  self->pywcsprm = (PyObject*)owner;
  self->size = size;
  self->array = array;
  return (PyObject*)self;
}

static Py_ssize_t
PyWcsprmListProxy_len(PyWcsprmListProxy* self) {
  return self->size;
}

static PyObject*
PyWcsprmListProxy_getitem(PyWcsprmListProxy* self, Py_ssize_t index) {
  if (index > self->size) {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    return NULL;
  }

  return PyString_FromString(self->array[index]);
}

static int
PyWcsprmListProxy_setitem(PyWcsprmListProxy* self, Py_ssize_t index, PyObject* arg) {
  char* value;
  Py_ssize_t value_length;

  if (index > self->size) {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    return -1;
  }

  if (PyString_AsStringAndSize(arg, &value, &value_length))
    return -1;

  if (value_length >= 68) {
    PyErr_SetString(PyExc_ValueError, "string must be less than 68 characters");
    return -1;
  }

  strncpy(self->array[index], value, 72);

  return 0;
}

static PyObject*
PyWcsprmListProxy_repr(PyWcsprmListProxy* self) {
  char*      buffer = NULL;
  char*      wp     = NULL;
  char*      rp     = NULL;
  Py_ssize_t i      = 0;
  Py_ssize_t j      = 0;
  PyObject*  result = NULL;

  /* Overallocating, of course, to be safe */
  buffer = malloc(self->size*80 + 2);
  if (buffer == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
    return NULL;
  }

  buffer[0] = '[';
  wp = buffer + 1;

  for (i = 0; i < self->size; ++i) {
    *wp++ = '\'';
    rp = self->array[i];
    for (j = 0; j < 68 && *rp != 0; ++j)
      *wp++ = *rp++;
    *wp++ = '\'';

    if (i != self->size - 1) {
      *wp++ = ',';
      *wp++ = ' ';
    }
  }

  *wp++ = ']';
  *wp++ = 0;

  result = PyString_FromString(buffer);
  free(buffer);
  return result;
}

static PySequenceMethods PyWcsprmListProxy_sequence_methods = {
  (lenfunc)PyWcsprmListProxy_len,
  NULL,
  NULL,
  (ssizeargfunc)PyWcsprmListProxy_getitem,
  NULL,
  (ssizeobjargproc)PyWcsprmListProxy_setitem,
  NULL,
  NULL,
  NULL
};

static PyTypeObject PyWcsprmListProxyType = {
  PyObject_HEAD_INIT(NULL)
  0,                          /*ob_size*/
  "pywcs.ListProxy",          /*tp_name*/
  sizeof(PyWcsprmListProxy),  /*tp_basicsize*/
  0,                          /*tp_itemsize*/
  (destructor)PyWcsprmListProxy_dealloc, /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  (reprfunc)PyWcsprmListProxy_repr, /*tp_repr*/
  0,                          /*tp_as_number*/
  &PyWcsprmListProxy_sequence_methods, /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash */
  0,                          /*tp_call*/
  (reprfunc)PyWcsprmListProxy_repr, /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
  0,                          /* tp_doc */
  (traverseproc)PyWcsprmListProxy_traverse, /* tp_traverse */
  (inquiry)PyWcsprmListProxy_clear, /* tp_clear */
  0,                          /* tp_richcompare */
  0,                          /* tp_weaklistoffset */
  0,                          /* tp_iter */
  0,                          /* tp_iternext */
  0,                          /* tp_methods */
  0,                          /* tp_members */
  0,                          /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tp_dictoffset */
  0,                          /* tp_init */
  0,                          /* tp_alloc */
  PyWcsprmListProxy_new,      /* tp_new */
};

/***************************************************************************
 * Module-level
 ***************************************************************************/

static PyMethodDef module_methods[] = {
  {NULL}
};

#define DEFINE_EXCEPTION(exc) \
  WcsExc_##exc = PyErr_NewException("_pywcs." #exc "Error", PyExc_ValueError, NULL); \
  if (WcsExc_##exc == NULL) \
    return; \
  PyModule_AddObject(m, #exc "Error", WcsExc_##exc); \


#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_pywcs(void)
{
  PyObject* m;

  import_array();

  if (PyType_Ready(&PyWcsprmType) < 0)
    return;

  if (PyType_Ready(&PyWcsprmListProxyType) < 0)
    return;

  m = Py_InitModule3("_pywcs", module_methods, doc_pywcs);

  if (m == NULL)
    return;

  Py_INCREF(&PyWcsprmType);
  PyModule_AddObject(m, "_WCS", (PyObject *)&PyWcsprmType);

  DEFINE_EXCEPTION(SingularMatrix);
  DEFINE_EXCEPTION(InconsistentAxisTypes);
  DEFINE_EXCEPTION(InvalidTransform);
  DEFINE_EXCEPTION(InvalidCoordinate);
  DEFINE_EXCEPTION(NoSolution);
  DEFINE_EXCEPTION(InvalidSubimageSpecification);
  DEFINE_EXCEPTION(NonseparableSubimageCoordinateSystem);

  PyModule_AddObject(m, "__docformat__", PyString_FromString("epytext"));
}
