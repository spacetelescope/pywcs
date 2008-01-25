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
#include <wcslib/wcs.h>
#include <wcslib/wcsfix.h>
#include <wcslib/wcshdr.h>

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

PyObject*
PyArray_SimpleNewFromDataCopy(int nd, const npy_intp* dims, int typenum,
                              const void* data) {
  PyArrayObject* result = NULL;
  int size = 0;
  int i = 0;

  result = (PyArrayObject*) PyArray_SimpleNew(nd, (npy_intp*)dims, typenum);
  if (result == NULL)
    return NULL;

  size = PyArray_ITEMSIZE(result);
  for (i = 0; i < PyArray_NDIM(result); ++i)
    size *= PyArray_DIM(result, i);

  memcpy(PyArray_DATA(result), data, size);

  return (PyObject*)result;
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

/* This is an array mapping the wcs status codes to Python
 * exception types.  The exception string is stored as part
 * of wcslib itself in wcs_errmsg.
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
  return (PyObject*)PyWcsprm_cnew(NULL);
}

static int
PyWcsprm_init(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  int            naxis;
  struct wcsprm* obj;
  int            status;

  if (self->x != NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Wcs object is not new during initialization");
    return -1;
  }

  if (!PyArg_ParseTuple(args, "i", &naxis))
    return -1;

  obj = malloc(sizeof(struct wcsprm));
  if (obj == NULL) {
    PyErr_SetString(PyExc_MemoryError,
                    "Could not allocate wcsprm object");
    return -1;
  }

  status = wcsini(1, naxis, obj);

  if (status != 0) {
    free(obj);
    PyErr_SetString(PyExc_MemoryError,
                    "Could not initialize wcsprm object");
    return -1;
  }

  self->x = obj;

  return 0;
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

  status = celfix(self->x);

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
    naxis_array = (PyArrayObject*)PyArray_ContiguousFromObject(naxis_obj, 1, 1,
                                                               PyArray_INT);
    if (naxis_array == NULL)
      return NULL;
    if (PyArray_DIM(naxis_array, 0) != self->x->naxis) {
      PyErr_Format(PyExc_ValueError,
                   "naxis must be same length as the number of axes of "
                   "the Wcs object (%d).",
                   self->x->naxis);
      Py_DECREF(naxis_array);
      return NULL;
    }
    naxis = (int*)PyArray_DATA(naxis_array);
  }

  status = cylfix(naxis, self->x);

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
    naxis_array = (PyArrayObject*)PyArray_ContiguousFromObject(naxis_obj, 1, 1,
                                                               PyArray_INT);
    if (naxis_array == NULL)
      return NULL;
    if (PyArray_DIM(naxis_array, 0) != self->x->naxis) {
      PyErr_Format(PyExc_ValueError,
                   "naxis must be same length as the number of axes of "
                   "the Wcs object (%d).",
                   self->x->naxis);
      Py_DECREF(naxis_array);
      return NULL;
    }
    naxis = (int*)PyArray_DATA(naxis_array);
  }

  status = wcsfix(ctrl, naxis, self->x, stat);
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
    }
  }

  return result;
}

static PyObject*
PyWcsprm_has_cdi_ja(PyWcsprm* self) {
  int result = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  result = self->x->altlin & 2;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_has_crotaia(PyWcsprm* self) {
  int result = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  result = self->x->altlin & 4;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_has_pci_ja(PyWcsprm* self) {
  int result = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  result = self->x->altlin & 1;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_mix(PyWcsprm* self, PyObject* args, PyObject* kwds) {
  int            mixpix     = 0;
  int            mixcel     = 0;
  double         vspan[2]   = {0, 0};
  double         vstep      = 0;
  int            viter      = 0;
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

  world = (PyArrayObject*)PyArray_ContiguousFromObject
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

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromObject
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
  phi = (PyArrayObject*)PyArray_SimpleNew
    (1, &self->x->naxis, PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto __PyWcsprm_mix_exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, &self->x->naxis, PyArray_DOUBLE);
  if (theta == NULL) {
    status = 2;
    goto __PyWcsprm_mix_exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (1, &self->x->naxis, PyArray_DOUBLE);
  if (imgcrd == NULL) {
    status = 2;
    goto __PyWcsprm_mix_exit;
  }

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
  PyArrayObject* pixcrd = NULL;
  PyArrayObject* imgcrd = NULL;
  PyArrayObject* phi    = NULL;
  PyArrayObject* theta  = NULL;
  PyArrayObject* world  = NULL;
  PyArrayObject* stat   = NULL;
  PyObject*      result = NULL;
  int            status = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromObject
    (arg, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL)
    return NULL;

  /* Now we allocate a bunch of numpy arrays to store the
   * results in.
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

  /* Make the call */
  status = wcsp2s(self->x,
                  PyArray_DIM(pixcrd, 0),
                  PyArray_DIM(pixcrd, 1),
                  (double*)PyArray_DATA(pixcrd),
                  (double*)PyArray_DATA(imgcrd),
                  (double*)PyArray_DATA(phi),
                  (double*)PyArray_DATA(theta),
                  (double*)PyArray_DATA(world),
                  (int*)PyArray_DATA(stat));

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
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return NULL;
  }
}

static PyObject*
PyWcsprm_s2p(PyWcsprm* self, PyObject* arg) {
  PyArrayObject* world  = NULL;
  PyArrayObject* phi    = NULL;
  PyArrayObject* theta  = NULL;
  PyArrayObject* imgcrd = NULL;
  PyArrayObject* pixcrd = NULL;
  PyArrayObject* stat   = NULL;

  PyObject* result = NULL;
  int status = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError,
                    "Underlying object is NULL.");
    return NULL;
  }

  world = (PyArrayObject*)PyArray_ContiguousFromObject
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
  status = wcss2p(self->x,
                  PyArray_DIM(world, 0),
                  PyArray_DIM(world, 1),
                  (double*)PyArray_DATA(world),
                  (double*)PyArray_DATA(phi),
                  (double*)PyArray_DATA(theta),
                  (double*)PyArray_DATA(imgcrd),
                  (double*)PyArray_DATA(pixcrd),
                  (int*)PyArray_DATA(stat));

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

static PyObject*
PyWcsprm_set(PyWcsprm* self) {
  int status = 0;
  status = wcsset(self->x);

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

/* TODO: This is convenient for debugging for now -- but it's
 * not very Pythonic.  It should probably be hooked into __str__
 * or something.
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
  wcsprt(self->x);
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

  status = spcfix(self->x);

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
  int   i = -1;
  char* py_ctype;
  char  ctype[9];
  int   status;

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

  status = wcssptr(self->x, &i, ctype);

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
PyWcsprm_get_cd(PyWcsprm* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if (self->x->cd == NULL || (self->x->altlin & 2) == 0) {
    PyErr_SetString(PyExc_AssertionError, "No CDi_ja is present.");
    return NULL;
  } else {
    return PyArray_SimpleNewFromDataCopy(2, dims,
                                         PyArray_DOUBLE, self->x->cd);
  }
}

static PyObject*
PyWcsprm_get_cdelt(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->cdelt == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyArray_SimpleNewFromDataCopy(1, &self->x->naxis,
                                       PyArray_DOUBLE, self->x->cdelt);
}

static PyObject*
PyWcsprm_get_crota(PyWcsprm* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  if (self->x->crota == NULL || (self->x->altlin & 4) == 0) {
    PyErr_SetString(PyExc_AssertionError, "No CROTAia is present.");
    return NULL;
  } else {
    return PyArray_SimpleNewFromDataCopy(2, dims,
                                         PyArray_DOUBLE, self->x->crota);
  }
}

static PyObject*
PyWcsprm_get_crpix(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->crpix == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyArray_SimpleNewFromDataCopy(1, &self->x->naxis,
                                       PyArray_DOUBLE, self->x->crpix);
}

static PyObject*
PyWcsprm_get_crval(PyWcsprm* self, void* closure) {
  if (self->x == NULL || self->x->crval == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyArray_SimpleNewFromDataCopy(1, &self->x->naxis,
                                       PyArray_DOUBLE, self->x->crval);
}

static PyObject*
PyWcsprm_get_ctype(PyWcsprm* self, void* closure) {
  PyObject* result;
  PyObject* subresult;
  int i = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  result = PyList_New(self->x->naxis);
  if (result == NULL)
    return NULL;

  for (i = 0; i < self->x->naxis; ++i) {
    subresult = PyString_FromString(self->x->ctype[i]);
    if (subresult == NULL ||
        PyList_SetItem(result, i, subresult)) {
      Py_XDECREF(subresult);
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

static PyObject*
PyWcsprm_get_cubeface(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyBool_FromLong(self->x->cubeface);
}

static PyObject*
PyWcsprm_get_cunit(PyWcsprm* self, void* closure) {
  PyObject* result;
  PyObject* subresult;
  int i = 0;

  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  result = PyList_New(self->x->naxis);
  if (result == NULL)
    return NULL;

  for (i = 0; i < self->x->naxis; ++i) {
    subresult = PyString_FromString(self->x->cunit[i]);
    if (subresult == NULL ||
        PyList_SetItem(result, i, subresult)) {
      Py_XDECREF(subresult);
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
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

static PyObject*
PyWcsprm_get_lng(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyInt_FromLong(self->x->lng);
}

static PyObject*
PyWcsprm_get_lonpole(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyFloat_FromDouble(self->x->lonpole);
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
PyWcsprm_get_pc(PyWcsprm* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x == NULL || self->x->pc == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyArray_SimpleNewFromDataCopy(2, dims,
                                       PyArray_DOUBLE, self->x->pc);
}

static PyObject*
PyWcsprm_get_ps(PyWcsprm* self, void* closure) {
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
PyWcsprm_get_pv(PyWcsprm* self, void* closure) {
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

static PyObject*
PyWcsprm_get_restfrq(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyFloat_FromDouble(self->x->restfrq);
}

static PyObject*
PyWcsprm_get_restwav(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyFloat_FromDouble(self->x->restwav);
}

static PyObject*
PyWcsprm_get_spec(PyWcsprm* self, void* closure) {
  if (self->x == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }
  return PyInt_FromLong(self->x->spec);
}

/***************************************************************************
 * PyWcsprm definition structures
 */

static PyMemberDef PyWcsprm_members[] = {
  {NULL}
};

static PyGetSetDef PyWcsprm_getset[] = {
  {"cd", (getter)PyWcsprm_get_cd, NULL, (char *)doc_cd},
  {"cdelt", (getter)PyWcsprm_get_cdelt, NULL, (char *)doc_cdelt},
  {"crota", (getter)PyWcsprm_get_crota, NULL, (char *)doc_crota},
  {"crpix", (getter)PyWcsprm_get_crpix, NULL, (char *)doc_crpix},
  {"crval", (getter)PyWcsprm_get_crval, NULL, (char *)doc_crval},
  {"ctype", (getter)PyWcsprm_get_ctype, NULL, (char *)doc_ctype},
  {"cubeface", (getter)PyWcsprm_get_cubeface, NULL, (char *)doc_cubeface},
  {"cunit", (getter)PyWcsprm_get_cunit, NULL, (char *)doc_cunit},
  {"lat", (getter)PyWcsprm_get_lat, NULL, (char *)doc_lat},
  {"latpole", (getter)PyWcsprm_get_latpole, NULL, (char *)doc_latpole},
  {"lng", (getter)PyWcsprm_get_lng, NULL, (char *)doc_lng},
  {"lonpole", (getter)PyWcsprm_get_lonpole, NULL, (char *)doc_lonpole},
  {"naxis", (getter)PyWcsprm_get_naxis, NULL, (char *)doc_naxis},
  {"pc", (getter)PyWcsprm_get_pc, NULL, (char *)doc_pc},
  {"ps", (getter)PyWcsprm_get_ps, NULL, (char *)doc_ps},
  {"pv", (getter)PyWcsprm_get_pv, NULL, (char *)doc_pv},
  {"restfrq", (getter)PyWcsprm_get_restfrq, NULL, (char *)doc_restfrq},
  {"restwav", (getter)PyWcsprm_get_restwav, NULL, (char *)doc_restwav},
  {"spec", (getter)PyWcsprm_get_spec, NULL, (char *)doc_spec},
  {NULL}
};

static PyMethodDef PyWcsprm_methods[] = {
  {"celfix", (PyCFunction)PyWcsprm_celfix, METH_NOARGS, doc_celfix},
  {"copy", (PyCFunction)PyWcsprm_copy, METH_NOARGS, doc_copy},
  {"__copy__", (PyCFunction)PyWcsprm_copy, METH_NOARGS, doc_copy},
  {"cylfix", (PyCFunction)PyWcsprm_cylfix, METH_VARARGS, doc_cylfix},
  {"datfix", (PyCFunction)PyWcsprm_datfix, METH_NOARGS, doc_datfix},
  {"fix", (PyCFunction)PyWcsprm_fix, METH_VARARGS, doc_fix},
  {"has_cdi_ja", (PyCFunction)PyWcsprm_has_cdi_ja, METH_NOARGS, doc_has_cdi_ja},
  {"has_crotaia", (PyCFunction)PyWcsprm_has_crotaia, METH_NOARGS, doc_has_crotaia},
  {"has_pci_ja", (PyCFunction)PyWcsprm_has_pci_ja, METH_NOARGS, doc_has_pci_ja},
  {"mix", (PyCFunction)PyWcsprm_mix, METH_VARARGS, doc_mix},
  {"p2s", (PyCFunction)PyWcsprm_p2s, METH_O, doc_p2s},
  {"print_contents", (PyCFunction)PyWcsprm_print_contents, METH_NOARGS, doc_print_contents},
  {"s2p", (PyCFunction)PyWcsprm_s2p, METH_O, doc_s2p},
  {"set", (PyCFunction)PyWcsprm_set, METH_NOARGS, doc_set},
  {"spcfix", (PyCFunction)PyWcsprm_spcfix, METH_NOARGS, doc_spcfix},
  {"sptr", (PyCFunction)PyWcsprm_sptr, METH_NOARGS, doc_sptr},
  {"unitfix", (PyCFunction)PyWcsprm_unitfix, METH_VARARGS, doc_unitfix},
  {NULL}
};

static PyTypeObject PyWcsprmType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "pywcs.Wcs",               /*tp_name*/
    sizeof(PyWcsprm),          /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PyWcsprm_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    doc_Wcs,                   /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    PyWcsprm_methods,          /* tp_methods */
    PyWcsprm_members,          /* tp_members */
    PyWcsprm_getset,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PyWcsprm_init,   /* tp_init */
    0,                         /* tp_alloc */
    PyWcsprm_new,              /* tp_new */
};

/***************************************************************************
 * Free functions                                                          *
 ***************************************************************************/

static PyObject*
pywcs_parse_image_header(PyObject* self, PyObject* args, PyObject* kwargs) {
  char*          header        = NULL;
  int            header_length = 0;
  int            nkeyrec       = 0;
  int            relax         = 0;
  int            ctrl          = 0;
  int            nreject       = 0;
  int            nwcs          = 0;
  struct wcsprm* wcs           = NULL;
  int            status        = 0;
  PyObject*      result        = NULL;
  struct wcsprm* wcsprm_obj    = NULL;
  PyObject*      pywcsprm_obj  = NULL;
  int            i             = 0;

  const char* keywords[] = {"header", "relax", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|i", (char **)keywords,
                                   &header, &header_length, &relax))
    return NULL;

  if (relax)
    relax = WCSHDR_all;

  nkeyrec = header_length / 80;

  // Make the call
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
    return NULL;
  }

  result = PyList_New(nwcs);
  if (result == NULL)
    return NULL;

  for (i = 0; i < nwcs; ++i) {
    // We have to copy all of the wcsprm objects here, since they
    // were all allocated as a single array (by wcspih), but they
    // need to show up as individual objects in Python, that can be
    // deallocated individually.  So this isn't optimal, but hopefully
    // these objects don't take long to copy.
    pywcsprm_obj = NULL;
    wcsprm_obj = malloc(sizeof(struct wcsprm));
    if (wcsprm_obj != NULL) {
      wcsprm_obj->flag = -1;
      if (wcscopy(1, wcs + i, wcsprm_obj) == 0) {
        pywcsprm_obj = (PyObject*)PyWcsprm_cnew(wcsprm_obj);
      } else {
        wcsfree(wcsprm_obj);
        free(wcsprm_obj);
      }
    }

    // Something above didn't work, so free everything
    if (pywcsprm_obj == NULL ||
        PyList_SetItem(result, (Py_ssize_t)i, pywcsprm_obj)) {
      Py_XDECREF(pywcsprm_obj);
      Py_DECREF(result);
      wcsvfree(&nwcs, &wcs);
      PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
      return NULL;
    }
  }

  wcsvfree(&nwcs, &wcs);
  return result;
}

static PyMethodDef module_methods[] = {
  {"parse_image_header", (PyCFunction)pywcs_parse_image_header, METH_KEYWORDS, doc_parse_image_header},
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

  m = Py_InitModule3("_pywcs", module_methods,
"The routines in this module implement the FITS World Coordinate System\n"
"(WCS) standard which defines methods to be used for computing world\n"
"coordinates from image pixel coordinates, and vice versa."
                     );

  if (m == NULL)
    return;

  Py_INCREF(&PyWcsprmType);
  PyModule_AddObject(m, "Wcs", (PyObject *)&PyWcsprmType);

  DEFINE_EXCEPTION(SingularMatrix);
  DEFINE_EXCEPTION(InconsistentAxisTypes);
  DEFINE_EXCEPTION(InvalidTransform);
  DEFINE_EXCEPTION(InvalidCoordinate);
  DEFINE_EXCEPTION(NoSolution);
  DEFINE_EXCEPTION(InvalidSubimageSpecification);
  DEFINE_EXCEPTION(NonseparableSubimageCoordinateSystem);

  PyModule_AddObject(m, "__docformat__", PyString_FromString("epytext"));
}
