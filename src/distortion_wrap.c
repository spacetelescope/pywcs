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

#include "distortion.h"
#include "docstrings.h"
#include "util.h"

#include <structmember.h>

static PyTypeObject PyDistLookupType;

typedef struct {
  PyObject_HEAD
  struct distortion_lookup_t x;
  PyArrayObject* py_data;
} PyDistLookup;

static void
PyDistLookup_dealloc(PyDistLookup* self) {
  distortion_lookup_t_free(&self->x);
  Py_XDECREF(self->py_data);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
PyDistLookup_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
  PyDistLookup* self;

  self = (PyDistLookup*)type->tp_alloc(type, 0);
  if (self != NULL) {
    distortion_lookup_t_init(&self->x);
    self->py_data = NULL;
  }
  return (PyObject*)self;
}

static int
PyDistLookup_init(PyDistLookup* self, PyObject* args, PyObject* kwds) {
  PyObject* py_array_obj = NULL;
  PyArrayObject* array_obj = NULL;

  if (!PyArg_ParseTuple(args, "O(dd)(dd)(dd)", &py_array_obj,
                        &(self->x.crpix[0]), &(self->x.crpix[1]),
                        &(self->x.crval[0]), &(self->x.crval[1]),
                        &(self->x.cdelt[0]), &(self->x.cdelt[1]))) {
    return -1;
  }

  array_obj = (PyArrayObject*)PyArray_ContiguousFromAny(py_array_obj, 2, 2, PyArray_DOUBLE);
  if (array_obj == NULL)
    return -1;

  self->py_data = array_obj;
  self->x.naxis[0] = PyArray_DIM(array_obj, 0);
  self->x.naxis[1] = PyArray_DIM(array_obj, 1);
  self->x.data = (double *)PyArray_DATA(array_obj);

  return 0;
}

static PyObject*
PyDistLookup_get_cdelt(PyDistLookup* self, void* closure) {
  Py_ssize_t naxis = 2;

  return PyArrayProxy_New((PyObject*)self, 1, &naxis, PyArray_DOUBLE,
                          self->x.cdelt);
}

static int
PyDistLookup_set_cdelt(PyDistLookup* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the cdelt attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);

  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2) {
    PyErr_SetString(PyExc_ValueError, "len(cdelt) != 2");
    Py_DECREF(value_array);
    return -1;
  }

  copy_array_to_c_double(value_array, self->x.cdelt);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyDistLookup_get_crpix(PyDistLookup* self, void* closure) {
  Py_ssize_t naxis = 2;

  return PyArrayProxy_New((PyObject*)self, 1, &naxis, PyArray_DOUBLE,
                          self->x.crpix);
}

static int
PyDistLookup_set_crpix(PyDistLookup* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the crpix attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2) {
    PyErr_SetString(PyExc_ValueError, "len(crpix) != 2");
    return -1;
  }

  copy_array_to_c_double(value_array, self->x.crpix);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyDistLookup_get_crval(PyDistLookup* self, void* closure) {
  Py_ssize_t naxis = 2;

  return PyArrayProxy_New((PyObject*)self, 1, &naxis, PyArray_DOUBLE,
                          self->x.crval);
}

static int
PyDistLookup_set_crval(PyDistLookup* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the crval attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2) {
    PyErr_SetString(PyExc_ValueError, "len(crval) != 2");
    Py_DECREF(value_array);
    return -1;
  }

  copy_array_to_c_double(value_array, self->x.crval);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyDistLookup_get_data(PyDistLookup* self, void* closure) {
  return (PyObject*)self->py_data;
}

static int
PyDistLookup_set_data(PyDistLookup* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the data attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, 2, 2, PyArray_DOUBLE);

  if (value_array == NULL)
    return -1;

  Py_XDECREF(self->py_data);

  self->py_data = value_array;
  self->x.naxis[0] = PyArray_DIM(value_array, 0);
  self->x.naxis[1] = PyArray_DIM(value_array, 1);
  self->x.data = (double *)PyArray_DATA(value_array);

  return 0;
}

static PyObject*
PyDistLookup_get_offset(PyDistLookup* self, PyObject* args, PyObject* kwds) {
  double coord[NAXES];
  double result;

  if (!PyArg_ParseTuple(args, "dd", &coord[0], &coord[1])) {
    return NULL;
  }

  result = get_distortion_offset(&self->x, coord);
  return PyFloat_FromDouble(result);
}

static PyMemberDef PyDistLookup_members[] = {
  {NULL}
};

static PyGetSetDef PyDistLookup_getset[] = {
  {"cdelt", (getter)PyDistLookup_get_cdelt, (setter)PyDistLookup_set_cdelt, (char *)doc_cdelt},
  {"crpix", (getter)PyDistLookup_get_crpix, (setter)PyDistLookup_set_crpix, (char *)doc_crpix},
  {"crval", (getter)PyDistLookup_get_crval, (setter)PyDistLookup_set_crval, (char *)doc_crval},
  {"data", (getter)PyDistLookup_get_data, (setter)PyDistLookup_set_data, (char *)doc_data},
  {NULL}
};

static PyMethodDef PyDistLookup_methods[] = {
  {"get_offset", (PyCFunction)PyDistLookup_get_offset, METH_VARARGS, doc_get_offset},
  {NULL}
};

static PyTypeObject PyDistLookupType = {
  PyObject_HEAD_INIT(NULL)
  0,                          /*ob_size*/
  "pywcs.DistortionLookupTable",  /*tp_name*/
  sizeof(PyDistLookup),           /*tp_basicsize*/
  0,                          /*tp_itemsize*/
  (destructor)PyDistLookup_dealloc, /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  0,                          /*tp_repr*/
  0,                          /*tp_as_number*/
  0,                          /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash */
  0,                          /*tp_call*/
  0,                          /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  doc_DistortionLookupTable,  /* tp_doc */
  0,                          /* tp_traverse */
  0,                          /* tp_clear */
  0,                          /* tp_richcompare */
  0,                          /* tp_weaklistoffset */
  0,                          /* tp_iter */
  0,                          /* tp_iternext */
  PyDistLookup_methods,           /* tp_methods */
  PyDistLookup_members,           /* tp_members */
  PyDistLookup_getset,            /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tp_dictoffset */
  (initproc)PyDistLookup_init,    /* tp_init */
  0,                          /* tp_alloc */
  PyDistLookup_new,               /* tp_new */
};

static PyTypeObject PyDistortionType;

typedef struct {
  PyObject_HEAD
  struct distortion_t x;
  /* TODO: The Python wrapper should store references to the distortion
     Python arrays so that memory isn't prematurely freed. */
  PyObject* py_pre_dist[NAXES];
  PyObject* py_post_dist[NAXES];
} PyDistortion;

static void
PyDistortion_dealloc(PyDistortion* self) {
  unsigned int i;

  distortion_t_free(&self->x);
  for (i = 0; i < NAXES; ++i) {
    Py_XDECREF(self->py_pre_dist[i]);
    Py_XDECREF(self->py_post_dist[i]);
  }
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
PyDistortion_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
  PyDistortion* self;
  unsigned int i;

  self = (PyDistortion*)type->tp_alloc(type, 0);
  if (self != NULL) {
    distortion_t_init(&self->x);
    for (i = 0; i < NAXES; ++i) {
      self->py_pre_dist[i] = NULL;
      self->py_post_dist[i] = NULL;
    }
  }
  return (PyObject*)self;
}

static int
PyDistortion_init(PyDistortion* self, PyObject* args, PyObject* kwds) {
  /* TODO: Write me */

  return 0;
}

static PyObject*
PyDistortion_get_cd(PyDistortion* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x.has_pc == 1) {
    PyErr_SetString(PyExc_AttributeError, "No cd is present.");
    return NULL;
  }

  return PyArrayProxy_New((PyObject*)self, 2, dims, PyArray_DOUBLE, self->x.pc);
}

static int
PyDistortion_set_cd(PyDistortion* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) { /* deletion */
    self->x.has_pc = 1;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          2, 2);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2 || PyArray_DIM(value_array, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "cd must be a 2x2 array");
    Py_DECREF(value_array);
    return -1;
  }

  copy_array_to_c_double(value_array, &(self->x.pc[0][0]));
  self->x.has_pc = 0;

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyDistortion_get_cdelt(PyDistortion* self, void* closure) {
  Py_ssize_t naxis = 2;

  if (self->x.has_pc == 0) {
    PyErr_SetString(PyExc_AttributeError, "No cdelt/pc is present.");
    return NULL;
  }

  return PyArrayProxy_New((PyObject*)self, 1, &naxis, PyArray_DOUBLE,
                          self->x.cdelt);
}

static int
PyDistortion_set_cdelt(PyDistortion* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) { /* deletion */
    self->x.has_pc = 0;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);

  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2) {
    PyErr_SetString(PyExc_ValueError, "len(cdelt) != 2");
    Py_DECREF(value_array);
    return -1;
  }

  copy_array_to_c_double(value_array, self->x.cdelt);
  self->x.has_pc = 1;

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyDistortion_get_cpdis(PyDistortion* self, void* closure) {
  return Py_BuildValue("OO", self->py_pre_dist[0], self->py_pre_dist[1]);
}

static int
PyDistortion_set_cpdis(PyDistortion* self, PyObject* value, void* closure) {
  Py_ssize_t i;
  PyObject* subvalue = NULL;

  if (!PySequence_Check(value) || PySequence_Size(value) != 2) {
    PyErr_SetString(PyExc_TypeError, "cpdis must be a 2-length sequence");
    return -1;
  }

  for (i = 0; i < NAXES; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (!PyObject_TypeCheck(subvalue, &PyDistLookupType)) {
      Py_XDECREF(subvalue);
      PyErr_SetString(PyExc_TypeError, "cpdis must be a 2-length sequence of DistortionLookupTable instances.");
      return -1;
    }

    Py_XDECREF(self->py_pre_dist[i]);
    self->py_pre_dist[i] = subvalue;
    self->x.pre_dist[i] = &(((PyDistLookup *)subvalue)->x);
  }

  return 0;
}

static PyObject*
PyDistortion_get_cqdis(PyDistortion* self, void* closure) {
  return Py_BuildValue("OO", self->py_post_dist[0], self->py_post_dist[1]);
}

static int
PyDistortion_set_cqdis(PyDistortion* self, PyObject* value, void* closure) {
  Py_ssize_t i;
  PyObject* subvalue = NULL;

  if (!PySequence_Check(value) || PySequence_Size(value) != 2) {
    PyErr_SetString(PyExc_TypeError, "cqdis must be a 2-length sequence");
    return -1;
  }

  for (i = 0; i < NAXES; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (!PyObject_TypeCheck(subvalue, &PyDistLookupType)) {
      Py_XDECREF(subvalue);
      PyErr_SetString(PyExc_TypeError, "cqdis must be a 2-length sequence of DistortionLookupTable instances.");
      return -1;
    }

    Py_XDECREF(self->py_post_dist[i]);
    self->py_post_dist[i] = subvalue;
    self->x.post_dist[i] = &(((PyDistLookup *)subvalue)->x);
  }

  return 0;
}

static PyObject*
PyDistortion_get_crpix(PyDistortion* self, void* closure) {
  Py_ssize_t naxis = 2;

  return PyArrayProxy_New((PyObject*)self, 1, &naxis, PyArray_DOUBLE,
                          self->x.crpix);
}

static int
PyDistortion_set_crpix(PyDistortion* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the crpix attribute");
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          1, 1);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2) {
    PyErr_SetString(PyExc_ValueError, "len(crpix) != 2");
    Py_DECREF(value_array);
    return -1;
  }

  copy_array_to_c_double(value_array, self->x.crpix);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyDistortion_get_crval(PyDistortion* self, void* closure) {
  Py_ssize_t naxis = 2;

  if (self->x.wcs.crval == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return NULL;
  }

  return PyArrayProxy_New((PyObject*)self, 1, &naxis, PyArray_DOUBLE,
                          self->x.wcs.crval);
}

static int
PyDistortion_set_crval(PyDistortion* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (self->x.wcs.crval == NULL) {
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

  if (PyArray_DIM(value_array, 0) != 2) {
    PyErr_SetString(PyExc_ValueError, "len(crval) != 2");
    Py_DECREF(value_array);
    return -1;
  }

  copy_array_to_c_double(value_array, self->x.wcs.crval);

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyDistortion_get_ctype(PyDistortion* self, void* closure) {
  return PyWcsprmListProxy_New((PyObject*)self, 2, self->x.wcs.ctype);
}

static int
PyDistortion_set_ctype(PyDistortion* self, PyObject* value, void* closure) {
  PyObject*  str      = NULL;
  char*      str_char = NULL;
  Py_ssize_t str_len  = 0;
  int        i        = 0;

  if (self->x.wcs.ctype == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return -1;
  }

  if (!PySequence_Check(value)) {
    PyErr_SetString(PyExc_TypeError, "ctype must be a sequence of strings");
    return -1;
  }

  if (PySequence_Size(value) != 2) {
    PyErr_SetString(PyExc_ValueError, "len(ctype) != 2");
    return -1;
  }

  for (i = 0; i < 2; ++i) {
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

    strncpy(self->x.wcs.ctype[i], str_char, 72);

    Py_DECREF(str);
  }

  return 0;
}

static PyObject*
PyDistortion_get_pc(PyDistortion* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x.has_pc == 0) {
    PyErr_SetString(PyExc_AttributeError, "No pc is present.");
    return NULL;
  }

  return PyArrayProxy_New((PyObject*)self, 2, dims, PyArray_DOUBLE, self->x.pc);
}

static int
PyDistortion_set_pc(PyDistortion* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) { /* deletion */
    self->x.has_pc = 0;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value,
                                                          PyArray_DOUBLE,
                                                          2, 2);
  if (value_array == NULL)
    return -1;

  if (PyArray_DIM(value_array, 0) != 2 || PyArray_DIM(value_array, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "pc must be a 2x2 array");
    Py_DECREF(value_array);
    return -1;
  }

  copy_array_to_c_double(value_array, &(self->x.pc[0][0]));
  self->x.has_pc = 1;

  Py_DECREF(value_array);

  return 0;
}

static PyObject*
PyDistortion_get_pv(PyDistortion* self, PyObject* args, PyObject* kwds) {
  PyObject* result    = NULL;
  PyObject* subresult = NULL;
  int       i         = 0;

  if (self->x.wcs.pv == NULL) {
    PyErr_SetString(PyExc_AssertionError, "No PVi_ma records present.");
    return NULL;
  }

  result = PyList_New(self->x.wcs.npv);
  if (result == NULL)
    return NULL;

  for (i = 0; i < self->x.wcs.npv; ++i) {
    subresult = Py_BuildValue("iid",
                              self->x.wcs.pv[i].i,
                              self->x.wcs.pv[i].m,
                              self->x.wcs.pv[i].value);
    if (subresult == NULL || PyList_SetItem(result, i, subresult)) {
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

static PyObject*
PyDistortion_set_pv(PyDistortion* self, PyObject* arg, PyObject* kwds) {
  PyObject*  subvalue  = NULL;
  int        i         = 0;
  Py_ssize_t size      = 0;
  int        ival      = 0;
  int        mval      = 0;
  double     value     = 0;

  if (self->x.wcs.ps == NULL) {
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

  if (size > self->x.wcs.npvmax) {
    free(self->x.wcs.pv);
    self->x.wcs.pv = malloc(sizeof(struct pscard) * size);
    if (self->x.wcs.pv == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
      return NULL;
    }
    self->x.wcs.npvmax = size;
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

    self->x.wcs.pv[i].i = ival;
    self->x.wcs.pv[i].m = mval;
    self->x.wcs.pv[i].value = value;
    self->x.wcs.npv = i + 1;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyMemberDef PyDistortion_members[] = {
  {NULL}
};

static PyGetSetDef PyDistortion_getset[] = {
  {"cd",    (getter)PyDistortion_get_cd,    (setter)PyDistortion_set_cd,    (char *)doc_cd},
  {"cdelt", (getter)PyDistortion_get_cdelt, (setter)PyDistortion_set_cdelt, (char *)doc_cdelt},
  {"cpdis", (getter)PyDistortion_get_cpdis, (setter)PyDistortion_set_cpdis, (char *)doc_cpdis},
  {"cqdis", (getter)PyDistortion_get_cqdis, (setter)PyDistortion_set_cqdis, (char *)doc_cqdis},
  {"crpix", (getter)PyDistortion_get_crpix, (setter)PyDistortion_set_crpix, (char *)doc_crpix},
  {"crval", (getter)PyDistortion_get_crval, (setter)PyDistortion_set_crval, (char *)doc_crval},
  {"ctype", (getter)PyDistortion_get_ctype, (setter)PyDistortion_set_ctype, (char *)doc_ctype},
  {"pc",    (getter)PyDistortion_get_pc,    (setter)PyDistortion_set_pc,    (char *)doc_pc},
  {NULL}
};

static PyMethodDef PyDistortion_methods[] = {
  {"get_pv", (PyCFunction)PyDistortion_get_pv, METH_NOARGS, doc_get_pv},
  {"set_pv", (PyCFunction)PyDistortion_set_pv, METH_O, doc_set_pv},
  {NULL}
};

static PyTypeObject PyDistortionType = {
  PyObject_HEAD_INIT(NULL)
  0,                          /*ob_size*/
  "pywcs._WCS",               /*tp_name*/
  sizeof(PyDistortion),           /*tp_basicsize*/
  0,                          /*tp_itemsize*/
  (destructor)PyDistortion_dealloc, /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  0,                          /*tp_repr*/
  0,                          /*tp_as_number*/
  0,                          /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash */
  0,                          /*tp_call*/
  0,                          /*tp_str*/
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
  PyDistortion_methods,           /* tp_methods */
  PyDistortion_members,           /* tp_members */
  PyDistortion_getset,            /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tp_dictoffset */
  (initproc)PyDistortion_init,    /* tp_init */
  0,                          /* tp_alloc */
  PyDistortion_new,               /* tp_new */
};

void _setup_distortion_type(PyObject* m) {
  if (PyType_Ready(&PyDistortionType) < 0)
    return;

  if (PyType_Ready(&PyDistLookupType) < 0)
    return;

  Py_INCREF(&PyDistortionType);
  PyModule_AddObject(m, "Distortion", (PyObject *)&PyDistortionType);

  Py_INCREF(&PyDistLookupType);
  PyModule_AddObject(m, "DistortionLookupTable", (PyObject *)&PyDistLookupType);
}
