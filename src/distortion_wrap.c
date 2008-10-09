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

#include "distortion_wrap.h"
#include "docstrings.h"

#include <structmember.h> /* From Python */

static int
PyDistLookup_traverse(
    PyDistLookup* self,
    visitproc visit,
    void* arg) {

  Py_VISIT(self->py_data);

  return 0;
}

static int
PyDistLookup_clear(
    PyDistLookup* self) {

  PyObject* tmp;

  tmp = (PyObject*)self->py_data;
  self->py_data = NULL;
  Py_XDECREF(tmp);

  return 0;
}

static void
PyDistLookup_dealloc(
    PyDistLookup* self) {

  distortion_lookup_t_free(&self->x);
  Py_XDECREF(self->py_data);
  self->ob_type->tp_free((PyObject*)self);
}

/*@null@*/ static PyObject *
PyDistLookup_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyDistLookup* self;

  self = (PyDistLookup*)type->tp_alloc(type, 0);
  if (self != NULL) {
    if (distortion_lookup_t_init(&self->x)) {
      return NULL;
    }
    self->py_data = NULL;
  }
  return (PyObject*)self;
}

static int
PyDistLookup_init(
    PyDistLookup* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyObject* py_array_obj = NULL;
  PyArrayObject* array_obj = NULL;

  if (!PyArg_ParseTuple(args, "O(dd)(dd)(dd):DistortionLookupTable.__init__",
                       &py_array_obj,
                       &(self->x.crpix[0]), &(self->x.crpix[1]),
                       &(self->x.crval[0]), &(self->x.crval[1]),
                        &(self->x.cdelt[0]), &(self->x.cdelt[1]))) {
    return -1;
  }

  array_obj = (PyArrayObject*)PyArray_ContiguousFromAny(py_array_obj, PyArray_FLOAT32, 2, 2);
  if (array_obj == NULL) {
    return -1;
  }

  self->py_data = array_obj;
  self->x.naxis[0] = (unsigned int)PyArray_DIM(array_obj, 0);
  self->x.naxis[1] = (unsigned int)PyArray_DIM(array_obj, 1);
  self->x.data = (float *)PyArray_DATA(array_obj);

  return 0;
}

static PyObject*
PyDistLookup_get_cdelt(
    PyDistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_cdelt(
    PyDistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return set_double_array("cdelt", value, 1, &naxis, self->x.cdelt);
}

static PyObject*
PyDistLookup_get_crpix(
    PyDistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_crpix(
    PyDistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return set_double_array("crpix", value, 1, &naxis, self->x.crpix);
}

static PyObject*
PyDistLookup_get_crval(
    PyDistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("crval", self->x.crval, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_crval(
    PyDistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return set_double_array("crval", value, 1, &naxis, self->x.crval);
}

/*@shared@*/ static PyObject*
PyDistLookup_get_data(
    PyDistLookup* self,
    /*@unused@*/ void* closure) {

  if (self->py_data == NULL) {
    Py_INCREF(Py_None);
    return Py_None;
  } else {
    Py_INCREF(self->py_data);
    return (PyObject*)self->py_data;
  }
}

static int
PyDistLookup_set_data(
    PyDistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  PyArrayObject* value_array = NULL;

  if (value == NULL) {
    Py_XDECREF(self->py_data);
    self->py_data = NULL;
    self->x.data = NULL;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_FLOAT32, 2, 2);

  if (value_array == NULL) {
    return -1;
  }

  Py_XDECREF(self->py_data);

  self->py_data = value_array;
  self->x.naxis[0] = (unsigned int)PyArray_DIM(value_array, 0);
  self->x.naxis[1] = (unsigned int)PyArray_DIM(value_array, 1);
  self->x.data = (float *)PyArray_DATA(value_array);

  return 0;
}

/*@null@*/ static PyObject*
PyDistLookup_get_offset(
    PyDistLookup* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  double coord[NAXES];
  double result;

  if (self->x.data == NULL) {
    PyErr_SetString(PyExc_RuntimeError,
                    "No data has been set for the lookup table");
    return NULL;
  }

  if (!PyArg_ParseTuple(args, "dd:get_offset", &coord[0], &coord[1])) {
    return NULL;
  }

  result = get_distortion_offset(&self->x, coord);
  return PyFloat_FromDouble(result);
}

static PyGetSetDef PyDistLookup_getset[] = {
  {"cdelt", (getter)PyDistLookup_get_cdelt, (setter)PyDistLookup_set_cdelt, (char *)doc_cdelt},
  {"crpix", (getter)PyDistLookup_get_crpix, (setter)PyDistLookup_set_crpix, (char *)doc_crpix},
  {"crval", (getter)PyDistLookup_get_crval, (setter)PyDistLookup_set_crval, (char *)doc_crval},
  {"data",  (getter)PyDistLookup_get_data,  (setter)PyDistLookup_set_data,  (char *)doc_data},
  {NULL}
};

static PyMethodDef PyDistLookup_methods[] = {
  {"get_offset", (PyCFunction)PyDistLookup_get_offset, METH_VARARGS, doc_get_offset},
  {NULL}
};

PyTypeObject PyDistLookupType = {
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  "pywcs.DistortionLookupTable", /*tp_name*/
  sizeof(PyDistLookup),         /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PyDistLookup_dealloc, /*tp_dealloc*/
  0,                            /*tp_print*/
  0,                            /*tp_getattr*/
  0,                            /*tp_setattr*/
  0,                            /*tp_compare*/
  0,                            /*tp_repr*/
  0,                            /*tp_as_number*/
  0,                            /*tp_as_sequence*/
  0,                            /*tp_as_mapping*/
  0,                            /*tp_hash */
  0,                            /*tp_call*/
  0,                            /*tp_str*/
  0,                            /*tp_getattro*/
  0,                            /*tp_setattro*/
  0,                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
  doc_DistortionLookupTable,    /* tp_doc */
  (traverseproc)PyDistLookup_traverse, /* tp_traverse */
  (inquiry)PyDistLookup_clear,  /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PyDistLookup_methods,         /* tp_methods */
  0,                            /* tp_members */
  PyDistLookup_getset,          /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc)PyDistLookup_init,  /* tp_init */
  0,                            /* tp_alloc */
  PyDistLookup_new,             /* tp_new */
};

int _setup_distortion_type(
    PyObject* m) {

  if (PyType_Ready(&PyDistLookupType) < 0) {
    return -1;
  }

  Py_INCREF(&PyDistLookupType);
  return PyModule_AddObject(m, "DistortionLookupTable", (PyObject *)&PyDistLookupType);
}
