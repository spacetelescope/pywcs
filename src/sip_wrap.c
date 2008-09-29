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

#include "sip_wrap.h"
#include "docstrings.h"

static void
PySip_dealloc(PySip* self) {
  sip_free(&self->x);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
PySip_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
  PySip* self;

  self = (PySip*)type->tp_alloc(type, 0);
  if (self != NULL) {
    sip_clear(&self->x);
  }
  return (PyObject*)self;
}

static int
convert_matrix(PyObject* pyobj, PyArrayObject** array) {
  *array = (PyArrayObject*)PyArray_ContiguousFromAny(pyobj, PyArray_DOUBLE, 2, 2);
  if (*array == NULL)
    return -1;

  if (PyArray_DIM(*array, 0) != PyArray_DIM(*array, 1)) {
    PyErr_SetString(PyExc_ValueError,
                    "Matrix must be square.");
    return -1;
  }

  return 0;
}

static int
PySip_init(PySip* self, PyObject* args, PyObject* kwds) {
  PyObject*      py_a     = NULL;
  PyObject*      py_b     = NULL;
  PyObject*      py_ap    = NULL;
  PyObject*      py_bp    = NULL;
  PyObject*      py_crpix = NULL;
  PyArrayObject* a        = NULL;
  PyArrayObject* b        = NULL;
  PyArrayObject* ap       = NULL;
  PyArrayObject* bp       = NULL;
  PyArrayObject* crpix    = NULL;
  int            result   = 0;

  if (!PyArg_ParseTuple(args, "OOOOO:Sip.__init__",
                        &py_a, &py_b, &py_ap, &py_bp, &py_crpix)) {
    return -1;
  }

  if (convert_matrix(py_a, &a) ||
      convert_matrix(py_b, &b) ||
      convert_matrix(py_ap, &ap) ||
      convert_matrix(py_bp, &bp)) {
    result = -1;
    goto PySip_init_exit;
  }

  crpix = (PyArrayObject*)PyArray_ContiguousFromAny(py_crpix, PyArray_DOUBLE, 1, 1);
  if (crpix == NULL) {
    result = -1;
    goto PySip_init_exit;
  }

  if (PyArray_DIM(crpix, 0) != 2) {
    PyErr_SetString(PyExc_ValueError,
                    "CRPIX wrong length");
    result = -1;
    goto PySip_init_exit;
  }

  if (sip_init(&self->x,
               PyArray_DIM(a, 0) - 1, PyArray_DATA(a),
               PyArray_DIM(b, 0) - 1, PyArray_DATA(b),
               PyArray_DIM(ap, 0) - 1, PyArray_DATA(ap),
               PyArray_DIM(bp, 0) - 1, PyArray_DATA(bp),
               PyArray_DATA(crpix))) {
    result = -1;
    goto PySip_init_exit;
  }

 PySip_init_exit:
  Py_XDECREF(a);
  Py_XDECREF(b);
  Py_XDECREF(ap);
  Py_XDECREF(bp);
  Py_XDECREF(crpix);

  return result;
}

static PyObject*
PySip_pix2foc_generic(PySip* self, PyObject* arg, int do_shift) {
  PyArrayObject* pixcrd = NULL;
  PyArrayObject* foccrd = NULL;
  PyObject* result = NULL;

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(arg, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    goto _exit;
  }

  if (PyArray_DIM(pixcrd, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto _exit;
  }

  foccrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd),
                                             PyArray_DOUBLE);
  if (foccrd == NULL) {
    goto _exit;
  }

  if (do_shift)
    offset_array(pixcrd, 1.0);

  if (sip_pix2foc(&self->x,
                  PyArray_DIM(pixcrd, 0),
                  PyArray_DIM(pixcrd, 1),
                  (double*)PyArray_DATA(pixcrd),
                  (double*)PyArray_DATA(foccrd))) {
    PyErr_SetString(PyExc_RuntimeError, "Error correcting distortion");
    Py_XDECREF(foccrd);
    goto _exit;
  }

  if (do_shift)
    offset_array(pixcrd, -1.0);

  result = (PyObject*)foccrd;

 _exit:

  Py_XDECREF(pixcrd);

  return result;
}

static PyObject*
PySip_pix2foc(PySip* self, PyObject* arg, PyObject* kwds) {
  return PySip_pix2foc_generic(self, arg, 1);
}

static PyObject*
PySip_pix2foc_fits(PySip* self, PyObject* arg, PyObject* kwds) {
  return PySip_pix2foc_generic(self, arg, 1);
}

static PyObject*
PySip_foc2pix_generic(PySip* self, PyObject* arg, int do_shift) {
  PyArrayObject* foccrd = NULL;
  PyArrayObject* pixcrd = NULL;
  PyObject* result = NULL;

  foccrd = (PyArrayObject*)PyArray_ContiguousFromAny(arg, PyArray_DOUBLE, 2, 2);
  if (foccrd == NULL) {
    goto _exit;
  }

  if (PyArray_DIM(foccrd, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto _exit;
  }

  pixcrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(foccrd),
                                             PyArray_DOUBLE);
  if (pixcrd == NULL) {
    goto _exit;
  }

  if (do_shift)
    offset_array(pixcrd, 1.0);

  if (sip_foc2pix(&self->x,
                  PyArray_DIM(pixcrd, 0),
                  PyArray_DIM(pixcrd, 1),
                  (double*)PyArray_DATA(foccrd),
                  (double*)PyArray_DATA(pixcrd))) {
    PyErr_SetString(PyExc_RuntimeError, "Error correcting distortion");
    Py_XDECREF(foccrd);
    goto _exit;
  }

  if (do_shift)
    offset_array(pixcrd, -1.0);

  result = (PyObject*)pixcrd;

 _exit:

  Py_XDECREF(foccrd);

  return result;
}

static PyObject*
PySip_foc2pix(PySip* self, PyObject* arg, PyObject* kwds) {
  return PySip_foc2pix_generic(self, arg, 1);
}

static PyObject*
PySip_foc2pix_fits(PySip* self, PyObject* arg, PyObject* kwds) {
  return PySip_foc2pix_generic(self, arg, 1);
}

static PyObject*
PySip_get_a(PySip* self, void* closure) {
  const npy_intp dims[2] = {
    self->x.a_order + 1,
    self->x.a_order + 1 };

  if (is_null(self->x.a)) {
    return NULL;
  }

  return get_double_array("a", self->x.a, 2, dims, (PyObject*)self);
}

static PyObject*
PySip_get_b(PySip* self, void* closure) {
  const npy_intp dims[2] = {
    self->x.b_order + 1,
    self->x.b_order + 1 };

  if (is_null(self->x.b)) {
    return NULL;
  }

  return get_double_array("b", self->x.b, 2, dims, (PyObject*)self);
}

static PyObject*
PySip_get_ap(PySip* self, void* closure) {
  const npy_intp dims[2] = {
    self->x.ap_order + 1,
    self->x.ap_order + 1 };

  if (is_null(self->x.ap)) {
    return NULL;
  }

  return get_double_array("ap", self->x.ap, 2, dims, (PyObject*)self);
}

static PyObject*
PySip_get_bp(PySip* self, void* closure) {
  const npy_intp dims[2] = {
    self->x.bp_order + 1,
    self->x.bp_order + 1 };

  if (is_null(self->x.bp)) {
    return NULL;
  }

  return get_double_array("bp", self->x.bp, 2, dims, (PyObject*)self);
}

static PyObject*
PySip_get_a_order(PySip* self, void* closure) {
  return get_int("a_order", self->x.a_order);
}

static PyObject*
PySip_get_b_order(PySip* self, void* closure) {
  return get_int("b_order", self->x.b_order);
}

static PyObject*
PySip_get_ap_order(PySip* self, void* closure) {
  return get_int("ap_order", self->x.ap_order);
}

static PyObject*
PySip_get_bp_order(PySip* self, void* closure) {
  return get_int("bp_order", self->x.bp_order);
}

static PyGetSetDef PySip_getset[] = {
  {"a", (getter)PySip_get_a, NULL, (char *)doc_a},
  {"a_order", (getter)PySip_get_a_order, NULL, (char *)doc_a_order},
  {"b", (getter)PySip_get_b, NULL, (char *)doc_b},
  {"b_order", (getter)PySip_get_b_order, NULL, (char *)doc_b_order},
  {"ap", (getter)PySip_get_ap, NULL, (char *)doc_ap},
  {"ap_order", (getter)PySip_get_ap_order, NULL, (char *)doc_ap_order},
  {"bp", (getter)PySip_get_bp, NULL, (char *)doc_bp},
  {"bp_order", (getter)PySip_get_bp_order, NULL, (char *)doc_bp_order},
  {NULL}
};

static PyMethodDef PySip_methods[] = {
  {"pix2foc", (PyCFunction)PySip_pix2foc, METH_VARARGS, doc_sip_pix2foc},
  {"pix2foc_fits", (PyCFunction)PySip_pix2foc_fits, METH_VARARGS, doc_sip_pix2foc_fits},
  {"foc2pix", (PyCFunction)PySip_foc2pix, METH_VARARGS, doc_sip_foc2pix},
  {"foc2pix_fits", (PyCFunction)PySip_foc2pix_fits, METH_VARARGS, doc_sip_foc2pix_fits},
  {NULL}
};

PyTypeObject PySipType = {
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  "pywcs.Sip",                  /*tp_name*/
  sizeof(PySip),                /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PySip_dealloc,    /*tp_dealloc*/
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
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  doc_Sip,                      /* tp_doc */
  0,                            /* tp_traverse */
  0,                            /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PySip_methods,                /* tp_methods */
  0,                            /* tp_members */
  PySip_getset,                 /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc)PySip_init,         /* tp_init */
  0,                            /* tp_alloc */
  PySip_new,                    /* tp_new */
};

int _setup_sip_type(PyObject* m) {
  if (PyType_Ready(&PySipType) < 0)
    return -1;

  Py_INCREF(&PySipType);
  return PyModule_AddObject(m, "Sip", (PyObject *)&PySipType);
}
