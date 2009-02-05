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

#include "sip_wrap.h"
#include "docstrings.h"
#include "wcs.h"

static void
PySip_dealloc(
    PySip* self) {

  sip_free(&self->x);
  self->ob_type->tp_free((PyObject*)self);
}

/*@null@*/ static PyObject *
PySip_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PySip* self;

  self = (PySip*)type->tp_alloc(type, 0);
  if (self != NULL) {
    sip_clear(&self->x);
  }
  return (PyObject*)self;
}

static int
convert_matrix(
    /*@null@*/ PyObject* pyobj,
    PyArrayObject** array,
    double** data,
    unsigned int* order) {

  if (pyobj == Py_None) {
    *array = NULL;
    *data = NULL;
    *order = 0;
    return 0;
  }

  *array = (PyArrayObject*)PyArray_ContiguousFromAny(
      pyobj, PyArray_DOUBLE, 2, 2);
  if (*array == NULL) {
    return -1;
  }

  if (PyArray_DIM(*array, 0) != PyArray_DIM(*array, 1)) {
    PyErr_SetString(PyExc_ValueError,
                    "Matrix must be square.");
    return -1;
  }

  *data = (double*)PyArray_DATA(*array);
  *order = (unsigned int)PyArray_DIM(*array, 0) - 1;

  return 0;
}

static int
PySip_init(
    PySip* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds) {

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
  double*        a_data   = NULL;
  double*        b_data   = NULL;
  double*        ap_data  = NULL;
  double*        bp_data  = NULL;
  unsigned int   a_order  = 0;
  unsigned int   b_order  = 0;
  unsigned int   ap_order = 0;
  unsigned int   bp_order = 0;
  int            status   = -1;

  if (!PyArg_ParseTuple(args, "OOOOO:Sip.__init__",
                        &py_a, &py_b, &py_ap, &py_bp, &py_crpix)) {
    return -1;
  }

  if (convert_matrix(py_a, &a, &a_data, &a_order) ||
      convert_matrix(py_b, &b, &b_data, &b_order) ||
      convert_matrix(py_ap, &ap, &ap_data, &ap_order) ||
      convert_matrix(py_bp, &bp, &bp_data, &bp_order)) {
    goto exit;
  }

  crpix = (PyArrayObject*)PyArray_ContiguousFromAny(py_crpix, PyArray_DOUBLE,
                                                    1, 1);
  if (crpix == NULL) {
    goto exit;
  }

  if (PyArray_DIM(crpix, 0) != 2) {
    PyErr_SetString(PyExc_ValueError, "CRPIX wrong length");
    goto exit;
  }

  status = sip_init(&self->x,
                    a_order, a_data,
                    b_order, b_data,
                    ap_order, ap_data,
                    bp_order, bp_data,
                    PyArray_DATA(crpix));

 exit:
  Py_XDECREF(a);
  Py_XDECREF(b);
  Py_XDECREF(ap);
  Py_XDECREF(bp);
  Py_XDECREF(crpix);

  if (status == 0) {
    return 0;
  } else if (status > 0 && status < WCS_ERRMSG_MAX) {
    PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
    return -1;
  } else if (status == -1) {
    return -1;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return -1;
  }
}

/*@null@*/ static PyObject*
PySip_pix2foc(
    PySip* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* pixcrd     = NULL;
  PyArrayObject* foccrd     = NULL;
  int            status     = -1;
  const char*    keywords[] = {
    "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:pix2foc", (char **)keywords,
                                   &pixcrd_obj, &origin)) {
    return NULL;
  }

  if (self->x.a == NULL || self->x.b == NULL) {
    PyErr_SetString(
        PyExc_ValueError,
        "SIP object does not have coefficients for pix2foc transformation (A and B)");
    return NULL;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(pixcrd_obj, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    goto exit;
  }

  if (PyArray_DIM(pixcrd, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto exit;
  }

  foccrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd),
                                             PyArray_DOUBLE);
  if (foccrd == NULL) {
    status = 2;
    goto exit;
  }

  preoffset_array(pixcrd, origin);
  status = sip_pix2foc(&self->x,
                       (unsigned int)PyArray_DIM(pixcrd, 1),
                       (unsigned int)PyArray_DIM(pixcrd, 0),
                       (const double*)PyArray_DATA(pixcrd),
                       (double*)PyArray_DATA(foccrd));
  unoffset_array(pixcrd, origin);

 exit:

  Py_XDECREF(pixcrd);

  if (status == 0) {
    return (PyObject*)foccrd;
  } else {
    Py_XDECREF(foccrd);
    if (status == -1) {
      return NULL;
    } else if (status > 0 && status < WCS_ERRMSG_MAX) {
      PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
      return NULL;
    } else {
      PyErr_SetString(
          PyExc_RuntimeError,
          "Unknown error occurred.  Something is seriously wrong.");
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
PySip_foc2pix(
    PySip* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      foccrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* foccrd     = NULL;
  PyArrayObject* pixcrd     = NULL;
  int            status     = -1;
  const char*    keywords[] = {
    "foccrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:foc2pix", (char **)keywords,
                                   &foccrd_obj, &origin)) {
    return NULL;
  }

  if (self->x.ap == NULL || self->x.bp == NULL) {
    PyErr_SetString(
        PyExc_ValueError,
        "SIP object does not have coefficients for foc2pix transformation (AP and BP)");
    return NULL;
  }

  foccrd = (PyArrayObject*)PyArray_ContiguousFromAny(foccrd_obj, PyArray_DOUBLE, 2, 2);
  if (foccrd == NULL) {
    goto exit;
  }

  if (PyArray_DIM(foccrd, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto exit;
  }

  pixcrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(foccrd),
                                             PyArray_DOUBLE);
  if (pixcrd == NULL) {
    status = 2;
    goto exit;
  }

  preoffset_array(foccrd, origin);
  status = sip_foc2pix(&self->x,
                       (unsigned int)PyArray_DIM(pixcrd, 1),
                       (unsigned int)PyArray_DIM(pixcrd, 0),
                       (double*)PyArray_DATA(foccrd),
                       (double*)PyArray_DATA(pixcrd));
  unoffset_array(foccrd, origin);

 exit:
  Py_XDECREF(foccrd);

  if (status == 0) {
    return (PyObject*)pixcrd;
  } else {
    Py_XDECREF(pixcrd);
    if (status == -1) {
      return NULL;
    } else if (status > 0 && status < WCS_ERRMSG_MAX) {
      PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
      return NULL;
    } else {
      PyErr_SetString(
          PyExc_RuntimeError,
          "Unknown error occurred.  Something is seriously wrong.");
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
PySip_get_a(
    PySip* self,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {
    (npy_intp)self->x.a_order + 1,
    (npy_intp)self->x.a_order + 1 };

  if (is_null(self->x.a)) {
    return NULL;
  }

  return get_double_array("a", self->x.a, 2, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PySip_get_b(
    PySip* self,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {
    (npy_intp)self->x.b_order + 1,
    (npy_intp)self->x.b_order + 1 };

  if (is_null(self->x.b)) {
    return NULL;
  }

  return get_double_array("b", self->x.b, 2, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PySip_get_ap(
    PySip* self,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {
    (npy_intp)self->x.ap_order + 1,
    (npy_intp)self->x.ap_order + 1 };

  if (is_null(self->x.ap)) {
    return NULL;
  }

  return get_double_array("ap", self->x.ap, 2, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PySip_get_bp(
    PySip* self,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {
    (npy_intp)self->x.bp_order + 1,
    (npy_intp)self->x.bp_order + 1 };

  if (is_null(self->x.bp)) {
    return NULL;
  }

  return get_double_array("bp", self->x.bp, 2, dims, (PyObject*)self);
}

static PyObject*
PySip_get_a_order(
    PySip* self,
    /*@unused@*/ void* closure) {

  return get_int("a_order", (long int)self->x.a_order);
}

static PyObject*
PySip_get_b_order(
    PySip* self,
    /*@unused@*/ void* closure) {

  return get_int("b_order", (long int)self->x.b_order);
}

static PyObject*
PySip_get_ap_order(
    PySip* self,
    /*@unused@*/ void* closure) {

  return get_int("ap_order", (long int)self->x.ap_order);
}

static PyObject*
PySip_get_bp_order(
    PySip* self,
    /*@unused@*/ void* closure) {

  return get_int("bp_order", (long int)self->x.bp_order);
}

static PyObject*
PySip_get_crpix(
    PySip* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
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
  {"crpix", (getter)PySip_get_crpix, NULL, (char *)doc_crpix},
  {NULL}
};

static PyMethodDef PySip_methods[] = {
  {"pix2foc", (PyCFunction)PySip_pix2foc, METH_VARARGS, doc_sip_pix2foc},
  {"foc2pix", (PyCFunction)PySip_foc2pix, METH_VARARGS, doc_sip_foc2pix},
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

int
_setup_sip_type(
    PyObject* m) {

  if (PyType_Ready(&PySipType) < 0)
    return -1;

  Py_INCREF(&PySipType);
  return PyModule_AddObject(m, "Sip", (PyObject *)&PySipType);
}
