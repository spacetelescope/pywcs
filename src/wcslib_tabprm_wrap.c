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

#include "wcslib_tabprm_wrap.h"

#include <wcs.h>
#include <tab.h>

/*
 It gets to be really tedious to type long docstrings in ANSI C syntax
 (since multi-line strings literals are not valid).  Therefore, the
 docstrings are written in doc/docstrings.py, which are then converted
 by setup.py into docstrings.h, which we include here.
*/
#include "docstrings.h"

/***************************************************************************
 * Helper functions                                                        *
 ***************************************************************************/

void
note_change(PyTabprm* self) {
  self->x->flag = 0;
}

int
make_fancy_dims(PyTabprm* self, npy_intp* ndims, npy_intp* dims) {
  npy_intp i, M;

  M = (npy_intp)self->x->M;
  if (M + 1 > NPY_MAXDIMS) {
    PyErr_SetString(PyExc_ValueError, "Too many dimensions");
    return -1;
  }

  *ndims = M + 1;

  for (i = 0; i < M; ++i) {
    dims[i] = self->x->K[M-1-i];
  }

  dims[M] = M;

  return 0;
}

PyObject** tab_errexc[6];

/***************************************************************************
 * PyTabprm methods
 */

static int
PyTabprm_traverse(
    PyTabprm* self, visitproc visit, void *arg) {
  int vret;

  vret = visit(self->owner, arg);
  if (vret != 0) {
    return vret;
  }

  return 0;
}

static int
PyTabprm_clear(
    PyTabprm* self) {
  PyObject* tmp;

  tmp = self->owner;
  self->owner = NULL;
  Py_XDECREF(tmp);

  return 0;
}

static void
PyTabprm_dealloc(
    PyTabprm* self) {

  PyTabprm_clear(self);
  self->ob_type->tp_free((PyObject*)self);
}

PyTabprm*
PyTabprm_cnew(PyObject* wcsprm, struct tabprm* x) {
  PyTabprm* self;
  self = (PyTabprm*)(&PyTabprmType)->tp_alloc(&PyTabprmType, 0);
  self->x = x;
  Py_INCREF(wcsprm);
  self->owner = wcsprm;
  return self;
}

static int
PyTabprm_cset(
    PyTabprm* self) {

  int status = 0;

  status = tabset(self->x);

  if (status == 0) {
    return 0;
  } else if (status > 0 && status < WCS_ERRMSG_MAX) {
    PyErr_SetString(*tab_errexc[status], tab_errmsg[status]);
    return -1;
  } else {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    return -1;
  }
}

/*@null@*/ static PyObject*
PyTabprm_set(
    PyTabprm* self) {

  if (PyTabprm_cset(self)) {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
PyTabprm_print_contents(
    PyTabprm* self) {

  int ignored;

  if (PyTabprm_cset(self)) {
    return NULL;
  }

  ignored = tabprt(self->x);

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
PyTabprm___str__(
    PyTabprm* self) {

  const int MAX_LEN = 1 << 18;
  char buffer[MAX_LEN+1];
  int out_pipe[2];
  int saved_stdout;
  int ignored;
  int len;

  buffer[0] = 0;

  if (PyTabprm_cset(self)) {
    return NULL;
  }

  saved_stdout = dup(STDOUT_FILENO);

  if (pipe(out_pipe) != 0) {
    PyErr_SetString(PyExc_RuntimeError, "Error creating pipe");
    return NULL;
  }

  dup2(out_pipe[1], STDOUT_FILENO);
  close(out_pipe[1]);

  ignored = tabprt(self->x);
  fflush(stdout);

  len = read(out_pipe[0], buffer, MAX_LEN);

  if (len >= MAX_LEN) {
    PyErr_SetString(PyExc_MemoryError, "Buffer overflow");
    return NULL;
  }

  dup2(saved_stdout, STDOUT_FILENO);

  return PyString_FromString(buffer);
}

/***************************************************************************
 * Member getters/setters (properties)
 */

/*@null@*/ static PyObject*
PyTabprm_get_coord(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  npy_intp ndims;
  npy_intp dims[NPY_MAXDIMS];

  if (is_null(self->x->coord)) {
    return NULL;
  }

  if (make_fancy_dims(self, &ndims, dims)) {
    return NULL;
  }

  return get_double_array("coord", self->x->coord, ndims, dims, (PyObject*)self);
}

/*@null@*/ static int
PyTabprm_set_coord(
    PyTabprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp ndims;
  npy_intp dims[NPY_MAXDIMS];

  if (is_null(self->x->coord)) {
    return -1;
  }

  if (make_fancy_dims(self, &ndims, dims)) {
    return -1;
  }

  return set_double_array("coord", value, ndims, dims, self->x->coord);
}

/*@null@*/ static PyObject*
PyTabprm_get_crval(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->crval)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_double_array("crval", self->x->crval, 1, &M, (PyObject*)self);
}

static int
PyTabprm_set_crval(
    PyTabprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp M = 0;

  if (is_null(self->x->crval)) {
    return -1;
  }

  M = (Py_ssize_t)self->x->M;

  note_change(self);

  return set_double_array("crval", value, 1, &M, self->x->crval);
}

/*@null@*/ static PyObject*
PyTabprm_get_delta(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->delta)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_double_array("delta", self->x->delta, 1, &M, (PyObject*)self);
}

/*@null@*/ static PyObject*
PyTabprm_get_extrema(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  npy_intp ndims;
  npy_intp dims[NPY_MAXDIMS];

  if (is_null(self->x->coord)) {
    return NULL;
  }

  if (make_fancy_dims(self, &ndims, dims)) {
    return NULL;
  }

  dims[ndims-2] = 2;

  return get_double_array("extrema", self->x->extrema, ndims, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PyTabprm_get_K(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->K)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_int_array("K", self->x->K, 1, &M, (PyObject*)self);
}

/*@null@*/ static PyObject*
PyTabprm_get_M(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  return get_int("M", self->x->M);
}

/*@null@*/ static PyObject*
PyTabprm_get_map(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->map)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_int_array("map", self->x->map, 1, &M, (PyObject*)self);
}

static int
PyTabprm_set_map(
    PyTabprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp M = 0;

  if (is_null(self->x->map)) {
    return -1;
  }

  M = (Py_ssize_t)self->x->M;

  note_change(self);

  return set_int_array("map", value, 1, &M, self->x->map);
}

/*@null@*/ static PyObject*
PyTabprm_get_nc(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  return get_int("nc", self->x->nc);
}

/*@null@*/ static PyObject*
PyTabprm_get_p0(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->p0)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_int_array("p0", self->x->p0, 1, &M, (PyObject*)self);
}

/*@null@*/ static PyObject*
PyTabprm_get_sense(
    PyTabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->sense)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_int_array("sense", self->x->sense, 1, &M, (PyObject*)self);
}

/***************************************************************************
 * PyTabprm definition structures
 */

static PyGetSetDef PyTabprm_getset[] = {
  {"coord", (getter)PyTabprm_get_coord, (setter)PyTabprm_set_coord, (char *)doc_coord},
  {"crval", (getter)PyTabprm_get_crval, (setter)PyTabprm_set_crval, (char *)doc_crval_tabprm},
  {"delta", (getter)PyTabprm_get_delta, NULL, (char *)doc_delta},
  {"extrema", (getter)PyTabprm_get_extrema, NULL, (char *)doc_extrema},
  {"K", (getter)PyTabprm_get_K, NULL, (char *)doc_K},
  {"M", (getter)PyTabprm_get_M, NULL, (char *)doc_M},
  {"map", (getter)PyTabprm_get_map, (setter)PyTabprm_set_map, (char *)doc_map},
  {"nc", (getter)PyTabprm_get_nc, NULL, (char *)doc_nc},
  {"p0", (getter)PyTabprm_get_p0, NULL, (char *)doc_p0},
  {"sense", (getter)PyTabprm_get_sense, NULL, (char *)doc_sense},
  {NULL}
};

static PyMethodDef PyTabprm_methods[] = {
  {"print_contents", (PyCFunction)PyTabprm_print_contents, METH_NOARGS, doc_print_contents_tabprm},
  {"set", (PyCFunction)PyTabprm_set, METH_NOARGS, doc_set_tabprm},
  {NULL}
};

PyTypeObject PyTabprmType = {
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  "pywcs.Tabprm",               /*tp_name*/
  sizeof(PyTabprm),             /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PyTabprm_dealloc, /*tp_dealloc*/
  0,                            /*tp_print*/
  0,                            /*tp_getattr*/
  0,                            /*tp_setattr*/
  0,                            /*tp_compare*/
  (reprfunc)PyTabprm___str__,   /*tp_repr*/
  0,                            /*tp_as_number*/
  0,                            /*tp_as_sequence*/
  0,                            /*tp_as_mapping*/
  0,                            /*tp_hash */
  0,                            /*tp_call*/
  (reprfunc)PyTabprm___str__,   /*tp_str*/
  0,                            /*tp_getattro*/
  0,                            /*tp_setattro*/
  0,                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  doc_Tabprm,                   /* tp_doc */
  (traverseproc)PyTabprm_traverse, /* tp_traverse */
  (inquiry)PyTabprm_clear,         /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PyTabprm_methods,             /* tp_methods */
  0,                            /* tp_members */
  PyTabprm_getset,              /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  0,                            /* tp_init */
  0,                            /* tp_alloc */
  0,                            /* tp_new */
};

int
_setup_tabprm_type(
    PyObject* m) {

  if (PyType_Ready(&PyTabprmType) < 0) {
    return -1;
  }

  Py_INCREF(&PyTabprmType);

  PyModule_AddObject(m, "Tabprm", (PyObject *)&PyTabprmType);

  tab_errexc[0] = NULL;                         /* Success */
  tab_errexc[1] = &PyExc_MemoryError;           /* Null wcsprm pointer passed */
  tab_errexc[2] = &PyExc_MemoryError;           /* Memory allocation failed */
  tab_errexc[3] = &WcsExc_InvalidTabularParameters;  /* Invalid tabular parameters */
  tab_errexc[4] = &WcsExc_InvalidCoordinate; /* One or more of the x coordinates were invalid */
  tab_errexc[5] = &WcsExc_InvalidCoordinate; /* One or more of the world coordinates were invalid */

  return 0;
}
