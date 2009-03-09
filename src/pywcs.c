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

/* util.h must be imported first */
#include "util.h"
#include "wcslib_wrap.h"
#include "distortion_wrap.h"
#include "sip_wrap.h"
#include "pipeline.h"
#include "docstrings.h"

#include <structmember.h> /* from Python */

/***************************************************************************
 * Pywcs type
 ***************************************************************************/

static PyTypeObject PyWcsType;

typedef struct {
  PyObject_HEAD
  pipeline_t x;
  /*@null@*/ /*@shared@*/ PyObject* py_sip;
  /*@shared@*/ PyObject*            py_distortion_lookup[2];
  /*@null@*/ /*@shared@*/ PyObject* py_wcsprm;
} PyWcs;

static int _setup_pywcs_type(PyObject* m);


/***************************************************************************
 * PyWcs methods
 */

static int
PyWcs_traverse(
    PyWcs* self,
    visitproc visit,
    void* arg) {

  Py_VISIT(self->py_sip);
  Py_VISIT(self->py_distortion_lookup[0]);
  Py_VISIT(self->py_distortion_lookup[1]);
  Py_VISIT(self->py_wcsprm);

  return 0;
}

static int
PyWcs_clear(
    PyWcs* self) {

  PyObject* tmp;

  tmp = self->py_sip;
  self->py_sip = NULL;
  Py_XDECREF(tmp);

  tmp = self->py_distortion_lookup[0];
  self->py_distortion_lookup[0] = NULL;
  Py_XDECREF(tmp);

  tmp = self->py_distortion_lookup[1];
  self->py_distortion_lookup[1] = NULL;
  Py_XDECREF(tmp);

  tmp = self->py_wcsprm;
  self->py_wcsprm = NULL;
  Py_XDECREF(tmp);

  return 0;
}

static void
PyWcs_dealloc(
    PyWcs* self) {

  int ignored;
  ignored = PyWcs_clear(self);
  pipeline_free(&self->x);
  self->ob_type->tp_free((PyObject*)self);
}

/*@null@*/ static PyObject *
PyWcs_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyWcs* self;
  self = (PyWcs*)type->tp_alloc(type, 0);
  if (self != NULL) {
    pipeline_clear(&self->x);
    self->py_sip = NULL;
    self->py_distortion_lookup[0] = NULL;
    self->py_distortion_lookup[1] = NULL;
    self->py_wcsprm = NULL;
  }
  return (PyObject*)self;
}

static int
PyWcs_init(
    PyWcs* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  size_t    i;
  PyObject* py_sip;
  PyObject* py_wcsprm;
  PyObject* py_distortion_lookup[2];

  if (!PyArg_ParseTuple(args, "O(OO)O:Wcs.__init__",
                        &py_sip,
                        &py_distortion_lookup[0],
                        &py_distortion_lookup[1],
                        &py_wcsprm)) {
    return -1;
  }

  /* Check and set SIP */
  if (py_sip != NULL && py_sip != Py_None) {
    if (!PyObject_TypeCheck(py_sip, &PySipType)) {
      PyErr_SetString(PyExc_TypeError,
                      "Arg 1 must be Sip object");
      return -1;
    }

    self->py_sip = py_sip;
    self->x.sip = &(((PySip*)py_sip)->x);
  }

  /* Check and set Distortion lookup tables */
  for (i = 0; i < 2; ++i) {
    if (py_distortion_lookup[i] != NULL && py_distortion_lookup[i] != Py_None) {
      if (!PyObject_TypeCheck(py_distortion_lookup[i], &PyDistLookupType)) {
        PyErr_SetString(PyExc_TypeError,
                        "Arg 2 must be a pair of DistortionLookupTable or None objects");
        return -1;
      }

      self->py_distortion_lookup[i] = py_distortion_lookup[i];
      self->x.cpdis[i] = &(((PyDistLookup*)py_distortion_lookup[i])->x);
    }
  }

  /* Set and lookup Wcsprm object */
  if (py_wcsprm != NULL && py_wcsprm != Py_None) {
    if (!PyObject_TypeCheck(py_wcsprm, &PyWcsprmType)) {
      PyErr_SetString(PyExc_TypeError,
                      "Arg 3 must be Wcsprm object");
      return -1;
    }

    self->py_wcsprm = py_wcsprm;
    self->x.wcs = &(((PyWcsprm*)py_wcsprm)->x);
  }

  Py_XINCREF(self->py_sip);
  Py_XINCREF(self->py_distortion_lookup[0]);
  Py_XINCREF(self->py_distortion_lookup[1]);
  Py_XINCREF(self->py_wcsprm);

  return 0;
}

/*@null@*/ static PyObject*
PyWcs_all_pix2sky(
    PyWcs* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* pixcrd     = NULL;
  PyArrayObject* world      = NULL;
  int            status     = -1;
  const char*    keywords[] = {
    "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:all_pix2sky", (char **)keywords,
                                   &pixcrd_obj, &origin)) {
    return NULL;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(pixcrd_obj, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  world = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (world == NULL) {
    goto exit;
  }


  /* Make the call */
  preoffset_array(pixcrd, origin);
  wcsprm_python2c(self->x.wcs);
  status = pipeline_all_pixel2world(&self->x,
                                    (unsigned int)PyArray_DIM(pixcrd, 0),
                                    (unsigned int)PyArray_DIM(pixcrd, 1),
                                    (double*)PyArray_DATA(pixcrd),
                                    (double*)PyArray_DATA(world));
  wcsprm_c2python(self->x.wcs);
  unoffset_array(pixcrd, origin);

 exit:
  Py_XDECREF(pixcrd);

  if (status == 0 || status == 8) {
    return (PyObject*)world;
  } else {
    Py_DECREF(world);
    if (status > 0 && status < WCS_ERRMSG_MAX) {
      PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
      return NULL;
    } else if (status == -1) {
      return NULL;
    } else {
      PyErr_SetString(PyExc_RuntimeError,
                      "Unknown error occurred.  Something is seriously wrong.");
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
PyWcs_p4_pix2foc(
    PyWcs* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* pixcrd     = NULL;
  PyArrayObject* foccrd     = NULL;
  int            status     = -1;
  const char*    keywords[] = {
    "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:p4_pix2foc", (char **)keywords,
                                   &pixcrd_obj, &origin)) {
    return NULL;
  }

  if (self->x.cpdis[0] == NULL && self->x.cpdis[1] == NULL) {
    Py_INCREF(pixcrd_obj);
    return pixcrd_obj;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(pixcrd_obj, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  if (PyArray_DIM(pixcrd, 1) != NAXES) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto exit;
  }

  foccrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (foccrd == NULL) {
    status = 2;
    goto exit;
  }

  preoffset_array(pixcrd, origin);
  status = p4_pix2foc(2, (void *)self->x.cpdis,
                      (unsigned int)PyArray_DIM(pixcrd, 0),
                      (double*)PyArray_DATA(pixcrd),
                      (double*)PyArray_DATA(foccrd));
  unoffset_array(pixcrd, origin);

 exit:

  Py_XDECREF(pixcrd);

  if (status == 0) {
    return (PyObject*)foccrd;
  } else {
    Py_XDECREF(foccrd);
    if (status > 0 && status < WCS_ERRMSG_MAX) {
      PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
      return NULL;
    } else if (status == -1) {
      return NULL;
    } else {
      PyErr_SetString(PyExc_RuntimeError,
                      "Unknown error occurred.  Something is seriously wrong.");
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
PyWcs_pix2foc(
    PyWcs* self,
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

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(pixcrd_obj, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  if (PyArray_DIM(pixcrd, 1) != NAXES) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto _exit;
  }

  foccrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (foccrd == NULL) {
    goto _exit;
  }

  preoffset_array(pixcrd, origin);
  status = pipeline_pix2foc(&self->x,
                            (unsigned int)PyArray_DIM(pixcrd, 0),
                            (unsigned int)PyArray_DIM(pixcrd, 1),
                            (double*)PyArray_DATA(pixcrd),
                            (double*)PyArray_DATA(foccrd));
  unoffset_array(pixcrd, origin);

 _exit:

  Py_XDECREF(pixcrd);

  if (status == 0) {
    return (PyObject*)foccrd;
  } else {
    Py_XDECREF(foccrd);
    if (status > 0 && status < WCS_ERRMSG_MAX) {
      PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
      return NULL;
    } else if (status == -1) {
      return NULL;
    } else {
      PyErr_SetString(PyExc_RuntimeError,
                      "Unknown error occurred.  Something is seriously wrong.");
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
PyWcs_get_wcs(
    PyWcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_wcsprm) {
    Py_INCREF(self->py_wcsprm);
    return self->py_wcsprm;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
PyWcs_set_wcs(
    PyWcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_XDECREF(self->py_wcsprm);
  self->py_wcsprm = NULL;
  self->x.wcs = NULL;

  if (value != Py_None) {
    if (!PyObject_TypeCheck(value, &PyWcsprmType)) {
      PyErr_SetString(PyExc_TypeError,
                      "wcs must be Wcsprm object");
      return -1;
    }

    Py_INCREF(value);
    self->py_wcsprm = value;
    self->x.wcs = &(((PyWcsprm*)value)->x);
  }

  return 0;
}

static PyObject*
PyWcs_get_cpdis1(
    PyWcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_distortion_lookup[0]) {
    Py_INCREF(self->py_distortion_lookup[0]);
    return self->py_distortion_lookup[0];
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
PyWcs_set_cpdis1(
    PyWcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_XDECREF(self->py_distortion_lookup[0]);
  self->py_distortion_lookup[0] = NULL;
  self->x.cpdis[0] = NULL;

  if (value != Py_None) {
    if (!PyObject_TypeCheck(value, &PyDistLookupType)) {
      PyErr_SetString(PyExc_TypeError,
                      "cpdis1 must be DistortionLookupTable object");
      return -1;
    }

    Py_INCREF(value);
    self->py_distortion_lookup[0] = value;
    self->x.cpdis[0] = &(((PyDistLookup*)value)->x);
  }

  return 0;
}

/*@shared@*/ static PyObject*
PyWcs_get_cpdis2(
    PyWcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_distortion_lookup[1]) {
    Py_INCREF(self->py_distortion_lookup[1]);
    return self->py_distortion_lookup[1];
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
PyWcs_set_cpdis2(
    PyWcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_XDECREF(self->py_distortion_lookup[1]);
  self->py_distortion_lookup[1] = NULL;
  self->x.cpdis[1] = NULL;

  if (value != Py_None) {
    if (!PyObject_TypeCheck(value, &PyDistLookupType)) {
      PyErr_SetString(PyExc_TypeError,
                      "cpdis2 must be DistortionLookupTable object");
      return -1;
    }

    Py_INCREF(value);
    self->py_distortion_lookup[1] = value;
    self->x.cpdis[1] = &(((PyDistLookup*)value)->x);
  }

  return 0;
}

/*@shared@*/ static PyObject*
PyWcs_get_sip(
    PyWcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_sip) {
    Py_INCREF(self->py_sip);
    return self->py_sip;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
PyWcs_set_sip(
    PyWcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_XDECREF(self->py_sip);
  self->py_sip = NULL;
  self->x.sip = NULL;

  if (value != Py_None) {
    if (!PyObject_TypeCheck(value, &PySipType)) {
      PyErr_SetString(PyExc_TypeError,
                      "sip must be Sip object");
      return -1;
    }

    Py_INCREF(value);
    self->py_sip = value;
    self->x.sip = &(((PySip*)value)->x);
  }

  return 0;
}

static PyObject*
PyWcs___copy__(
    PyWcs* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyObject* copy = NULL;

  copy = PyWcs_new(&PyWcsType, NULL, NULL);
  if (copy == NULL) {
    return NULL;
  }

  if (self->py_sip) {
    PyWcs_set_sip((PyWcs*)copy, self->py_sip, NULL);
  }

  if (self->py_distortion_lookup[0]) {
    PyWcs_set_cpdis1((PyWcs*)copy, self->py_distortion_lookup[0], NULL);
  }

  if (self->py_distortion_lookup[1]) {
    PyWcs_set_cpdis2((PyWcs*)copy, self->py_distortion_lookup[1], NULL);
  }

  if (self->py_wcsprm) {
    PyWcs_set_wcs((PyWcs*)copy, self->py_wcsprm, NULL);
  }

  return copy;
}

static PyObject*
PyWcs___deepcopy__(
    PyWcs* self,
    PyObject* memo,
    /*@unused@*/ PyObject* kwds) {

  PyObject* copy;
  PyObject* obj_copy;

  copy = PyWcs_new(&PyWcsType, NULL, NULL);
  if (copy == NULL) {
    return NULL;
  }

  if (self->py_sip) {
    obj_copy = get_deepcopy(self->py_sip, memo);
    if (obj_copy == NULL) {
      Py_DECREF(copy);
      return NULL;
    }
    PyWcs_set_sip((PyWcs*)copy, obj_copy, NULL);
    Py_DECREF(obj_copy);
  }

  if (self->py_distortion_lookup[0]) {
    obj_copy = get_deepcopy(self->py_distortion_lookup[0], memo);
    if (obj_copy == NULL) {
      Py_DECREF(copy);
      return NULL;
    }
    PyWcs_set_cpdis1((PyWcs*)copy, obj_copy, NULL);
    Py_DECREF(obj_copy);
  }

  if (self->py_distortion_lookup[1]) {
    obj_copy = get_deepcopy(self->py_distortion_lookup[1], memo);
    if (obj_copy == NULL) {
      Py_DECREF(copy);
      return NULL;
    }
    PyWcs_set_cpdis2((PyWcs*)copy, obj_copy, NULL);
    Py_DECREF(obj_copy);
  }

  if (self->py_wcsprm) {
    obj_copy = get_deepcopy(self->py_wcsprm, memo);
    if (obj_copy == NULL) {
      Py_DECREF(copy);
      return NULL;
    }
    PyWcs_set_wcs((PyWcs*)copy, obj_copy, NULL);
    Py_DECREF(obj_copy);
  }

  return copy;
}


/***************************************************************************
 * PyWcs definition structures
 */

static PyGetSetDef PyWcs_getset[] = {
  {"cpdis1", (getter)PyWcs_get_cpdis1, (setter)PyWcs_set_cpdis1, (char *)doc_cpdis1},
  {"cpdis2", (getter)PyWcs_get_cpdis2, (setter)PyWcs_set_cpdis2, (char *)doc_cpdis2},
  {"sip", (getter)PyWcs_get_sip, (setter)PyWcs_set_sip, (char *)doc_sip},
  {"wcs", (getter)PyWcs_get_wcs, (setter)PyWcs_set_wcs, (char *)doc_wcs},
  {NULL}
};

static PyMethodDef PyWcs_methods[] = {
  {"_all_pix2sky", (PyCFunction)PyWcs_all_pix2sky, METH_VARARGS, doc_all_pix2sky},
  {"__copy__", (PyCFunction)PyWcs___copy__, METH_NOARGS, NULL},
  {"__deepcopy__", (PyCFunction)PyWcs___deepcopy__, METH_O, NULL},
  {"_p4_pix2foc", (PyCFunction)PyWcs_p4_pix2foc, METH_VARARGS, doc_p4_pix2foc},
  {"_pix2foc", (PyCFunction)PyWcs_pix2foc, METH_VARARGS, doc_pix2foc},
  {NULL}
};

static PyTypeObject PyWcsType = {
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  "pywcs._Wcs",                 /*tp_name*/
  sizeof(PyWcs),                /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PyWcs_dealloc,    /*tp_dealloc*/
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
  doc_Wcs,                      /* tp_doc */
  (traverseproc)PyWcs_traverse, /* tp_traverse */
  (inquiry)PyWcs_clear,         /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PyWcs_methods,                /* tp_methods */
  0,                            /* tp_members */
  PyWcs_getset,                 /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc)PyWcs_init,         /* tp_init */
  0,                            /* tp_alloc */
  PyWcs_new,                    /* tp_new */
};


/***************************************************************************
 * Module-level
 ***************************************************************************/

int _setup_pywcs_type(
    PyObject* m) {

  if (PyType_Ready(&PyWcsType) < 0)
    return -1;

  Py_INCREF(&PyWcsType);
  return PyModule_AddObject(m, "_Wcs", (PyObject *)&PyWcsType);
}

struct module_state {

};

#if PY_MAJOR_VERSION >= 3
    #define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
    #define GETSTATE(m) (&_state)
    static struct module_state _state;
#endif

#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_pywcs",
        NULL,
        sizeof(struct module_state),
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
    };

    #define INITERROR return NULL

    PyObject *
    PyInit__pywcs(void)

#else
    #define INITERROR return

    void
    init_pywcs(void)
#endif

{
#if PY_MAJOR_VERSION >= 3
    PyObject *m = PyModule_Create(&moduledef);
#else
    PyObject *m = Py_InitModule3("_pywcs", NULL, NULL);
#endif

    if (m == NULL)
        INITERROR;

    import_array();

    if (_setup_str_list_proxy_type(m) ||
        _setup_wcsprm_type(m)         ||
        _setup_distortion_type(m)     ||
        _setup_sip_type(m)            ||
        _setup_pywcs_type(m)          ||
        _define_exceptions(m)) {
        Py_DECREF(m);
        INITERROR;
    }

    PyModule_AddObject(m, "__docformat__", PyString_FromString("epytext"));

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
