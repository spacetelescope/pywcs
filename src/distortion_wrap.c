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

/*
 MGDTODO: It's possible that we want to do the pipelining of SIP +
 Paper IV distortion + WCS in Python.  In which case, it's much easier
 to drop the Distortion class (which does the whole pipeline of Paper
 IV, Figure 1), in favor of just doing the first box in that figure.

 The complication is that it is hard to insert steps into the wcslib
 processing, so CQDIS becomes harder (or impossible) to handle.  For
 now, we are assuming we don't need CQDIS.  If that really is the
 case, the whole Distortion class (but not DistortionLookupTable)
 could be removed.
*/

/* util.h must be imported first */
#include "util.h"

#include "distortion.h"
#include "docstrings.h"
#include "str_list_proxy.h"

#include <structmember.h> /* From Python */

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
    if (distortion_lookup_t_init(&self->x))
      return NULL;
    self->py_data = NULL;
  }
  return (PyObject*)self;
}

static int
PyDistLookup_init(PyDistLookup* self, PyObject* args, PyObject* kwds) {
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
  if (array_obj == NULL)
    return -1;

  self->py_data = array_obj;
  self->x.naxis[0] = PyArray_DIM(array_obj, 0);
  self->x.naxis[1] = PyArray_DIM(array_obj, 1);
  self->x.data = (float *)PyArray_DATA(array_obj);

  return 0;
}

static PyObject*
PyDistLookup_get_cdelt(PyDistLookup* self, void* closure) {
  Py_ssize_t naxis = 2;

  return get_double_array("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_cdelt(PyDistLookup* self, PyObject* value, void* closure) {
  Py_ssize_t naxis = 2;

  return set_double_array("cdelt", value, 1, &naxis, self->x.cdelt);
}

static PyObject*
PyDistLookup_get_crpix(PyDistLookup* self, void* closure) {
  Py_ssize_t naxis = 2;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_crpix(PyDistLookup* self, PyObject* value, void* closure) {
  Py_ssize_t naxis = 2;

  return set_double_array("crpix", value, 1, &naxis, self->x.crpix);
}

static PyObject*
PyDistLookup_get_crval(PyDistLookup* self, void* closure) {
  Py_ssize_t naxis = 2;

  return get_double_array("crval", self->x.crval, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_crval(PyDistLookup* self, PyObject* value, void* closure) {
  Py_ssize_t naxis = 2;

  return set_double_array("crval", value, 1, &naxis, self->x.crval);
}

static PyObject*
PyDistLookup_get_data(PyDistLookup* self, void* closure) {
  Py_INCREF(self->py_data);
  return (PyObject*)self->py_data;
}

static int
PyDistLookup_set_data(PyDistLookup* self, PyObject* value, void* closure) {
  PyArrayObject* value_array = NULL;

  if (value == NULL) {
    Py_XDECREF(self->py_data);
    self->py_data = NULL;
    self->x.data = NULL;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, PyArray_FLOAT32, 2, 2);

  if (value_array == NULL)
    return -1;

  Py_XDECREF(self->py_data);

  self->py_data = value_array;
  self->x.naxis[0] = PyArray_DIM(value_array, 0);
  self->x.naxis[1] = PyArray_DIM(value_array, 1);
  self->x.data = (float *)PyArray_DATA(value_array);

  return 0;
}

static PyObject*
PyDistLookup_get_offset(PyDistLookup* self, PyObject* args, PyObject* kwds) {
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

static PyMemberDef PyDistLookup_members[] = {
  {NULL}
};

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
    if (distortion_t_init(&self->x)) {
      PyErr_SetString(PyExc_MemoryError,
                      "Could not initialize wcsprm object");
      return NULL;
    }

    for (i = 0; i < NAXES; ++i) {
      self->py_pre_dist[i] = NULL;
      self->py_post_dist[i] = NULL;
    }
  }
  return (PyObject*)self;
}

static int
PyDistortion_init(PyDistortion* self, PyObject* args, PyObject* kwds) {
  if (PySequence_Size(args) != 0) {
    PyErr_SetString(PyExc_IndexError, "Distortion.__init__ takes exactly 0 arguments");
    return -1;
  }

  return 0;
}

static PyObject*
PyDistortion_get_cd(PyDistortion* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x.has_pc == 1) {
    PyErr_SetString(PyExc_AttributeError, "No cd is present.");
    return NULL;
  }

  return get_double_array("cd", self->x.pc, 2, dims, (PyObject*)self);
}

static int
PyDistortion_set_cd(PyDistortion* self, PyObject* value, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (value == NULL) { /* deletion */
    self->x.has_pc = 1;
    return 0;
  }

  if (set_double_array("cd", value, 2, dims, self->x.pc)) {
    return -1;
  }

  self->x.has_pc = 0;

  return 0;
}

static PyObject*
PyDistortion_get_cdelt(PyDistortion* self, void* closure) {
  Py_ssize_t naxis = 2;

  if (self->x.has_pc == 0) {
    PyErr_SetString(PyExc_AttributeError, "No cdelt/pc is present.");
    return NULL;
  }

  return get_double_array("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
}

static int
PyDistortion_set_cdelt(PyDistortion* self, PyObject* value, void* closure) {
  Py_ssize_t naxis = 2;

  if (value == NULL) { /* deletion */
    self->x.has_pc = 0;
    return 0;
  }

  if (set_double_array("cdelt", value, 1, &naxis, self->x.cdelt)) {
    return -1;
  }

  self->x.has_pc = 1;

  return 0;
}

static inline PyObject*
get_lookup_tables(const char* propname, PyObject **tables) {
  PyObject* tables0 = tables[0];
  PyObject* tables1 = tables[1];
  if (tables0 == NULL)
    tables0 = Py_None;
  if (tables1 == NULL)
    tables1 = Py_None;
  return Py_BuildValue("OO", tables0, tables1);
}

static inline int
set_lookup_tables(const char* propname, PyObject* value,
                  PyObject** py_tables,
                  struct distortion_lookup_t* tables[]) {
  Py_ssize_t i;
  PyObject* subvalue = NULL;

  if (check_delete(propname, value)) {
    return -1;
  }

  if (!PySequence_Check(value) || PySequence_Size(value) != 2) {
    PyErr_Format(PyExc_TypeError, "'%s' must be a 2-length sequence", propname);
    return -1;
  }

  for (i = 0; i < NAXES; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (!(subvalue == Py_None || PyObject_TypeCheck(subvalue, &PyDistLookupType))) {
      Py_XDECREF(subvalue);
      PyErr_Format(PyExc_TypeError, "'%s' must be a 2-length sequence of DistortionLookupTable instances or None.", propname);
      return -1;
    }

    Py_DECREF(subvalue);
  }

  for (i = 0; i < NAXES; ++i) {
    subvalue = PySequence_GetItem(value, i);
    Py_XDECREF(py_tables[i]);
    if (subvalue == Py_None) {
      py_tables[i] = NULL;
      tables[i] = NULL;
    } else {
      py_tables[i] = subvalue;
      tables[i] = &(((PyDistLookup *)subvalue)->x);
    }
  }

  return 0;
}

static PyObject*
PyDistortion_get_cpdis(PyDistortion* self, void* closure) {
  return get_lookup_tables("cpdis", self->py_pre_dist);
}

static int
PyDistortion_set_cpdis(PyDistortion* self, PyObject* value, void* closure) {
  return set_lookup_tables("cpdis", value, self->py_pre_dist, self->x.pre_dist);
}

static PyObject*
PyDistortion_get_cqdis(PyDistortion* self, void* closure) {
  return get_lookup_tables("cqdis", self->py_post_dist);
}

static int
PyDistortion_set_cqdis(PyDistortion* self, PyObject* value, void* closure) {
  return set_lookup_tables("cqdis", value, self->py_post_dist, self->x.post_dist);
}

static PyObject*
PyDistortion_get_crpix(PyDistortion* self, void* closure) {
  Py_ssize_t naxis = 2;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static int
PyDistortion_set_crpix(PyDistortion* self, PyObject* value, void* closure) {
  Py_ssize_t naxis = 2;

  return set_double_array("crpix", value, 1, &naxis, self->x.crpix);
}

static PyObject*
PyDistortion_get_crval(PyDistortion* self, void* closure) {
  Py_ssize_t naxis = 2;

  if (is_null(self->x.wcs.crval)) {
    return NULL;
  }

  return get_double_array("crval", self->x.wcs.crval, 1, &naxis, (PyObject*)self);
}

static int
PyDistortion_set_crval(PyDistortion* self, PyObject* value, void* closure) {
  Py_ssize_t naxis = 2;

  if (is_null(self->x.wcs.crval)) {
    return -1;
  }

  return set_double_array("crval", value, 1, &naxis, self->x.wcs.crval);
}

static PyObject*
PyDistortion_get_ctype(PyDistortion* self, void* closure) {
  if (is_null(self->x.wcs.ctype)) {
    return NULL;
  }

  return get_str_list("ctype", self->x.wcs.ctype, NAXES, (PyObject*)self);
}

static int
PyDistortion_set_ctype(PyDistortion* self, PyObject* value, void* closure) {
  if (is_null(self->x.wcs.ctype)) {
    return -1;
  }

  return set_str_list("ctype", value, NAXES, 0, self->x.wcs.ctype);
}

static PyObject*
PyDistortion_get_pc(PyDistortion* self, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (self->x.has_pc == 0) {
    PyErr_SetString(PyExc_AttributeError, "No pc is present.");
    return NULL;
  }

  return get_double_array("pc", self->x.pc, 2, dims, (PyObject*)self);
}

static int
PyDistortion_set_pc(PyDistortion* self, PyObject* value, void* closure) {
  const npy_intp dims[2] = {2, 2};

  if (value == NULL) { /* deletion */
    self->x.has_pc = 0;
    return 0;
  }

  if (set_double_array("pc", value, 2, dims, self->x.pc)) {
    return -1;
  }

  self->x.has_pc = 1;

  return 0;
}

static PyObject*
PyDistortion_has_pc(PyDistortion* self, PyObject* arg) {
  if (self->x.has_pc) {
    Py_INCREF(Py_True);
    return Py_True;
  } else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}

static PyObject*
PyDistortion_get_pv(PyDistortion* self, PyObject* args, PyObject* kwds) {
  if (self->x.wcs.pv == NULL) {
    PyErr_SetString(PyExc_AssertionError, "No PVi_ma records present.");
    return NULL;
  }

  return get_pvcards("pv", self->x.wcs.pv, self->x.wcs.npv);
}

static PyObject*
PyDistortion_set_pv(PyDistortion* self, PyObject* arg) {
  if (is_null(self->x.wcs.ps)) {
    return NULL;
  }

  if (set_pvcards("pv", arg, &self->x.wcs.pv, &self->x.wcs.npv, &self->x.wcs.npvmax)) {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
PyDistortion_p2s_generic(PyDistortion* self, PyObject* arg, int do_shift) {
  PyArrayObject* pixcrd = NULL;
  PyArrayObject* world = NULL;
  int status = 0;
  unsigned int i = 0;

  for (i = 0; i < NAXES; ++i) {
    if (self->x.pre_dist[i] && self->x.pre_dist[i]->data == NULL) {
      PyErr_Format(PyExc_ValueError, "CPDIS DistortionLookupTable %d does not actually have a data array",
                   i);
      return NULL;
    }
  }

  for (i = 0; i < NAXES; ++i) {
    if (self->x.post_dist[i] && self->x.post_dist[i]->data == NULL) {
      PyErr_Format(PyExc_ValueError, "CQDIS DistortionLookupTable %d does not actually have a data array",
                   i);
      return NULL;
    }
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(arg, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  if (PyArray_DIM(pixcrd, 0) != NAXES) {
    PyErr_SetString(PyExc_ValueError, "Input pixel array must be of size (2, n)");
    status = -1;
    goto __PyDistortion_p2s_exit;
  }

  /* Now allocate an array for the results */
  world = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (world == NULL) {
    goto __PyDistortion_p2s_exit;
  }

  if (do_shift)
    offset_array(pixcrd, 1.0);
  wcsprm_python2c(&self->x.wcs);
  status = distortion_pipeline(&self->x,
                               PyArray_DIM(pixcrd, 1),
                               PyArray_DATA(pixcrd),
                               PyArray_DATA(world));
  wcsprm_c2python(&self->x.wcs);
  if (do_shift)
    offset_array(pixcrd, -1.0);

 __PyDistortion_p2s_exit:
  Py_XDECREF(pixcrd);

  if (status == 0 || status == 8) {
    return (PyObject*)world;
  }

  Py_XDECREF(world);

  if (status > 0 && status < WCS_ERRMSG_MAX) {
    PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
  } else if (status != -1) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
  }

  return NULL;
}

static PyObject*
PyDistortion_p2s(PyDistortion* self, PyObject* arg) {
  return PyDistortion_p2s_generic(self, arg, 1);
}

static PyObject*
PyDistortion_p2s_fits(PyDistortion* self, PyObject* arg) {
  return PyDistortion_p2s_generic(self, arg, 0);
}

static PyMemberDef PyDistortion_members[] = {
  {NULL}
};

static PyGetSetDef PyDistortion_getset[] = {
  {"cd",    (getter)PyDistortion_get_cd,    (setter)PyDistortion_set_cd,    (char *)doc_distortion_cd},
  {"cdelt", (getter)PyDistortion_get_cdelt, (setter)PyDistortion_set_cdelt, (char *)doc_cdelt},
  {"cpdis", (getter)PyDistortion_get_cpdis, (setter)PyDistortion_set_cpdis, (char *)doc_cpdis},
  {"cqdis", (getter)PyDistortion_get_cqdis, (setter)PyDistortion_set_cqdis, (char *)doc_cqdis},
  {"crpix", (getter)PyDistortion_get_crpix, (setter)PyDistortion_set_crpix, (char *)doc_distortion_crpix},
  {"crval", (getter)PyDistortion_get_crval, (setter)PyDistortion_set_crval, (char *)doc_crval},
  {"ctype", (getter)PyDistortion_get_ctype, (setter)PyDistortion_set_ctype, (char *)doc_ctype},
  {"pc",    (getter)PyDistortion_get_pc,    (setter)PyDistortion_set_pc,    (char *)doc_pc},
  {NULL}
};

static PyMethodDef PyDistortion_methods[] = {
  {"get_pv",      (PyCFunction)PyDistortion_get_pv,   METH_NOARGS, doc_get_pv},
  {"has_pc",      (PyCFunction)PyDistortion_has_pc,   METH_NOARGS, doc_distortion_has_pc},
  {"p2s",         (PyCFunction)PyDistortion_p2s,      METH_O,      doc_distortion_p2s},
  {"p2s_fits",    (PyCFunction)PyDistortion_p2s_fits, METH_O,      doc_distortion_p2s_fits},
  {"pixel2world", (PyCFunction)PyDistortion_p2s,      METH_O,      doc_distortion_pixel2world}, /* alias for p2s */
  {"pixel2world_fits", (PyCFunction)PyDistortion_p2s, METH_O,      doc_distortion_pixel2world_fits}, /* alias for p2s_fits */
  {"set_pv",      (PyCFunction)PyDistortion_set_pv,   METH_O,      doc_set_pv},
  {NULL}
};

static PyTypeObject PyDistortionType = {
  PyObject_HEAD_INIT(NULL)
  0,                          /*ob_size*/
  "pywcs.Distortion",         /*tp_name*/
  sizeof(PyDistortion),       /*tp_basicsize*/
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
  doc_Distortion,             /* tp_doc */
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

PyObject*
PyWcs_do_distortion(PyObject* self, PyObject* args, PyObject* kwds) {
  PyObject* lookups_obj = NULL;
  PyObject* lookup_obj = NULL;
  PyDistLookup* lookup = NULL;
  const struct distortion_lookup_t* lookups[] = { NULL, NULL };
  PyObject* pix_obj = NULL;
  PyArrayObject* pix_array = NULL;
  PyArrayObject* foc_array = NULL;
  Py_ssize_t i = 0;
  PyObject* result = NULL;

  if (!PyArg_ParseTuple(args, "OO:do_distortion",
                        &lookups_obj, &pix_obj)) {
    return NULL;
  }

  if (!PySequence_Check(lookups_obj) || PySequence_Size(lookups_obj) != 2) {
    PyErr_SetString(PyExc_ValueError, "First arg must be a 2-length sequence of DistortionLookupTable objects.");
    return NULL;
  }

  for (i = 0; i < 2; ++i) {
    lookup_obj = PySequence_GetItem(lookups_obj, i);
    if (!lookup_obj) {
      return NULL;
    }

    if (!PyObject_TypeCheck(lookup_obj, &PyDistLookupType)) {
      PyErr_SetString(PyExc_ValueError, "First arg must be a 2-length sequence of DistortionLookupTable objects.");
      Py_XDECREF(lookup_obj);
      return NULL;
    }

    lookup = (PyDistLookup*)lookup_obj;
    lookups[i] = &(lookup->x);
    Py_DECREF(lookup_obj);
  }

  pix_array = (PyArrayObject*)PyArray_ContiguousFromAny(pix_obj, PyArray_DOUBLE, 2, 2);
  if (pix_array == NULL) {
    goto _exit;
  }

  if (PyArray_DIM(pix_array, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto _exit;
  }

  foc_array = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pix_array), PyArray_DOUBLE);
  if (foc_array == NULL) {
    goto _exit;
  }

  if (do_distortion(2, lookups, PyArray_DIM(pix_array, 0),
                    PyArray_DATA(pix_array), PyArray_DATA(foc_array))) {
    PyErr_SetString(PyExc_RuntimeError, "Error performing distortion");
    Py_XDECREF(foc_array);
    goto _exit;
  }

  result = (PyObject*)foc_array;

 _exit:

  Py_XDECREF(pix_array);

  return result;
}

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
