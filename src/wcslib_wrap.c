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

#include "wcslib_wrap.h"
#include <structmember.h> /* from Python */

#include <wcs.h>
#include <wcsfix.h>
#include <wcshdr.h>
#include <wcsmath.h>

#include "isnan.h"
#include "str_list_proxy.h"
#include "distortion.h"

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

enum e_altlin {
  has_pc = 1,
  has_cd = 2,
  has_crota = 4
};

static int
parse_unsafe_unit_conversion_spec(
    const char* arg) {

  int ctrl = 0;
  const char* p = NULL;

  p = arg;
  for (p = arg; *p != '\0'; ++p) {
    switch (*p) {
    case 's':
    case 'S':
      ctrl |= 1;
      break;
    case 'h':
    case 'H':
      ctrl |= 2;
      break;
    case 'd':
    case 'D':
      ctrl |= 4;
      break;
    default:
      break;
    }
  }

  return ctrl;
}

static int
is_valid_alt_key(
    const char* key) {

  if (key[1] != '\0' ||
      !(key[0] == ' ' ||
        (key[0] >= 'A' && key[0] <= 'Z'))) {
    PyErr_SetString(PyExc_ValueError, "key must be ' ' or 'A'-'Z'");
    return 0;
  }

  return 1;
}

/***************************************************************************
 * PyWcsprm methods
 */

static inline void
note_change(PyWcsprm* self) {
  self->x.flag = 0;
}

static void
PyWcsprm_dealloc(
    PyWcsprm* self) {

  int ignored;
  ignored = wcsfree(&self->x);
  self->ob_type->tp_free((PyObject*)self);
}

static PyWcsprm*
PyWcsprm_cnew(void) {
  PyWcsprm* self;
  self = (PyWcsprm*)(&PyWcsprmType)->tp_alloc(&PyWcsprmType, 0);
  return self;
}

static PyObject *
PyWcsprm_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyWcsprm* self;
  self = (PyWcsprm*)type->tp_alloc(type, 0);
  return (PyObject*)self;
}

static int
PyWcsprm_init(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

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
  const char*    keywords[]    = {"header", "key", "relax", "naxis", NULL};
  PyObject*      ignored       = NULL;
  int            ignored_int;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Osii:WCSBase.__init__",
                                  (char **)keywords, &header_obj, &key,
                                   &relax, &naxis)) {
    return -1;
  }

  if (header_obj == NULL || header_obj == Py_None) {
    if (relax || key[0] != ' ' || key[1] != '\0') {
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

    self->x.flag = -1;
    status = wcsini(1, naxis, &self->x);

    if (status != 0) {
      PyErr_SetString(PyExc_MemoryError,
                      "Could not initialize wcsprm object");
      return -1;
    }

    wcsprm_c2python(&self->x);

    return 0;
  } else { /* header != NULL */
    if (PyString_AsStringAndSize(header_obj, &header, &header_length)) {
      return -1;
    }

    if (relax) {
      relax = WCSHDR_all;
    }

    if (!is_valid_alt_key(key)) {
      return -1;
    }

    if (naxis != 2) {
      PyErr_SetString(PyExc_ValueError,
                      "naxis may not be provided if a header is provided.");
      return -1;
    }

    nkeyrec = header_length / 80;
    if (nkeyrec > 0x7fffffff) {
      return -1;
    }

    status = wcspih(header,
                    (int)nkeyrec,
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

    if (nwcs == 0) {
      PyErr_SetString(WcsExc_NoWcsKeywordsFound,
                      "No WCS keywords found in the given header");
      return -1;
    }

    /* Find the desired WCS */
    for (i = 0; i < nwcs; ++i) {
      if (wcs[i].alt[0] == key[0]) {
        break;
      }
    }

    if (i >= nwcs) {
      ignored_int = wcsvfree(&nwcs, &wcs);
      ignored = PyErr_Format(
          PyExc_KeyError,
          "No WCS with key '%s' was found in the given header",
          key);
      return -1;
    }

    self->x.flag = -1;
    if (wcscopy(1, wcs + i, &self->x) != 0) {
      ignored_int = wcsfree(&self->x);
      ignored_int = wcsvfree(&nwcs, &wcs);
      PyErr_SetString(PyExc_MemoryError,
                      "Could not initialize wcsprm object");
      return -1;
    }

    wcsprm_c2python(&self->x);
    ignored_int = wcsvfree(&nwcs, &wcs);
    return 0;
  }
}

/*@null@*/ static PyObject*
PyWcsprm_copy(
    PyWcsprm* self) {

  PyWcsprm*      copy      = NULL;
  int            status;

  copy = PyWcsprm_cnew();
  if (copy == NULL) {
    return NULL;
  }

  wcsprm_python2c(&self->x);
  status = wcscopy(1, &self->x, &copy->x);
  wcsprm_c2python(&self->x);

  if (status == 0) {
    wcsprm_c2python(&copy->x);
    return (PyObject*)copy;
  } else {
    Py_XDECREF(copy);
    if (status > 0 && status < WCS_ERRMSG_MAX) {
      PyErr_SetString(*wcs_errexc[status], wcscopy_errmsg[status]);
      return NULL;
    } else {
      PyErr_SetString(PyExc_RuntimeError,
                      "Unknown error occurred.  Something is seriously wrong.");
      return NULL;
    }
  }
}

static PyObject*
PyWcsprm_celfix(
    PyWcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = celfix(&self->x);
  wcsprm_c2python(&self->x);

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

/*@null@*/ static PyObject*
PyWcsprm_cylfix(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      naxis_obj   = NULL;
  PyArrayObject* naxis_array = NULL;
  int*           naxis       = NULL;
  int            status      = 0;
  const char*    keywords[]  = {"naxis", NULL};
  PyObject*      ignored     = NULL;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O:cylfix", (char **)keywords,
                                   &naxis_obj)) {
    return NULL;
  }

  if (naxis_obj != NULL) {
    naxis_array = (PyArrayObject*)PyArray_ContiguousFromAny(naxis_obj, 1, 1,
                                                            PyArray_INT);
    if (naxis_array == NULL) {
      return NULL;
    }
    if (PyArray_DIM(naxis_array, 0) != self->x.naxis) {
      ignored = PyErr_Format(
          PyExc_ValueError,
          "naxis must be same length as the number of axes of "
          "the Wcsprm object (%d).",
          self->x.naxis);
      Py_DECREF(naxis_array);
      return NULL;
    }
    naxis = (int*)PyArray_DATA(naxis_array);
  }

  wcsprm_python2c(&self->x);
  status = cylfix(naxis, &self->x);
  wcsprm_c2python(&self->x);

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
PyWcsprm_datfix(
    PyWcsprm* self) {

  int status = 0;

  status = datfix(&self->x);

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

/*@null@*/ static PyObject*
PyWcsprm_fix(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

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
  PyObject*      ignored         = NULL;

  struct message_map_entry {
    const char* name;
    const int index;
  };
  const struct message_map_entry message_map[NWCSFIX] = {
    {"datfix", DATFIX},
    {"unitfix", UNITFIX},
    {"celfix", CELFIX},
    {"spcfix", SPCFIX},
    {"cylfix", CYLFIX},
    {NULL}
  };
  const char* keywords[] = {"translate_units", "naxis", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|sO:fix", (char **)keywords,
                                   &translate_units, &naxis_obj)) {
    return NULL;
  }

  if (translate_units != NULL) {
    ctrl = parse_unsafe_unit_conversion_spec(translate_units);
  }

  if (naxis_obj != NULL) {
    naxis_array = (PyArrayObject*)PyArray_ContiguousFromAny(naxis_obj, 1, 1,
                                                            PyArray_INT);
    if (naxis_array == NULL) {
      return NULL;
    }
    if (PyArray_DIM(naxis_array, 0) != self->x.naxis) {
      ignored = PyErr_Format(
          PyExc_ValueError,
          "naxis must be same length as the number of axes of "
          "the Wcprm object (%d).",
          self->x.naxis);
      Py_DECREF(naxis_array);
      return NULL;
    }
    naxis = (int*)PyArray_DATA(naxis_array);
  }

  wcsprm_python2c(&self->x);
  status = wcsfix(ctrl, naxis, &self->x, stat);
  wcsprm_c2python(&self->x);

  /* We're done with this already, so deref now so we don't have to remember
     later */
  Py_XDECREF(naxis_array);

  result = PyDict_New();
  if (result == NULL) {
    return NULL;
  }

  for (i = 0; i < NWCSFIX; ++i) {
    msg_index = stat[message_map[i].index];
    if (msg_index >= 0 && msg_index < 11) {
      subresult = PyString_FromString(wcsfix_errmsg[msg_index]);
      if (subresult == NULL ||
          PyDict_SetItemString(result, message_map[i].name, subresult)) {
        Py_XDECREF(subresult);
        Py_XDECREF(result);
        return NULL;
      }
      Py_XDECREF(subresult);
    }
  }

  return result;
}

/*@null@*/ static PyObject*
PyWcsprm_get_ps(
    PyWcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  if (self->x.ps == NULL) {
    PyErr_SetString(PyExc_AssertionError, "No PSi_ma records present.");
    return NULL;
  }

  return get_pscards("ps", self->x.ps, self->x.nps);
}

/*@null@*/ static PyObject*
PyWcsprm_get_pv(
    PyWcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  if (self->x.pv == NULL) {
    PyErr_SetString(PyExc_AssertionError, "No PVi_ma records present.");
    return NULL;
  }

  return get_pvcards("pv", self->x.pv, self->x.npv);
}

static PyObject*
PyWcsprm_has_cdi_ja(
    PyWcsprm* self) {

  int result = 0;

  result = self->x.altlin & has_cd;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_has_crotaia(
    PyWcsprm* self) {

  int result = 0;

  result = self->x.altlin & has_crota;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_has_pci_ja(
    PyWcsprm* self) {

  int result = 0;

  result = self->x.altlin & has_pc;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_is_unity(
    PyWcsprm* self) {

  return PyBool_FromLong(self->x.lin.unity);
}

/*@null@*/ static PyObject*
PyWcsprm_mix_generic(
    PyWcsprm* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds,
    int do_shift) {

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
  PyObject*      ignored    = NULL;

  if (!PyArg_ParseTuple(args, "ii(dd)diOO:mix",
                       &mixpix, &mixcel, &vspan[0], &vspan[1],
                        &vstep, &viter, &world_obj, &pixcrd_obj)) {
    return NULL;
  }

  if (viter < 5 || viter > 10) {
    PyErr_SetString(PyExc_ValueError,
                    "viter must be in the range 5 - 10");
    goto exit;
  }

  world = (PyArrayObject*)PyArray_ContiguousFromAny
    (world_obj, PyArray_DOUBLE, 1, 1);
  if (world == NULL) {
    PyErr_SetString(PyExc_TypeError,
                    "Argument 6 (world) must be a 1-dimensional numpy array");
    goto exit;
  }
  if ((int)PyArray_DIM(world, 0) != self->x.naxis) {
    ignored = PyErr_Format(
        PyExc_TypeError,
        "Argument 6 (world) must be the same length as the number "
        "of axes (%d)",
        self->x.naxis);
    goto exit;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny
    (pixcrd_obj, PyArray_DOUBLE, 1, 1);
  if (pixcrd == NULL) {
    PyErr_SetString(PyExc_TypeError,
                    "Argument 7 (pixcrd) must be a 1-dimensional numpy array");
    goto exit;
  }
  if ((int)PyArray_DIM(pixcrd, 0) != self->x.naxis) {
    ignored = PyErr_Format(PyExc_TypeError,
                           "Argument 7 (pixcrd) must be the same length as the "
                           "number of axes (%d)",
                           self->x.naxis);
    goto exit;
  }

  if (mixpix < 1 || mixpix > 2) {
    PyErr_SetString(PyExc_ValueError,
                    "Argument 1 (mixpix) must specify a pixel coordinate "
                    "axis number");
    goto exit;
  }

  if (mixcel < 1 || mixcel > 2) {
    PyErr_SetString(PyExc_ValueError,
                    "Argument 2 (mixcel) must specify a celestial coordinate "
                    "axis number (1 for latitude, 2 for longitude)");
    goto exit;
  }

  /* Now we allocate a bunch of numpy arrays to store the
   * results in.
   */
  naxis = (Py_ssize_t)self->x.naxis;
  phi = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, PyArray_DOUBLE);
  if (theta == NULL) {
    status = 2;
    goto exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, PyArray_DOUBLE);
  if (imgcrd == NULL) {
    status = 2;
    goto exit;
  }

  /* Convert pixel coordinates to 1-based */
  if (do_shift) {
    offset_array(pixcrd, 1.0);
  }
  wcsprm_python2c(&self->x);
  status = wcsmix(&self->x,
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
  wcsprm_c2python(&self->x);
  /* Convert pixel coordinates back to 0-based) */
  if (do_shift) {
    offset_array(pixcrd, -1.0);
  }

  if (status == 0) {
    result = PyDict_New();
    if (result == NULL ||
        PyDict_SetItemString(result, "imgcrd", (PyObject*)imgcrd) ||
        PyDict_SetItemString(result, "phi", (PyObject*)phi) ||
        PyDict_SetItemString(result, "theta", (PyObject*)theta)) {
      status = 2;
    }
  }

 exit:
  Py_XDECREF(world);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(imgcrd);
  Py_XDECREF(pixcrd);

  if (status == 0) {
    return result;
  } else {
    Py_XDECREF(result);
    if (status == -1) {
      /* The error message has already been set */
      return NULL;
    } else if (status > 0 && status < WCS_ERRMSG_MAX) {
      PyErr_SetString(*wcs_errexc[status], wcsmix_errmsg[status]);
      return NULL;
    } else {
      PyErr_SetString(PyExc_RuntimeError,
                      "Unknown error occurred.  Something is seriously wrong.");
      return NULL;
    }
  }
}

static PyObject*
PyWcsprm_mix(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {
  return PyWcsprm_mix_generic(self, args, kwds, 1);
}

static PyObject*
PyWcsprm_mix_fits(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {
  return PyWcsprm_mix_generic(self, args, kwds, 0);
}

/*@null@*/ static PyObject*
PyWcsprm_p2s_generic(
    PyWcsprm* self,
    PyObject* arg,
    int do_shift) {

  PyArrayObject* pixcrd  = NULL;
  PyArrayObject* imgcrd  = NULL;
  PyArrayObject* phi     = NULL;
  PyArrayObject* theta   = NULL;
  PyArrayObject* world   = NULL;
  PyArrayObject* stat    = NULL;
  PyObject*      result  = NULL;
  int            status  = 0;

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny
    (arg, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  /* Now we allocate a bunch of numpy arrays to store the results in.
   */
  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (imgcrd == NULL) {
    status = 2;
    goto exit;
  }

  phi = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (theta == NULL) {
    status = 2;
    goto exit;
  }

  world = (PyArrayObject*)PyArray_SimpleNew
    (2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (world == NULL) {
    status = 2;
    goto exit;
  }

  stat = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(pixcrd), PyArray_INT);
  if (stat == NULL) {
    status = 2;
    goto exit;
  }

  /* Adjust pixel coordinates to be 1-based */
  if (do_shift) {
    offset_array(pixcrd, 1.0);
  }

  /* Make the call */
  wcsprm_python2c(&self->x);
  status = wcsp2s(&self->x,
                  (int)PyArray_DIM(pixcrd, 0),
                  (int)PyArray_DIM(pixcrd, 1),
                  (double*)PyArray_DATA(pixcrd),
                  (double*)PyArray_DATA(imgcrd),
                  (double*)PyArray_DATA(phi),
                  (double*)PyArray_DATA(theta),
                  (double*)PyArray_DATA(world),
                  (int*)PyArray_DATA(stat));
  wcsprm_c2python(&self->x);
  /* Adjust pixel coordinates back to 0-based */
  if (do_shift) {
    offset_array(pixcrd, -1.0);
  }

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

 exit:
  Py_XDECREF(pixcrd);
  Py_XDECREF(imgcrd);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(world);
  Py_XDECREF(stat);

  if (status == 0 || status == 8) {
    return result;
  } else {
    Py_XDECREF(result);
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
PyWcsprm_p2s(
    PyWcsprm* self,
    PyObject* arg) {

  return PyWcsprm_p2s_generic(self, arg, 1);
}

/*@null@*/ static PyObject*
PyWcsprm_p2s_fits(
    PyWcsprm* self,
    PyObject* arg) {

  return PyWcsprm_p2s_generic(self, arg, 0);
}

/*@null@*/ static PyObject*
PyWcsprm_s2p_generic(
    PyWcsprm* self,
    PyObject* arg,
    int do_shift) {

  PyArrayObject* world   = NULL;
  PyArrayObject* phi     = NULL;
  PyArrayObject* theta   = NULL;
  PyArrayObject* imgcrd  = NULL;
  PyArrayObject* pixcrd  = NULL;
  PyArrayObject* stat    = NULL;
  PyObject*      result  = NULL;
  int            status  = 0;

  world = (PyArrayObject*)PyArray_ContiguousFromAny
    (arg, PyArray_DOUBLE, 2, 2);
  if (world == NULL) {
    return NULL;
  }

  /* Now we allocate a bunch of numpy arrays to store the
   * results in.
   */
  phi = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(world), PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(world), PyArray_DOUBLE);
  if (phi == NULL) {
    status = 2;
    goto exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (2, PyArray_DIMS(world), PyArray_DOUBLE);
  if (theta == NULL) {
    status = 2;
    goto exit;
  }

  pixcrd = (PyArrayObject*)PyArray_SimpleNew
    (2, PyArray_DIMS(world), PyArray_DOUBLE);
  if (pixcrd == NULL) {
    status = 2;
    goto exit;
  }

  stat = (PyArrayObject*)PyArray_SimpleNew
    (1, PyArray_DIMS(world), PyArray_INT);
  if (stat == NULL) {
    status = 2;
    goto exit;
  }

  /* Adjust pixel coordinates to be zero-based */
  if (do_shift) {
    offset_array(pixcrd, 1.0);
  }

  /* Make the call */
  wcsprm_python2c(&self->x);
  status = wcss2p(&self->x,
                  (int)PyArray_DIM(world, 0),
                  (int)PyArray_DIM(world, 1),
                  (double*)PyArray_DATA(world),
                  (double*)PyArray_DATA(phi),
                  (double*)PyArray_DATA(theta),
                  (double*)PyArray_DATA(imgcrd),
                  (double*)PyArray_DATA(pixcrd),
                  (int*)PyArray_DATA(stat));
  wcsprm_c2python(&self->x);

  /* Adjust pixel coordinates to be zero-based */
  if (do_shift) {
    offset_array(pixcrd, -1.0);
  }

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

 exit:
  Py_XDECREF(pixcrd);
  Py_XDECREF(imgcrd);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(world);
  Py_XDECREF(stat);

  if (status == 0 || status == 9) {
    return result;
  } else {
    Py_XDECREF(result);
    if (status > 0 && status < WCS_ERRMSG_MAX) {
      PyErr_SetString(*wcs_errexc[status], wcsp2s_errmsg[status]);
      return NULL;
    } else {
      PyErr_SetString(PyExc_RuntimeError,
                      "Unknown error occurred.  Something is seriously wrong.");
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
PyWcsprm_s2p(
    PyWcsprm* self,
    PyObject* arg) {

  return PyWcsprm_s2p_generic(self, arg, 1);
}

/*@null@*/ static PyObject*
PyWcsprm_s2p_fits(
    PyWcsprm* self,
    PyObject* arg) {

  return PyWcsprm_s2p_generic(self, arg, 0);
}

static int
PyWcsprm_cset(
    PyWcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = wcsset(&self->x);
  wcsprm_c2python(&self->x);

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

/*@null@*/ static PyObject*
PyWcsprm_set(
    PyWcsprm* self) {

  if (PyWcsprm_cset(self)) {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
PyWcsprm_set_ps(
    PyWcsprm* self,
    PyObject* arg,
    /*@unused@*/ PyObject* kwds) {

  if (is_null(self->x.ps)) {
    return NULL;
  }

  if (set_pscards("ps", arg, &self->x.ps, &self->x.nps, &self->x.npsmax)) {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
PyWcsprm_set_pv(
    PyWcsprm* self,
    PyObject* arg,
    /*@unused@*/ PyObject* kwds) {

  if (is_null(self->x.pv)) {
    return NULL;
  }

  if (set_pvcards("pv", arg, &self->x.pv, &self->x.npv, &self->x.npvmax)) {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/* TODO: This is convenient for debugging for now -- but it's not very
 * Pythonic.  It should probably be hooked into __str__ or something.
 */
/*@null@*/ static PyObject*
PyWcsprm_print_contents(
    PyWcsprm* self) {

  int ignored;

  if (PyWcsprm_set(self) == NULL) {
    return NULL;
  }

  wcsprm_python2c(&self->x);
  ignored = wcsprt(&self->x);
  wcsprm_c2python(&self->x);

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
PyWcsprm_spcfix(
    PyWcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = spcfix(&self->x);
  wcsprm_c2python(&self->x);

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

/*@null@*/ static PyObject*
PyWcsprm_sptr(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int   i                = -1;
  char* py_ctype         = NULL;
  char  ctype[9];
  int   status           = 0;
  const char* keywords[] = {"i", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|i:sptr", (char **)keywords,
                                   &py_ctype, &i)) {
    return NULL;
  }

  if (strlen(py_ctype) > 8) {
    PyErr_SetString(PyExc_ValueError,
                    "ctype string has more than 8 characters.");
  }

  strncpy(ctype, py_ctype, 9);

  wcsprm_python2c(&self->x);
  status = wcssptr(&self->x, &i, ctype);
  wcsprm_c2python(&self->x);

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

/*@null@*/ static PyObject*
PyWcsprm_to_header(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int       relax        = 0;
  int       nkeyrec      = 0;
  char*     header       = NULL;
  int       status       = 0;
  PyObject* result       = NULL;
  const char* keywords[] = {"relax", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i:to_header",
                                   (char **)keywords, &relax)) {
    return NULL;
  }

  if (relax) {
    relax = -1;
  }

  wcsprm_python2c(&self->x);
  status = wcshdo(relax, &self->x, &nkeyrec, &header);
  wcsprm_c2python(&self->x);

  if (status != 0) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Unknown error occurred.  Something is seriously wrong.");
    free(header);
    return NULL;
  }

  /* Just return the raw header string.  PyFITS on the Python side will help
     to parse and use this information. */
  result = PyString_FromStringAndSize(header, (Py_ssize_t)nkeyrec * 80);

  free(header);

  return result;
}

/*@null@*/ static PyObject*
PyWcsprm_unitfix(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  char* translate_units  = NULL;
  int   ctrl             = 0;
  int   status           = 0;
  const char* keywords[] = {"translate_units", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|s:unitfix", (char **)keywords,
                                   &translate_units)) {
    return NULL;
  }

  if (translate_units != NULL) {
    ctrl = parse_unsafe_unit_conversion_spec(translate_units);
  }

  status = unitfix(ctrl, &self->x);

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
 * Member getters/setters (properties)
 */
/*@null@*/ static PyObject*
PyWcsprm_get_alt(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.alt)) {
    return NULL;
  }

  /* Force a null-termination of this single-character string */
  self->x.alt[1] = '\0';
  return get_string("alt", self->x.alt);
}

static int
PyWcsprm_set_alt(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  char value_string[2];

  if (is_null(self->x.alt)) {
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x.alt[0] = ' ';
    self->x.alt[1] = '\0';
    return 0;
  }

  if (set_string("propname", value, value_string, 2)) {
    return -1;
  }

  if (!is_valid_alt_key(value_string)) {
    return -1;
  }

  strncpy(self->x.alt, value_string, 2);

  note_change(self);

  return 0;
}

/*@null@*/ static PyObject*
PyWcsprm_get_cd(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {2, 2};

  if (is_null(self->x.cd)) {
    return NULL;
  }

  if ((self->x.altlin & has_cd) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No cd is present.");
    return NULL;
  }

  return get_double_array("cd", self->x.cd, 2, dims, (PyObject*)self);
}

static int
PyWcsprm_set_cd(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {2, 2};

  if (is_null(self->x.cd)) {
    return -1;
  }

  if (value == NULL) {
    self->x.altlin &= ~has_cd;
    return 0;
  }

  if (set_double_array("cd", value, 2, dims, self->x.cd)) {
    return -1;
  }

  self->x.altlin |= has_cd;

  note_change(self);

  return 0;
}

/*@null@*/ static PyObject*
PyWcsprm_get_cname(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.cname)) {
    return NULL;
  }

  return get_str_list("cname", self->x.cname, (Py_ssize_t)self->x.naxis, (PyObject*)self);
}

/*@null@*/ static int
PyWcsprm_set_cname(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {
  if (is_null(self->x.cname)) {
    return -1;
  }

  note_change(self);

  return set_str_list("cname", value, (Py_ssize_t)self->x.naxis, 0, self->x.cname);
}

/*@null@*/ static PyObject*
PyWcsprm_get_cdelt(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.cdelt)) {
    return NULL;
  }

  naxis = self->x.naxis;

  return get_double_array("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
}

/*@null@*/ static int
PyWcsprm_set_cdelt(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp dims;

  if (is_null(self->x.cdelt)) {
    return -1;
  }

  dims = (npy_int)self->x.naxis;

  note_change(self);

  return set_double_array("cdelt", value, 1, &dims, self->x.cdelt);
}

static PyObject*
PyWcsprm_get_cel_offset(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return PyBool_FromLong(self->x.cel.offset);
}

static int
PyWcsprm_set_cel_offset(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {
  return set_bool("cel_offset", value, &self->x.cel.offset);
}

static PyObject*
PyWcsprm_get_colnum(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("colnum", self->x.colnum);
}

static int
PyWcsprm_set_colnum(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_int("colnum", value, &self->x.colnum);
}

/*@null@*/ static PyObject*
PyWcsprm_get_colax(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.colax)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_int_array("colax", self->x.colax, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_colax(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 0;

  if (is_null(self->x.colax)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  note_change(self);

  return set_int_array("colax", value, 1, &naxis, self->x.colax);
}

/*@null@*/ static PyObject*
PyWcsprm_get_crder(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crder)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crder", self->x.crder, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_crder(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 0;

  if (is_null(self->x.crder)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  note_change(self);

  return set_double_array("crder", value, 1, &naxis, self->x.crder);
}

/*@null@*/ static PyObject*
PyWcsprm_get_crota(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crota)) {
    return NULL;
  }

  if ((self->x.altlin & has_crota) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No crota is present.");
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crota", self->x.crota, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_crota(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 0;

  if (is_null(self->x.crota)) {
    return -1;
  }

  if (value == NULL) { /* Deletion */
    self->x.altlin &= ~has_crota;
    return 0;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  if (set_double_array("crota", value, 1, &naxis, self->x.crota)) {
    return -1;
  }

  self->x.altlin |= has_crota;

  note_change(self);

  return 0;
}

/*@null@*/ static PyObject*
PyWcsprm_get_crpix(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crpix)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_crpix(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 0;

  if (is_null(self->x.crpix)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  note_change(self);

  return set_double_array("crpix", value, 1, &naxis, self->x.crpix);
}

/*@null@*/ static PyObject*
PyWcsprm_get_crval(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crval)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crval", self->x.crval, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_crval(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis;

  if (is_null(self->x.crval)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  note_change(self);

  return set_double_array("crval", value, 1, &naxis, self->x.crval);
}

/*@null@*/ static PyObject*
PyWcsprm_get_csyer(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis;

  if (is_null(self->x.csyer)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("csyer", self->x.csyer, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_csyer(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis;

  if (is_null(self->x.csyer)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  note_change(self);

  return set_double_array("csyer", value, 1, &naxis, self->x.csyer);
}

/*@null@*/ static PyObject*
PyWcsprm_get_ctype(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ctype)) {
    return NULL;
  }

  return get_str_list("ctype", self->x.ctype, self->x.naxis, (PyObject*)self);
}

static int
PyWcsprm_set_ctype(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ctype)) {
    return -1;
  }

  note_change(self);

  return set_str_list("ctype", value, (Py_ssize_t)self->x.naxis, 0, self->x.ctype);
}

static PyObject*
PyWcsprm_get_cubeface(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_bool("cubeface", self->x.cubeface);
}

static int
PyWcsprm_set_cubeface(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_bool("cubeface", value, &self->x.cubeface);
}

/*@null@*/ static PyObject*
PyWcsprm_get_cunit(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.cunit)) {
    return NULL;
  }

  return get_str_list("cunit", self->x.cunit, (Py_ssize_t)self->x.naxis, (PyObject*)self);
}

static int
PyWcsprm_set_cunit(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.cunit)) {
    return -1;
  }

  note_change(self);

  return set_str_list("cunit", value, (Py_ssize_t)self->x.naxis, 0, self->x.cunit);
}

/*@null@*/ static PyObject*
PyWcsprm_get_dateavg(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateavg)) {
    return NULL;
  }

  return get_string("dateavg", self->x.dateavg);
}

static int
PyWcsprm_set_dateavg(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateavg)) {
    return -1;
  }

  note_change(self);

  return set_string("dateavg", value, self->x.dateavg, 72);
}

/*@null@*/ static PyObject*
PyWcsprm_get_dateobs(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateobs)) {
    return NULL;
  }

  return get_string("dateobs", self->x.dateobs);
}

static int
PyWcsprm_set_dateobs(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateobs)) {
    return -1;
  }

  note_change(self);

  return set_string("dateobs", value, self->x.dateobs, 72);
}

static PyObject*
PyWcsprm_get_equinox(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("equinox", self->x.equinox);
}

static int
PyWcsprm_set_equinox(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.equinox = UNDEFINED;
    return 0;
  }

  note_change(self);

  return set_double("equinox", value, &self->x.equinox);
}

/*@null@*/ static PyObject*
PyWcsprm_get_imgpix_matrix(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {2, 2};

  if (is_null(self->x.lin.imgpix)) {
    return NULL;
  }

  return get_double_array("imgpix_matrix", self->x.lin.imgpix, 2, dims,
                          (PyObject*)self);
}

static PyObject*
PyWcsprm_get_lat(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("lat", self->x.lat);
}

static PyObject*
PyWcsprm_get_latpole(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("latpole", self->x.latpole);
}

static int
PyWcsprm_set_latpole(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_double("latpole", value, &self->x.latpole);
}

/*@null@*/ static PyObject*
PyWcsprm_get_lattyp(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.lattyp)) {
    return NULL;
  }

  return get_string("lattyp", self->x.lattyp);
}

static PyObject*
PyWcsprm_get_lng(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("lng", self->x.lng);
}

/*@null@*/ static PyObject*
PyWcsprm_get_lngtyp(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.lngtyp)) {
    return NULL;
  }

  return get_string("lngtyp", self->x.lngtyp);
}

static PyObject*
PyWcsprm_get_lonpole(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("lonpole", self->x.lonpole);
}

static int
PyWcsprm_set_lonpole(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_double("lonpole", value, &self->x.lonpole);
}

static PyObject*
PyWcsprm_get_mjdavg(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("mjdavg", self->x.mjdavg);
}

static int
PyWcsprm_set_mjdavg(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_double("mjdavg", value, &self->x.mjdavg);
}

static PyObject*
PyWcsprm_get_mjdobs(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("mjdobs", self->x.mjdobs);
}

static int
PyWcsprm_set_mjdobs(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_double("mjdobs", value, &self->x.mjdobs);
}

/*@null@*/ static PyObject*
PyWcsprm_get_name(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.wcsname)) {
    return NULL;
  }

  return get_string("name", self->x.wcsname);
}

static int
PyWcsprm_set_name(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.wcsname)) {
    return -1;
  }

  note_change(self);

  return set_string("name", value, self->x.wcsname, 72);
}

static PyObject*
PyWcsprm_get_naxis(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("naxis", self->x.naxis);
}

/*@null@*/ static PyObject*
PyWcsprm_get_obsgeo(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t size = 3;

  if (is_null(self->x.obsgeo)) {
    return NULL;
  }

  return get_double_array("obsgeo", self->x.obsgeo, 1, &size, (PyObject*)self);
}

static int
PyWcsprm_set_obsgeo(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp size = 3;

  if (is_null(self->x.obsgeo)) {
    return -1;
  }

  note_change(self);

  return set_double_array("obsgeo", value, 1, &size, self->x.obsgeo);
}

/*@null@*/ static PyObject*
PyWcsprm_get_pc(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {2, 2};

  if (is_null(self->x.pc)) {
    return NULL;
  }

  if ((self->x.altlin & has_pc) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No pc is present.");
    return NULL;
  }

  return get_double_array("pc", self->x.pc, 2, dims, (PyObject*)self);
}

static int
PyWcsprm_set_pc(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {2, 2};

  if (is_null(self->x.pc)) {
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x.altlin &= ~has_pc;
    return 0;
  }

  if (set_double_array("pc", value, 2, dims, self->x.pc)) {
    return -1;
  }

  self->x.altlin |= has_pc;

  note_change(self);

  return 0;
}

static PyObject*
PyWcsprm_get_phi0(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("phi0", self->x.cel.phi0);
}

static int
PyWcsprm_set_phi0(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_double("phi0", value, &(self->x.cel.phi0));
}

static PyObject*
PyWcsprm_get_piximg_matrix(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  const npy_intp dims[2] = {2, 2};

  if (is_null(self->x.lin.piximg)) {
    return NULL;
  }

  return get_double_array("piximg_matrix", self->x.lin.piximg, 2, dims,
                          (PyObject*)self);
}

static PyObject*
PyWcsprm_get_radesys(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.radesys)) {
    return NULL;
  }

  return get_string("radesys", self->x.radesys);
}

static int
PyWcsprm_set_radesys(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.radesys)) {
    return -1;
  }

  note_change(self);

  return set_string("radesys", value, self->x.radesys, 72);
}

static PyObject*
PyWcsprm_get_restfrq(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("restfrq", self->x.restfrq);
}

static int
PyWcsprm_set_restfrq(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.restfrq = UNDEFINED;
    return 0;
  }

  note_change(self);

  return set_double("restfrq", value, &self->x.restfrq);
}

static PyObject*
PyWcsprm_get_restwav(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("restwav", self->x.restwav);
}

static int
PyWcsprm_set_restwav(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.restwav = UNDEFINED;
    return 0;
  }

  note_change(self);

  return set_double("restwav", value, &self->x.restwav);
}

static PyObject*
PyWcsprm_get_spec(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("spec", self->x.spec);
}

/*@null@*/ static PyObject*
PyWcsprm_get_specsys(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.specsys)) {
    return NULL;
  }

  return get_string("specsys", self->x.specsys);
}

static int
PyWcsprm_set_specsys(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.specsys)) {
    return -1;
  }

  note_change(self);

  return set_string("specsys", value, self->x.specsys, 72);
}

/*@null@*/ static PyObject*
PyWcsprm_get_ssysobs(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssysobs)) {
    return NULL;
  }

  return get_string("ssysobs", self->x.ssysobs);
}

static int
PyWcsprm_set_ssysobs(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssysobs)) {
    return -1;
  }

  note_change(self);

  return set_string("ssysobs", value, self->x.ssysobs, 72);
}

/*@null@*/ static PyObject*
PyWcsprm_get_ssyssrc(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssyssrc)) {
    return NULL;
  }

  return get_string("ssyssrc", self->x.ssyssrc);
}

static int
PyWcsprm_set_ssyssrc(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssyssrc)) {
    return -1;
  }

  note_change(self);

  return set_string("ssyssrc", value, self->x.ssyssrc, 72);
}

static PyObject*
PyWcsprm_get_theta0(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("theta0", self->x.cel.theta0);
}

static int
PyWcsprm_set_theta0(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_double("theta0", value, &self->x.cel.theta0);
}

static PyObject*
PyWcsprm_get_velangl(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("velangl", self->x.velangl);
}

static int
PyWcsprm_set_velangl(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.velangl = UNDEFINED;
    return 0;
  }

  note_change(self);

  return set_double("velangl", value, &self->x.velangl);
}

static PyObject*
PyWcsprm_get_velosys(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("velosys", self->x.velosys);
}

static int
PyWcsprm_set_velosys(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.velosys = UNDEFINED;
    return 0;
  }

  note_change(self);

  return set_double("velosys", value, &self->x.velosys);
}

static PyObject*
PyWcsprm_get_zsource(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("zsource", self->x.zsource);
}

static int
PyWcsprm_set_zsource(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.zsource = UNDEFINED;
    return 0;
  }

  note_change(self);

  return set_double("zsource", value, &self->x.zsource);
}

/***************************************************************************
 * PyWcsprm definition structures
 */

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
  {"mix_fits", (PyCFunction)PyWcsprm_mix_fits, METH_VARARGS, doc_mix_fits},
  {"p2s", (PyCFunction)PyWcsprm_p2s, METH_O, doc_p2s},
  {"p2s_fits", (PyCFunction)PyWcsprm_p2s_fits, METH_O, doc_p2s_fits},
  {"print_contents", (PyCFunction)PyWcsprm_print_contents, METH_NOARGS, doc_print_contents},
  {"s2p", (PyCFunction)PyWcsprm_s2p, METH_O, doc_s2p},
  {"s2p_fits", (PyCFunction)PyWcsprm_s2p_fits, METH_O, doc_s2p_fits},
  {"set", (PyCFunction)PyWcsprm_set, METH_NOARGS, doc_set},
  {"set_ps", (PyCFunction)PyWcsprm_set_ps, METH_O, doc_set_ps},
  {"set_pv", (PyCFunction)PyWcsprm_set_pv, METH_O, doc_set_pv},
  {"spcfix", (PyCFunction)PyWcsprm_spcfix, METH_NOARGS, doc_spcfix},
  {"sptr", (PyCFunction)PyWcsprm_sptr, METH_NOARGS, doc_sptr},
  {"to_header", (PyCFunction)PyWcsprm_to_header, METH_VARARGS, doc_to_header},
  {"unitfix", (PyCFunction)PyWcsprm_unitfix, METH_VARARGS, doc_unitfix},
  {NULL}
};

PyTypeObject PyWcsprmType = {
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  "pywcs._Wcsprm",              /*tp_name*/
  sizeof(PyWcsprm),             /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PyWcsprm_dealloc, /*tp_dealloc*/
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
  doc_Wcsprm,                   /* tp_doc */
  0,                            /* tp_traverse */
  0,                            /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PyWcsprm_methods,             /* tp_methods */
  0,                            /* tp_members */
  PyWcsprm_getset,              /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc)PyWcsprm_init,      /* tp_init */
  0,                            /* tp_alloc */
  PyWcsprm_new,                 /* tp_new */
};

int
_setup_wcsprm_type(
    PyObject* m) {

  if (PyType_Ready(&PyWcsprmType) < 0) {
    return -1;
  }

  Py_INCREF(&PyWcsprmType);
  return PyModule_AddObject(m, "_Wcsprm", (PyObject *)&PyWcsprmType);
}
