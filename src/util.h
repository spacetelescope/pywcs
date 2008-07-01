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

#ifndef __UTIL_H__
#define __UTIL_H__

#include <Python.h>
#include <numpy/arrayobject.h>
#include <wcsmath.h>

#include "isnan.h"

PyObject*
PyArrayProxy_New(PyObject* self, int nd, const npy_intp* dims,
                 int typenum, const void* data);

PyObject *
PyWcsprmListProxy_New(PyObject* owner, Py_ssize_t size, char (*array)[72]);

static inline void
offset_c_array(double* value, size_t size, double offset) {
  double* end = value + size;

  for ( ; value != end; ++value)
    *value += offset;
}

void
offset_array(PyArrayObject* array, double value);

void
copy_array_to_c_double(PyArrayObject* array, double* dest);

void
copy_array_to_c_int(PyArrayObject* array, int* dest);

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


#endif /* __UTIL_H__ */
