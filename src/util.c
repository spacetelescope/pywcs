#include "util.h"
#include "isnan.h"

PyObject*
PyArrayProxy_New(PyObject* self, int nd, const npy_intp* dims,
                 int typenum, const void* data) {
  PyArray_Descr* type_descr = NULL;
  PyObject*      result     = NULL;

  type_descr = (PyArray_Descr*)PyArray_DescrFromType(typenum);
  if (type_descr == NULL)
    return NULL;
  result = (PyObject*)PyArray_NewFromDescr(&PyArray_Type, type_descr,
                                           nd, (npy_intp*)dims, NULL, (void*)data,
                                           NPY_CONTIGUOUS | NPY_WRITEABLE,
                                           NULL);

  if (result == NULL)
    return NULL;
  Py_INCREF(self);
  PyArray_BASE(result) = (PyObject*)self;
  return result;
}

void
offset_array(PyArrayObject* array, double value) {
  int     size = 1;
  int     i    = 0;
  double* data = NULL;

  for (i = 0; i < PyArray_NDIM(array); ++i)
    size *= PyArray_DIM(array, i);

  data = (double*)PyArray_DATA(array);

  offset_c_array(data, size, value);
}

void
copy_array_to_c_double(PyArrayObject* array, double* dest) {
  int     size = 1;
  int     i    = 0;
  double* data = NULL;

  data = (double*)PyArray_DATA(array);

  for (i = 0; i < PyArray_NDIM(array); ++i)
    size *= PyArray_DIM(array, i);

  for (i = 0; i < size; ++i, ++dest, ++data) {
    if (isnan64(*data))
      *dest = UNDEFINED;
    else
      *dest = *data;
  }
}

void
copy_array_to_c_int(PyArrayObject* array, int* dest) {
  int  size = 1;
  int  i    = 0;
  int* data = NULL;

  data = (int*)PyArray_DATA(array);

  for (i = 0; i < PyArray_NDIM(array); ++i)
    size *= PyArray_DIM(array, i);

  for (i = 0; i < size; ++i, ++dest, ++data)
    *dest = *data;
}
