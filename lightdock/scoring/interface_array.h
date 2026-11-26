/**
 *
 * Helper shared by the C scoring functions
 *
 **/
#ifndef LIGHTDOCK_INTERFACE_ARRAY_H
#define LIGHTDOCK_INTERFACE_ARRAY_H

#include <stdlib.h>
#include <string.h>


/**
 * Moves a malloc'ed buffer of interface indexes into a new NumPy array.
 *
 * The returned array owns its memory and the buffer is always freed, so it must
 * not be used after this call. Returns NULL (with a Python exception set) if
 * the array cannot be created.
 **/
static PyObject * interface_to_array(unsigned int *buffer, npy_intp length) {
    npy_intp dims[1];
    PyObject *array;

    dims[0] = length;
    array = PyArray_SimpleNew(1, dims, NPY_UINT);
    if (array != NULL && length > 0) {
        memcpy(PyArray_DATA((PyArrayObject *)array), buffer, length*sizeof(unsigned int));
    }
    free(buffer);

    return array;
}

#endif
