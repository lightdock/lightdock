#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include "structmember.h"
#include "numpy/arrayobject.h"


/**
 *
 * DFIRE distances
 *
 **/
static unsigned int dist_to_bins[50] = {
         1,  1,  1,  2,  3,  4,  5,  6,  7,  8,
         9, 10, 11, 12, 13, 14, 14, 15, 15, 16,
        16, 17, 17, 18, 18, 19, 19, 20, 20, 21,
        21, 22, 22, 23, 23, 24, 24, 25, 25, 26,
        26, 27, 27, 28, 28, 29, 29, 30, 30, 31};


/**
 *
 * Auxiliary function
 *
 **/
int compare(const void *a, const void *b) {
    const unsigned int *da = (const unsigned int *) a;
    const unsigned int *db = (const unsigned int *) b;

    return (*da > *db) - (*da < *db);
}


/**
 *
 * Computation of Euclidean distances and selection of nearest atoms
 *
 **/
void euclidean_dist(PyObject *receptor_coordinates, PyObject *ligand_coordinates,
                    unsigned int **indexes, unsigned int *indexes_len) {
    PyObject *tmp0, *tmp1;
    unsigned int rec_len, lig_len, i, j, n;
    double dist, **rec_array, **lig_array;
    npy_intp dims[2];
    PyArray_Descr *descr;

    descr = PyArray_DescrFromType(NPY_DOUBLE);

    *indexes_len = 0;
    n = 0;

    tmp0 = PyObject_GetAttrString(receptor_coordinates, "coordinates");
    tmp1 = PyObject_GetAttrString(ligand_coordinates, "coordinates");

    rec_len = PySequence_Size(tmp0);
    lig_len = PySequence_Size(tmp1);

    dims[1] = 3;
    dims[0] = rec_len;
    PyArray_AsCArray((PyObject **)&tmp0, (void **)&rec_array, dims, 2, descr);

    dims[0] = lig_len;
    PyArray_AsCArray((PyObject **)&tmp1, (void **)&lig_array, dims, 2, descr);

    *indexes = malloc(3*rec_len*lig_len*sizeof(unsigned int));

    for (i = 0; i < rec_len; i++) {
        for (j = 0; j < lig_len; j++) {
            dist = pow((rec_array[i][0] - lig_array[j][0]), 2.0) +
                   pow((rec_array[i][1] - lig_array[j][1]), 2.0) +
                   pow((rec_array[i][2] - lig_array[j][2]), 2.0);
            if (dist <= 225.) {
                (*indexes)[n++] = i;
                (*indexes)[n++] = j;
                (*indexes)[n++] = (sqrt(dist)*2.0 - 1.0);
                (*indexes_len)++;
            }
        }
    }
    *indexes = realloc(*indexes, n*sizeof(unsigned int));
    PyArray_Free(tmp0, rec_array);
    PyArray_Free(tmp1, lig_array);
    Py_DECREF(tmp0);
    Py_DECREF(tmp1);
}


/**
 *
 * calculate_dfire C implementation
 *
 **/
static PyObject * cdfire_calculate_dfire(PyObject *self, PyObject *args) {
    PyObject *receptor, *ligand, *dfire_energy, *receptor_coordinates, *ligand_coordinates = NULL;
    PyObject *take, *array_tuple, *array_object, *tmp0, *tmp1, **rec_objects, **lig_objects, *result = NULL;
    PyArray_Descr *descr;
    PyArrayObject *df_en_array;
    unsigned int n, m, i, j, d, dfire_bin, atoma, atomb, indexes_len, *array, *indexes;
    double energy, *dfire_en_array;
    npy_intp dims[1];
    unsigned long elapsed_dist, elapsed_indexes, elapsed_read;

    energy = 0.;
    elapsed_dist = elapsed_indexes = elapsed_read = 0;

    if (PyArg_ParseTuple(args, "OOOOO", &receptor, &ligand, &dfire_energy, &receptor_coordinates, &ligand_coordinates)) {
        euclidean_dist(receptor_coordinates, ligand_coordinates, &indexes, &indexes_len);
        array = malloc(indexes_len*sizeof(unsigned int));

        // Do not need to free rec_objects and lig_objects
        tmp0 = PyObject_GetAttrString(receptor, "objects");
        tmp1 = PySequence_Fast(tmp0, "");
        Py_DECREF(tmp0);
        rec_objects = PySequence_Fast_ITEMS(tmp1);
        Py_DECREF(tmp1);

        tmp0 = PyObject_GetAttrString(ligand, "objects");
        tmp1 = PySequence_Fast(tmp0, "");
        Py_DECREF(tmp0);
        lig_objects = PySequence_Fast_ITEMS(tmp1);
        Py_DECREF(tmp1);

        for (n = m = 0; n < indexes_len; n++) {

            i = indexes[m++];
            j = indexes[m++];
            d = indexes[m++];

            atoma = PyInt_AsUnsignedLongMask(rec_objects[i]);
            atomb = PyInt_AsUnsignedLongMask(lig_objects[j]);

            dfire_bin = dist_to_bins[d] - 1;

            array[n] = atoma*167*20 + atomb*20 + dfire_bin;
        }
        descr = PyArray_DescrFromType(NPY_DOUBLE);
        dims[0] = indexes_len;
        tmp0 = PyImport_ImportModule("numpy");
        take = PyObject_GetAttrString(tmp0, "take");
        Py_DECREF(tmp0);
        array_object = PyArray_SimpleNewFromData(1, dims, NPY_UINT, array);
        array_tuple = PyTuple_New(2);
        Py_INCREF(dfire_energy);
        PyTuple_SET_ITEM(array_tuple, 0, dfire_energy);
        PyTuple_SET_ITEM(array_tuple, 1, array_object);
        df_en_array = (PyArrayObject *)PyObject_Call(take, array_tuple, NULL);
        Py_DECREF(array_tuple);
        dfire_en_array = PyArray_GETPTR1(df_en_array, 0);

        for (n = 0; n < dims[0]; n++) {
            energy += dfire_en_array[n];
        }
        
        free(array);
        Py_DECREF(df_en_array);
        free(indexes);

        Py_DECREF(take);

        result = PyFloat_FromDouble((energy*0.0157 - 4.7)*-1);
    }
    return result;
}


/**
 *
 * Module methods table
 *
 **/
static PyMethodDef module_methods[] = {
    {"calculate_dfire", (PyCFunction)cdfire_calculate_dfire, METH_VARARGS, "calculate_dfire C implementation"},
    {NULL}
};


/**
 *
 * Initialization function
 *
 **/
PyMODINIT_FUNC initcdfire(void) {

    Py_InitModule3("cdfire", module_methods, "cdfire object");
    import_array();
}

