#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PyInt_AsUnsignedLongMask PyLong_AsUnsignedLongMask
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

    *indexes_len = 0;
    n = 0;

    tmp0 = PyObject_GetAttrString(receptor_coordinates, "coordinates");
    tmp1 = PyObject_GetAttrString(ligand_coordinates, "coordinates");

    rec_len = PySequence_Size(tmp0);
    lig_len = PySequence_Size(tmp1);

    dims[1] = 3;
    dims[0] = rec_len;
    PyArray_AsCArray((PyObject **)&tmp0, (void **)&rec_array, dims, 2, PyArray_DescrFromType(NPY_DOUBLE));

    dims[0] = lig_len;
    PyArray_AsCArray((PyObject **)&tmp1, (void **)&lig_array, dims, 2, PyArray_DescrFromType(NPY_DOUBLE));

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
    PyObject *receptor, *ligand, *dfire_energy, *receptor_coordinates, *ligand_coordinates;
    PyObject *take, *intf_rec_array_object, *intf_lig_array_object, *array_object, *tmp0, *tmp1, **rec_objects, **lig_objects, *result = NULL;
    PyArrayObject *df_en_array;
    unsigned int n, m, i, j, d, dfire_bin, atoma, atomb, indexes_len, interface_len, *array, *interface_receptor, *interface_ligand, *indexes;
    double interface_cutoff, energy, *dfire_en_array;
    npy_intp dims[1];

    interface_cutoff = 3.9;
    energy = 0.;
    interface_len = 0;

    if (PyArg_ParseTuple(args, "OOOOO|d", &receptor, &ligand, &dfire_energy, &receptor_coordinates, &ligand_coordinates, &interface_cutoff)) {
        euclidean_dist(receptor_coordinates, ligand_coordinates, &indexes, &indexes_len);

        array = malloc(indexes_len*sizeof(unsigned int));
        interface_receptor = malloc(indexes_len*sizeof(unsigned int));
        interface_ligand = malloc(indexes_len*sizeof(unsigned int));

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

            if (d <= interface_cutoff) {
                interface_receptor[interface_len] = i;
                interface_ligand[interface_len++] = j;
            }

            atoma = PyInt_AsUnsignedLongMask(rec_objects[i]);
            atomb = PyInt_AsUnsignedLongMask(lig_objects[j]);

            dfire_bin = dist_to_bins[d] - 1;

            array[n] = atoma*168*20 + atomb*20 + dfire_bin;
        }

        dims[0] = indexes_len;
        tmp0 = PyImport_ImportModule("numpy");
        take = PyObject_GetAttrString(tmp0, "take");
        Py_DECREF(tmp0);
        array_object = PyArray_SimpleNewFromData(1, dims, NPY_UINT, array);
        df_en_array = (PyArrayObject *)PyObject_CallFunctionObjArgs(take, dfire_energy, array_object, NULL);
        dfire_en_array = PyArray_GETPTR1(df_en_array, 0);

        for (n = 0; n < dims[0]; n++) {
            energy += dfire_en_array[n];
        }
        
        free(array);
        free(indexes);

        Py_DECREF(df_en_array);
        Py_DECREF(array_object);
        Py_DECREF(take);

    }

    dims[0] = interface_len;

    interface_receptor = realloc(interface_receptor, interface_len*sizeof(unsigned int));
    interface_ligand = realloc(interface_ligand, interface_len*sizeof(unsigned int));

    result = PyTuple_New(3);
    PyTuple_SET_ITEM(result, 0, PyFloat_FromDouble((energy*0.0157 - 4.7)*-1));
    PyTuple_SET_ITEM(result, 1, PyArray_SimpleNewFromData(1, dims, NPY_UINT, interface_receptor));
    PyTuple_SET_ITEM(result, 2, PyArray_SimpleNewFromData(1, dims, NPY_UINT, interface_ligand));

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
static struct PyModuleDef cdfire =
{
    PyModuleDef_HEAD_INIT,
    "cdfire",
    "",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_cdfire(void) {
    import_array();
    return PyModule_Create(&cdfire);
}
