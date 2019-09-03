#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PyInt_AsUnsignedLongMask PyLong_AsUnsignedLongMask
#include <Python.h>
#include <numpy/arrayobject.h>
#include "structmember.h"


/**
 *
 * calculate_dfire2 C implementation
 *
 **/
static PyObject * cdfire2_calculate_dfire2(PyObject *self, PyObject *args) {
    PyObject *res_index, *atom_index, *coordinates, *dfire2_energy, *result = NULL;
    unsigned int mol_length, i, j, b, atom1, atom2, index, interface_len, *interface_receptor = NULL, *interface_ligand = NULL;
    double energy, dist, **mol_array, *energies, interface_cutoff;
    int *atom_indexes, *res_indexes;

    npy_intp dims[2];

    PyObject *res_array;
    PyObject *atom_array;
    PyObject *energy_array;

    interface_cutoff = 3.9;
    energy = 0.;
    interface_len = 0;

    if (PyArg_ParseTuple(args, "OOOOi|d", &res_index, &atom_index, &coordinates, &dfire2_energy, &mol_length, &interface_cutoff)) {
        PyArray_Descr *descr;

        // Coordinates to C Array
        descr = PyArray_DescrFromType(NPY_DOUBLE);
        dims[1] = 3;
        dims[0] = mol_length;
        PyArray_AsCArray((PyObject **)&coordinates, (void **)&mol_array, dims, 2, descr);

        // Residue indexes
        res_array = PyArray_FROM_OTF(res_index, NPY_INT32, NPY_ARRAY_IN_ARRAY);
        res_indexes = (int*)PyArray_DATA(res_array);

        // Atom type indexes
        atom_array = PyArray_FROM_OTF(atom_index, NPY_INT32, NPY_ARRAY_IN_ARRAY);
        atom_indexes = (int*)PyArray_DATA(atom_array);

        // DFIRE2 potentials
        energy_array = PyArray_FROM_OTF(dfire2_energy, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        //int N = (int)PyArray_DIM(energy_array, 0);
        //printf("%d", N);
        energies = (double*)PyArray_DATA(energy_array);

        interface_receptor = malloc(mol_length*sizeof(unsigned int));
        interface_ligand = malloc(mol_length*sizeof(unsigned int));

        // Calculate euclidean distance
        for (i = 0; i < mol_length; i++) {
            for (j = i+1; j < mol_length; j++) {
                if (res_indexes[i] != res_indexes[j]) {
                    // Euclidean distance * 2
                    dist = sqrt(pow((mol_array[i][0] - mol_array[j][0]), 2.0) +
                                pow((mol_array[i][1] - mol_array[j][1]), 2.0) +
                                pow((mol_array[i][2] - mol_array[j][2]), 2.0))*2;
                    if (dist <= interface_cutoff) {
                        interface_receptor[interface_len] = i;
                        interface_ligand[interface_len++] = j;
                    }
                    // Distance to bin
                    b = (int)dist;
                    if (b < 30) {
                        atom1 = atom_indexes[i];
                        atom2 = atom_indexes[j];
                        // 1D energies array, original shape was (167, 167, 30)
                        index = atom1*167*30 + atom2*30 + b;
                        //printf("%d %d %f %d %f\n", res_indexes[i], res_indexes[j], dist, b, energies[index]);
                        energy += energies[index];
                    }
                }
            }
        }

        // Free memory
        Py_DECREF(res_array);
        Py_DECREF(atom_array);
        Py_DECREF(energy_array);
        PyArray_Free(coordinates, mol_array);
    }

    interface_receptor = realloc(interface_receptor, interface_len*sizeof(unsigned int));
    interface_ligand = realloc(interface_ligand, interface_len*sizeof(unsigned int));

    result = PyTuple_New(3);
    PyTuple_SET_ITEM(result, 0, PyFloat_FromDouble(energy/100.));
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
    {"calculate_dfire2", (PyCFunction)cdfire2_calculate_dfire2, METH_VARARGS, "calculate_dfire2 C implementation"},
    {NULL}
};


/**
 *
 * Initialization function
 *
 **/
static struct PyModuleDef cdfire2 =
{
    PyModuleDef_HEAD_INIT,
    "cdfire2",
    "",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_cdfire2(void) {
    import_array();
    return PyModule_Create(&cdfire2);
}

