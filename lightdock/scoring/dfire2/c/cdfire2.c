#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PyInt_AsUnsignedLongMask PyLong_AsUnsignedLongMask
#include <Python.h>
#include <numpy/arrayobject.h>
#include "structmember.h"
#include "../../interface_array.h"


/**
 *
 * calculate_dfire2 C implementation
 *
 **/
static PyObject * cdfire2_calculate_dfire2(PyObject *self, PyObject *args) {
    PyObject *res_index, *atom_index, *coordinates, *dfire2_energy, *result = NULL;
    unsigned int mol_length, rec_length, i, j, b, atom1, atom2, index, interface_len, interface_size, *interface_receptor = NULL, *interface_ligand = NULL;
    int is_intermolecular, is_scored;
    double energy, dist, **mol_array, *energies, interface_cutoff;
    int *atom_indexes, *res_indexes;

    npy_intp dims[2];

    PyObject *res_array;
    PyObject *atom_array;
    PyObject *energy_array;

    interface_cutoff = 3.9;
    energy = 0.;
    interface_len = 0;

    if (PyArg_ParseTuple(args, "OOOOii|d", &res_index, &atom_index, &coordinates, &dfire2_energy, &mol_length, &rec_length,
                         &interface_cutoff)) {
        PyArray_Descr *descr;

        // Coordinates to C Array
        descr = PyArray_DescrFromType(NPY_DOUBLE);
        dims[1] = 3;
        dims[0] = mol_length;
        PyArray_AsCArray((PyObject **)&coordinates, (void **)&mol_array, dims, 2, descr);

        // Residue indexes
        res_array = PyArray_FROM_OTF(res_index, NPY_INT32, NPY_ARRAY_IN_ARRAY);
        res_indexes = (int*)PyArray_DATA((PyArrayObject*)res_array);

        // Atom type indexes
        atom_array = PyArray_FROM_OTF(atom_index, NPY_INT32, NPY_ARRAY_IN_ARRAY);
        atom_indexes = (int*)PyArray_DATA((PyArrayObject*)atom_array);

        // DFIRE2 potentials
        energy_array = PyArray_FROM_OTF(dfire2_energy, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        //int N = (int)PyArray_DIM(energy_array, 0);
        //printf("%d", N);
        energies = (double*)PyArray_DATA((PyArrayObject*)energy_array);

        // Initial size of the interface buffers, they grow if required
        interface_size = mol_length + 1;
        interface_receptor = malloc(interface_size*sizeof(unsigned int));
        interface_ligand = malloc(interface_size*sizeof(unsigned int));

        // Calculate euclidean distance. Receptor atoms come first, followed by ligand atoms
        for (i = 0; i < mol_length; i++) {
            for (j = i+1; j < mol_length; j++) {
                is_intermolecular = (i < rec_length && j >= rec_length);
                is_scored = (res_indexes[i] != res_indexes[j]);
                if (!is_intermolecular && !is_scored) continue;

                dist = sqrt(pow((mol_array[i][0] - mol_array[j][0]), 2.0) +
                            pow((mol_array[i][1] - mol_array[j][1]), 2.0) +
                            pow((mol_array[i][2] - mol_array[j][2]), 2.0));

                // Interface: receptor-ligand contacts, using each molecule own indexes
                if (is_intermolecular && dist <= interface_cutoff) {
                    if (interface_len == interface_size) {
                        interface_size *= 2;
                        interface_receptor = realloc(interface_receptor, interface_size*sizeof(unsigned int));
                        interface_ligand = realloc(interface_ligand, interface_size*sizeof(unsigned int));
                    }
                    interface_receptor[interface_len] = i;
                    interface_ligand[interface_len++] = j - rec_length;
                }

                if (is_scored) {
                    // Distance to bin, which is Euclidean distance * 2
                    b = (int)(dist*2);
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

    result = PyTuple_New(3);
    PyTuple_SET_ITEM(result, 0, PyFloat_FromDouble(energy/100.));
    PyTuple_SET_ITEM(result, 1, interface_to_array(interface_receptor, interface_len));
    PyTuple_SET_ITEM(result, 2, interface_to_array(interface_ligand, interface_len));

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

