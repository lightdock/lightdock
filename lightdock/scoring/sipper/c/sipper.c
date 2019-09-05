#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PyInt_AsUnsignedLongMask PyLong_AsUnsignedLongMask
#include <Python.h>
#include <numpy/arrayobject.h>
#include "structmember.h"

#define DISTANCE2_CUTOFF 25.0


/**
 *
 * calculate_sipper C implementation
 *
 **/
static PyObject * calculate_sipper(PyObject *self, PyObject *args) {
    PyObject *receptor_coordinates, *ligand_coordinates = NULL;
    PyObject *tmp0, *tmp1 = NULL;
    PyArrayObject *sipper_energy, *receptor_indexes, *ligand_indexes, *rec_res_atoms, *lig_res_atoms = NULL;
    PyArrayObject *receptor_oda, *ligand_oda = NULL;
    double total_sipper, total_oda, total_energy, interface_cutoff, interface_cutoff2;;
    unsigned int rec_len, lig_len, i, j, rec_res_len, lig_res_len, atom_i, atom_j;
    unsigned int interface_len, intf_array_size, *interface_receptor = NULL, *interface_ligand = NULL;
    double **rec_array, **lig_array, x, y, z;
    npy_intp dims[2];
    PyArray_Descr *descr;
    int *receptor_c_indexes, *ligand_c_indexes, *rec_c_res_atoms, *lig_c_res_atoms = NULL;
    double *receptor_c_oda, *ligand_c_oda = NULL;
    unsigned int receptor_init, ligand_init = 0;
    PyObject *result = PyTuple_New(3);

    total_sipper = 0.0;
    total_oda = 0.0;
    total_energy = 0.0;
    interface_cutoff = 3.9;
    interface_len = 0;
    intf_array_size = 1;

    if (PyArg_ParseTuple(args, "OOOOOOOiiOO|d",
            &receptor_coordinates, &ligand_coordinates, &sipper_energy, &receptor_indexes, &ligand_indexes,
            &rec_res_atoms, &lig_res_atoms, &rec_res_len, &lig_res_len, &receptor_oda, &ligand_oda, &interface_cutoff)) {

        interface_cutoff2 = interface_cutoff*interface_cutoff;

        descr = PyArray_DescrFromType(NPY_DOUBLE);

        tmp0 = PyObject_GetAttrString(receptor_coordinates, "coordinates");
        tmp1 = PyObject_GetAttrString(ligand_coordinates, "coordinates");

        rec_len = PySequence_Size(tmp0);
        lig_len = PySequence_Size(tmp1);

        dims[1] = 3;
        dims[0] = rec_len;
        PyArray_AsCArray((PyObject **)&tmp0, (void **)&rec_array, dims, 2, descr);

        dims[0] = lig_len;
        PyArray_AsCArray((PyObject **)&tmp1, (void **)&lig_array, dims, 2, descr);

        // Get pointers to the Python array structures
        receptor_c_indexes = PyArray_GETPTR1(receptor_indexes, 0);
        ligand_c_indexes = PyArray_GETPTR1(ligand_indexes, 0);
        rec_c_res_atoms = PyArray_GETPTR1(rec_res_atoms, 0);
        lig_c_res_atoms = PyArray_GETPTR1(lig_res_atoms, 0);
        receptor_c_oda = PyArray_GETPTR1(receptor_oda, 0);
        ligand_c_oda = PyArray_GETPTR1(ligand_oda, 0);

        interface_receptor = malloc(lig_len*sizeof(unsigned int));
        interface_ligand  = malloc(lig_len*sizeof(unsigned int));

        // For all residues in receptor
        receptor_init = 0;
        for (i = 0; i < rec_res_len; i++) {
            // For all residues in ligand
            ligand_init = 0;
            for (j = 0; j < lig_res_len; j++) {
                for(atom_i=receptor_init; atom_i<(receptor_init+rec_c_res_atoms[i]); atom_i++) {
                    for(atom_j=ligand_init; atom_j<(ligand_init+lig_c_res_atoms[j]); atom_j++) {
                        x = rec_array[atom_i][0] - lig_array[atom_j][0];
                        x *= x;
                        if (x > DISTANCE2_CUTOFF) continue;
                        y = rec_array[atom_i][1] - lig_array[atom_j][1];
                        y *= y;
                        if (y > DISTANCE2_CUTOFF) continue;
                        z = rec_array[atom_i][2] - lig_array[atom_j][2];
                        z *= z;
                        if (z > DISTANCE2_CUTOFF) continue;

                        if (DISTANCE2_CUTOFF > (x+y+z)) {
                            total_sipper += *((double *) PyArray_GETPTR2(sipper_energy,
                                                                        receptor_c_indexes[i], ligand_c_indexes[j]));
                            total_oda += receptor_c_oda[i] + ligand_c_oda[j];
                            break;
                        }

                        if ((x+y+z) <= interface_cutoff2) {
                           interface_receptor[interface_len] = i;
                           interface_ligand[interface_len++] = j;
                        }
                    }

                    if (((interface_len + lig_len - 1)/lig_len + 1) > intf_array_size) {
                        intf_array_size++;
                        interface_receptor = realloc(interface_receptor, intf_array_size*lig_len*sizeof(unsigned int));
                        interface_ligand = realloc(interface_ligand, intf_array_size*lig_len*sizeof(unsigned int));
                    }
                }
                ligand_init += lig_c_res_atoms[j];
            }
            receptor_init += rec_c_res_atoms[i];
        }

        // Free structures
        PyArray_Free(tmp0, rec_array);
        PyArray_Free(tmp1, lig_array);
    }
    //printf("sipper=%5.3f oda=%5.3f\n", total_sipper, total_oda);
    total_energy = -1.0*(total_sipper - 0.019 * total_oda);

    interface_receptor = realloc(interface_receptor, interface_len*sizeof(unsigned int));
    interface_ligand = realloc(interface_ligand, interface_len*sizeof(unsigned int));
    dims[0] = interface_len;

    // Return a tuple with the following values for calculated energies:
    PyTuple_SetItem(result, 0, PyFloat_FromDouble(total_energy));
    PyTuple_SetItem(result, 1, PyArray_SimpleNewFromData(1, dims, NPY_UINT, interface_receptor));
    PyTuple_SetItem(result, 2, PyArray_SimpleNewFromData(1, dims, NPY_UINT, interface_ligand));
    return result;
}


/**
 *
 * Module methods table
 *
 **/
static PyMethodDef module_methods[] = {
    {"calculate_sipper", (PyCFunction)calculate_sipper, METH_VARARGS, "calculate_sipper C implementation"},
    {NULL}
};


/**
 *
 * Initialization function
 *
 **/
static struct PyModuleDef sipper =
{
    PyModuleDef_HEAD_INIT,
    "sipper",
    "",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_sipper(void) {
    import_array();
    return PyModule_Create(&sipper);
}

