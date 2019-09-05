#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PyInt_AsUnsignedLongMask PyLong_AsUnsignedLongMask
#include <Python.h>
#include "structmember.h"
#include "numpy/arrayobject.h"


#define EPSILON 4.0
#define FACTOR 332.0
#define CUTON 7.0
#define CUTOFF 9.0
#define CUTON2 CUTON*CUTON
#define CUTOFF2 CUTOFF*CUTOFF
#define VDW_CUTOFF 5000000.0


/**
 *
 * calculate_energy C implementation
 *
 **/
static PyObject * sd_calculate_energy(PyObject *self, PyObject *args) {
    PyObject *receptor_coordinates, *ligand_coordinates = NULL;
    PyObject *tmp0, *tmp1 = NULL;
    PyArrayObject *rec_charges, *lig_charges, *rec_vdw, *lig_vdw, *rec_vdw_radii, *lig_vdw_radii = NULL;
    double energy, atom_elec, total_elec, atom_vdw, total_vdw, vdw_energy,vdw_radius, p6, k;
    unsigned int rec_len, lig_len, i, j, interface_len, intf_array_size;
    unsigned int *interface_receptor = NULL, *interface_ligand = NULL;
    double **rec_array, **lig_array, x, y, z, distance, interface_cutoff, interface_cutoff2;
    npy_intp dims[2];
    PyArray_Descr *descr;
    double *rec_c_charges, *lig_c_charges, *rec_c_vdw, *lig_c_vdw, *rec_c_vdw_radii, *lig_c_vdw_radii = NULL;
    PyObject *result = PyTuple_New(3);

    energy = 0.;
    interface_cutoff = 3.9;
    interface_len = 0;
    intf_array_size = 1;

    if (PyArg_ParseTuple(args, "OOOOOOOO|d",
            &receptor_coordinates, &ligand_coordinates, &rec_charges, &lig_charges,
            &rec_vdw, &lig_vdw, &rec_vdw_radii, &lig_vdw_radii, &interface_cutoff)) {

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

        total_elec = 0.0;
        atom_elec = 0.0;
        total_vdw = 0.0;

        // Get pointers to the Python array structures
        rec_c_charges = PyArray_GETPTR1(rec_charges, 0);
        lig_c_charges = PyArray_GETPTR1(lig_charges, 0);
        rec_c_vdw = PyArray_GETPTR1(rec_vdw, 0);
        lig_c_vdw = PyArray_GETPTR1(lig_vdw, 0);
        rec_c_vdw_radii = PyArray_GETPTR1(rec_vdw_radii, 0);
        lig_c_vdw_radii = PyArray_GETPTR1(lig_vdw_radii, 0);

        interface_receptor = malloc(lig_len*sizeof(unsigned int));
        interface_ligand  = malloc(lig_len*sizeof(unsigned int));

        for (i = 0; i < rec_len; i++) {
            atom_vdw = 0.0;
            for (j = 0; j < lig_len; j++) {
                x = rec_array[i][0] - lig_array[j][0];
                y = rec_array[i][1] - lig_array[j][1];
                z = rec_array[i][2] - lig_array[j][2];
                distance = x*x + y*y + z*z;
                if (distance < CUTOFF2)
                {
                    // Electrostatics
                    atom_elec = (rec_c_charges[i] * lig_c_charges[j]) / distance;
                    // Convert total electrostatics to:
                    // Transform to Kcal/mol:
                    //      - coordinates are in Ang
                    //      - charges are in e (elementary charge units)
                    atom_elec *= FACTOR/EPSILON;

                    // VdW
                    vdw_energy = sqrt(rec_c_vdw[i] * lig_c_vdw[j]);
                    vdw_radius = rec_c_vdw_radii[i] + lig_c_vdw_radii[j];
                    p6 = pow(vdw_radius, 6) / pow(distance, 3);
                    k = vdw_energy * (p6*p6 - 2.0 * p6);
                    atom_vdw += k;
                    if (atom_vdw > VDW_CUTOFF) atom_vdw = VDW_CUTOFF;

                    if (distance < CUTON2)
                    {
                        energy += atom_elec + atom_vdw;
                    } else {
                        energy += (atom_elec + atom_vdw) * ( (CUTOFF2 - distance)*(CUTOFF2 - distance) *
                                    (CUTOFF2 + 2.*distance - 3.0*CUTON2) / ((CUTOFF2-CUTON2)*(CUTOFF2-CUTON2)*(CUTOFF2-CUTON2)) );
                    }
                }

                if (distance <= interface_cutoff2) {
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

        // Free structures
        Py_DECREF(rec_c_charges);
        Py_DECREF(lig_c_charges);
        Py_DECREF(rec_c_vdw);
        Py_DECREF(lig_c_vdw);
        Py_DECREF(rec_c_vdw_radii);
        Py_DECREF(lig_c_vdw_radii);
        Py_DECREF(descr);
        PyArray_Free(tmp0, rec_array);
        PyArray_Free(tmp1, lig_array);
    }

    interface_receptor = realloc(interface_receptor, interface_len*sizeof(unsigned int));
    interface_ligand = realloc(interface_ligand, interface_len*sizeof(unsigned int));
    dims[0] = interface_len;

    PyTuple_SetItem(result, 0, PyFloat_FromDouble(energy * -1.));
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
    {"calculate_energy", (PyCFunction)sd_calculate_energy, METH_VARARGS, "SD energy C implementation"},
    {NULL}
};


/**
 *
 * Initialization function
 *
 **/
static struct PyModuleDef sd =
{
    PyModuleDef_HEAD_INIT,
    "sd",
    "",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_sd(void) {
    import_array();
    return PyModule_Create(&sd);
}
