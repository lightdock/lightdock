#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PyInt_AsUnsignedLongMask PyLong_AsUnsignedLongMask
#include <Python.h>
#include "structmember.h"
#include "numpy/arrayobject.h"


#define EPSILON 4.0
#define FACTOR 332.0
#define MAX_ES_CUTOFF 1.0
#define MIN_ES_CUTOFF -1.0
#define VDW_CUTOFF 1.0
#define HUGE_DISTANCE 10000.0
#define ELEC_DIST_CUTOFF 30.0
#define ELEC_DIST_CUTOFF2 ELEC_DIST_CUTOFF*ELEC_DIST_CUTOFF
#define VDW_DIST_CUTOFF 10.0
#define VDW_DIST_CUTOFF2 VDW_DIST_CUTOFF*VDW_DIST_CUTOFF


/**
 *
 * calculate_energy pyDock C implementation
 *
 **/
static PyObject * cdna_calculate_energy(PyObject *self, PyObject *args) {
    PyObject *receptor_coordinates, *ligand_coordinates = NULL;
    PyObject *tmp0, *tmp1 = NULL;
    PyArrayObject *rec_charges, *lig_charges, *rec_vdw, *lig_vdw, *rec_vdw_radii, *lig_vdw_radii = NULL;
    double atom_elec, total_elec, total_vdw, vdw_energy, vdw_radius, p6, k;
    unsigned int rec_len, lig_len, i, j;
    unsigned int interface_len, intf_array_size;
    unsigned int *interface_receptor = NULL, *interface_ligand = NULL;
    double **rec_array, **lig_array, x, y, z, distance2, interface_cutoff, interface_cutoff2;
    npy_intp dims[2];
    PyArray_Descr *descr;
    double *rec_c_charges, *lig_c_charges, *rec_c_vdw, *lig_c_vdw, *rec_c_vdw_radii, *lig_c_vdw_radii = NULL;
    PyObject *result = PyTuple_New(4);

    total_elec = 0.0;
    atom_elec = 0.0;
    total_vdw = 0.0;
    interface_len = 0;
    intf_array_size = 1;
    interface_cutoff = 3.9;

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

        // Get pointers to the Python array structures
        rec_c_charges = PyArray_GETPTR1(rec_charges, 0);
        lig_c_charges = PyArray_GETPTR1(lig_charges, 0);
        rec_c_vdw = PyArray_GETPTR1(rec_vdw, 0);
        lig_c_vdw = PyArray_GETPTR1(lig_vdw, 0);
        rec_c_vdw_radii = PyArray_GETPTR1(rec_vdw_radii, 0);
        lig_c_vdw_radii = PyArray_GETPTR1(lig_vdw_radii, 0);

        // Store interface
        interface_receptor = malloc(lig_len*sizeof(unsigned int));
        interface_ligand  = malloc(lig_len*sizeof(unsigned int));

        // For all atoms in receptor
        for (i = 0; i < rec_len; i++) {
            // For all atoms in ligand
            for (j = 0; j < lig_len; j++) {
                // Euclidean^2 distance
                x = rec_array[i][0] - lig_array[j][0];
                y = rec_array[i][1] - lig_array[j][1];
                z = rec_array[i][2] - lig_array[j][2];
                distance2 = x*x + y*y + z*z;

                // Electrostatics energy
                if (distance2 <= ELEC_DIST_CUTOFF2) {
                    atom_elec = (rec_c_charges[i] * lig_c_charges[j]) / distance2;
                    if (atom_elec >= (MAX_ES_CUTOFF*EPSILON/FACTOR)) atom_elec = MAX_ES_CUTOFF*EPSILON/FACTOR;
                    if (atom_elec <= (MIN_ES_CUTOFF*EPSILON/FACTOR)) atom_elec = MIN_ES_CUTOFF*EPSILON/FACTOR;
                    total_elec += atom_elec;
                }

                // Van der Waals energy
                if (distance2 <= VDW_DIST_CUTOFF2){
                    vdw_energy = sqrt(rec_c_vdw[i] * lig_c_vdw[j]);
                    vdw_radius = rec_c_vdw_radii[i] + lig_c_vdw_radii[j];
                    p6 = pow(vdw_radius, 6) / pow(distance2, 3);
                    k = vdw_energy * (p6*p6 - 2.0 * p6);
                    if (k > VDW_CUTOFF) k = VDW_CUTOFF;
                    total_vdw += k;
                }

                if (distance2 <= interface_cutoff2) {
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
        // Convert total electrostatics to Kcal/mol:
        //      - coordinates are in Ang
        //      - charges are in e (elementary charge units)
        total_elec = total_elec * FACTOR / EPSILON;

        // Free structures
        PyArray_Free(tmp0, rec_array);
        PyArray_Free(tmp1, lig_array);
        Py_DECREF(tmp0);
        Py_DECREF(tmp1);
    }
    
    interface_receptor = realloc(interface_receptor, interface_len*sizeof(unsigned int));
    interface_ligand = realloc(interface_ligand, interface_len*sizeof(unsigned int));
    dims[0] = interface_len;

    // Return a tuple with the following values for calculated energies:
    PyTuple_SetItem(result, 0, PyFloat_FromDouble(total_elec));
    PyTuple_SetItem(result, 1, PyFloat_FromDouble(total_vdw));
    PyTuple_SetItem(result, 2, PyArray_SimpleNewFromData(1, dims, NPY_UINT, interface_receptor));
    PyTuple_SetItem(result, 3, PyArray_SimpleNewFromData(1, dims, NPY_UINT, interface_ligand));
    return result;
}


/**
 *
 * Module methods table
 *
 **/
static PyMethodDef module_methods[] = {
    {"calculate_energy", (PyCFunction)cdna_calculate_energy, METH_VARARGS, "pyDockDNA C implementation"},
    {NULL}
};


/**
 *
 * Initialization function
 *
 **/
static struct PyModuleDef cdna =
{
    PyModuleDef_HEAD_INIT,
    "cdna",
    "",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_cdna(void) {
    import_array();
    return PyModule_Create(&cdna);
}

