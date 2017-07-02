#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
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
    unsigned int rec_len, lig_len, i, j;
    double **rec_array, **lig_array, x, y, z, distance;
    npy_intp dims[2];
    PyArray_Descr *descr;
    double *rec_c_charges, *lig_c_charges, *rec_c_vdw, *lig_c_vdw, *rec_c_vdw_radii, *lig_c_vdw_radii = NULL;

    energy = 0.;

    if (PyArg_ParseTuple(args, "OOOOOOOO",
            &receptor_coordinates, &ligand_coordinates, &rec_charges, &lig_charges,
            &rec_vdw, &lig_vdw, &rec_vdw_radii, &lig_vdw_radii)) {

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

    return PyFloat_FromDouble(energy * -1.);
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
PyMODINIT_FUNC initsd(void) {

    Py_InitModule3("sd", module_methods, "sd object");
    import_array();
}
