#include <Python.h>
#include "structmember.h"
#include <math.h>


/**
 *
 * Exception object
 *
 **/
static PyObject *QuaternionError;


/**
 *
 * Quaternion object
 *
 **/
static PyTypeObject QuaternionType;

typedef struct {
    PyObject_HEAD
    double w;
    double x;
    double y;
    double z;
} Quaternion;


/**
 *
 * tp_dealloc
 *
 **/
static void Quaternion_dealloc(Quaternion* self) {
    self->ob_type->tp_free((PyObject*)self);
}


/**
 *
 * __init__ implementation
 *
 **/
static int Quaternion_init(Quaternion *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = {"w", "x", "y", "z", NULL};

    self->w = 1.;
    self->x = 0.;
    self->y = 0.;
    self->z = 0.;

    if (PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist, &self->w, &self->x, &self->y, &self->z)) {
        return 0;
    }

    return -1;
}


/**
 *
 * Auxiliary functions
 *
 **/
Quaternion * Quaternion_create(double w, double x, double y, double z) {
    Quaternion *new;
    
    new = PyObject_New(Quaternion, &QuaternionType);
    new->w = w;
    new->x = x;
    new->y = y;
    new->z = z;

    return new;
}

unsigned int compare(double a, double b) {
    return fabs(a-b) < 1e-7;
}

double dot(double w, double x, double y, double z, double w2, double x2, double y2, double z2) {
    return w*w2 + x*x2 + y*y2 + z*z2;
}

double norm2(double w, double x, double y, double z) {
    return w*w + x*x + y*y + z*z;
}

double norm(double w, double x, double y, double z) {
    return sqrt(norm2(w, x, y, z));
}


/**
 *
 * Quaternion type members
 *
 **/
static PyMemberDef Quaternion_members[] = {
    {"w", T_DOUBLE, offsetof(Quaternion, w), 0, "w"},
    {"x", T_DOUBLE, offsetof(Quaternion, x), 0, "x"},
    {"y", T_DOUBLE, offsetof(Quaternion, y), 0, "y"},
    {"z", T_DOUBLE, offsetof(Quaternion, z), 0, "z"},
    {NULL}  /* Sentinel */
};


/**
 *
 * __eq__ and __ne__ implementation
 *
 **/
static PyObject * Quaternion_richcompare(Quaternion *self, PyObject *operand, int op) {
    PyObject *result = NULL;
    Quaternion *other;

    if (PyObject_TypeCheck(operand, &QuaternionType)) {
        other = (Quaternion *)operand;
        if (op == Py_EQ) {
            result = (compare(self->w, other->w) && compare(self->x, other->x) && compare(self->y, other->y) && compare(self->z, other->z)) ? Py_True : Py_False;
        } else if (op == Py_NE) {
            result = (compare(self->w, other->w) && compare(self->x, other->x) && compare(self->y, other->y) && compare(self->z, other->z)) ? Py_False : Py_True;
        }

        return result; 
    }
    PyErr_SetString(PyExc_TypeError, "Expecting Quaternion type");

    return NULL;
}


/**
 *
 * clone implementation
 *
 **/
static PyObject * Quaternion_clone(Quaternion *self) {
    return (PyObject *)Quaternion_create(self->w, self->x, self->y, self->z); 
}


/**
 *
 * __neg__ implementation
 *
 **/
static PyObject * Quaternion_neg(Quaternion *self) {
    return (PyObject *)Quaternion_create(-self->w, -self->x, -self->y, -self->z);
}


/**
 *
 * __add__ implementation
 *
 **/
static PyObject * Quaternion_add(Quaternion *self, PyObject *operand) {
    Quaternion *other;

    if (PyObject_TypeCheck(operand, &QuaternionType)) {
        other = (Quaternion *)operand;

        return (PyObject *)Quaternion_create(self->w + other->w, self->x + other->x, self->y + other->y, self->z + other->z);
    }
    PyErr_SetString(PyExc_TypeError, "Expecting Quaternion type");

    return NULL;
}


/**
 *
 * __sub__ implementation
 *
 **/
static PyObject * Quaternion_sub(Quaternion *self, PyObject *operand) {
     Quaternion *other;

     if (PyObject_TypeCheck(operand, &QuaternionType)) {
        other = (Quaternion *)operand;

        return (PyObject *)Quaternion_create(self->w - other->w, self->x - other->x, self->y - other->y, self->z - other->z);
    }
    PyErr_SetString(PyExc_TypeError, "Expecting Quaternion type");

    return NULL;
}


/**
 *
 * __mul__ and __rmul__ implementation
 *
 **/
static PyObject * Quaternion_mult(PyObject *operand1, PyObject *operand2) {
    Quaternion *self, *other;
    double scalar, w, x, y, z;

    if (PyObject_TypeCheck(operand1, &QuaternionType) && PyObject_TypeCheck(operand2, &QuaternionType)) {
        self = (Quaternion *)operand1;
        other = (Quaternion *)operand2;
        w = self->w*other->w - self->x*other->x - self->y*other->y - self->z*other->z;
        x = self->w*other->x + self->x*other->w + self->y*other->z - self->z*other->y;
        y = self->w*other->y - self->x*other->z + self->y*other->w + self->z*other->x;
        z = self->w*other->z + self->x*other->y - self->y*other->x + self->z*other->w;
    } else if (PyObject_TypeCheck(operand1, &QuaternionType)) {
        scalar = PyFloat_AsDouble(operand2);
        self = (Quaternion *)operand1;
        w = self->w*scalar;
        x = self->x*scalar;
        y = self->y*scalar;
        z = self->z*scalar;
    } else {
        scalar = PyFloat_AsDouble(operand1);
        self = (Quaternion *)operand2;
        w = self->w*scalar;
        x = self->x*scalar;
        y = self->y*scalar;
        z = self->z*scalar;
    }

    return (PyObject *)Quaternion_create(w, x, y, z);
}


/**
 *
 * conjugate implementation
 *
 **/
static PyObject * Quaternion_conjugate(Quaternion *self) {
    return (PyObject *)Quaternion_create(self->w, -self->x, -self->y, -self->z);
}


/**
 *
 * __div__ implementation
 *
 **/
static PyObject * Quaternion_div(Quaternion *self, PyObject *operand) {
    double scalar;

    if (! PyObject_TypeCheck(operand, &QuaternionType)) {
        scalar = PyFloat_AsDouble(operand);

        return (PyObject *)Quaternion_create(self->w/scalar, self->x/scalar, self->y/scalar, self->z/scalar);
    }
    PyErr_SetString(PyExc_TypeError, "Quaternion type not valid");

    return NULL;
}


/**
 *
 * dot implementation
 *
 **/
static PyObject * Quaternion_dot(Quaternion *self, PyObject *args) {
    PyObject *result = NULL;
    Quaternion *other;
    double res;

    if (PyArg_ParseTuple(args, "O!", &QuaternionType, &other)) {
        res = dot(self->w, self->x, self->y, self->z, other->w, other->x, other->y, other->z);
        result = PyFloat_FromDouble(res);
    }

    return result;
}


/**
 *
 * norm2 implementation
 *
 **/
static PyObject * Quaternion_norm2(Quaternion *self) {
    return PyFloat_FromDouble(norm2(self->w, self->x, self->y, self->z));
}


/**
 *
 * norm implementation
 *
 **/
static PyObject * Quaternion_norm(Quaternion *self) {
    return PyFloat_FromDouble(norm(self->w, self->x, self->y, self->z));
}


/**
 *
 * normalize implementation
 *
 **/
static PyObject * Quaternion_normalize(Quaternion *self) {
    double normal = norm(self->w, self->x, self->y, self->z);

    return (PyObject *)Quaternion_create(self->w/normal, self->x/normal, self->y/normal, self->z/normal);
} 


/**
 *
 * inverse implementation
 *
 **/
static PyObject * Quaternion_inverse(Quaternion *self) {
    double normal2 = norm2(self->w, self->x, self->y, self->z);

    return (PyObject *)Quaternion_create(self->w/normal2, -self->x/normal2, -self->y/normal2, -self->z/normal2);
}


/**
 *
 * rotate implementation
 *
 **/
static PyObject * Quaternion_rotate(Quaternion *self, PyObject *args) {
    PyObject *tmp0, *tmp1, *result = NULL;
    Quaternion *v;
    double x, y, z;

    if (PyArg_ParseTuple(args, "(ddd)", &x, &y, &z)) {
        v = (Quaternion *)Quaternion_create(0., x, y, z);
        tmp0 = Quaternion_inverse(self);
        tmp1 = Quaternion_mult((PyObject *)v, tmp0);
        Py_DECREF(tmp0);
        Py_DECREF(v);
        v = (Quaternion *)Quaternion_mult((PyObject *)self, tmp1);
        Py_DECREF(tmp1);
        result = Py_BuildValue("(ddd)", v->x, v->y, v->z);
        Py_DECREF(v);
    }

    return result;
}


/**
 *
 * __repr__ implementation
 *
 **/
static PyObject * Quaternion_repr(Quaternion *self) {
    char repr[512];

    sprintf(repr, "(%10.8F, %10.8F, %10.8F, %10.8F)", self->w, self->x, self->y, self->z);

    return Py_BuildValue("s", repr);
}


/**
 *
 * lerp implementation
 *
 **/
static PyObject * Quaternion_lerp(Quaternion *self, PyObject *args) {
    PyObject *tmp0, *tmp1, *tmp2;
    Quaternion *other, *result = NULL;
    double t;

    if (PyArg_ParseTuple(args, "O!d", &QuaternionType, &other, &t)) {
        tmp0 = PyFloat_FromDouble(t);
        tmp1 = Quaternion_mult(tmp0, (PyObject *)other);
        Py_DECREF(tmp0);
        tmp0 = PyFloat_FromDouble(1.0 - t);
        tmp2 = Quaternion_mult(tmp0, (PyObject *)self);
        Py_DECREF(tmp0);
        result = (Quaternion *)Quaternion_add((Quaternion *)tmp2, tmp1);
        Py_DECREF(tmp1);
        Py_DECREF(tmp2);
    }

    return (PyObject *)result;
}


/**
 *
 * slerp implementation
 *
 **/
static PyObject * Quaternion_slerp(Quaternion *self, PyObject *args) {
    PyObject *DEFAULT_ROTATION_STEP, *LINEAR_THRESHOLD, *tmp0, *tmp1, *tmp2;
    double q_dot, omega, so, t, thres;
    Quaternion *other, *result = NULL;

    tmp0 = PyImport_ImportModule("lightdock.constants");
    DEFAULT_ROTATION_STEP = PyObject_GetAttrString(tmp0, "DEFAULT_ROTATION_STEP");
    Py_DECREF(tmp0);
    tmp0 = PyImport_ImportModule("lightdock.mathutil.constants");
    LINEAR_THRESHOLD = PyObject_GetAttrString(tmp0, "LINEAR_THRESHOLD");
    Py_DECREF(tmp0);
    t = PyFloat_AsDouble(DEFAULT_ROTATION_STEP);
    Py_DECREF(DEFAULT_ROTATION_STEP);
    thres = PyFloat_AsDouble(LINEAR_THRESHOLD);
    Py_DECREF(LINEAR_THRESHOLD);

    if (PyArg_ParseTuple(args, "O!|d", &QuaternionType, &other, &t)) {
        self = (Quaternion *)Quaternion_normalize(self);
        other = (Quaternion *)Quaternion_normalize(other);
        tmp0 = Py_BuildValue("(O)", other);
        q_dot = PyFloat_AsDouble(Quaternion_dot(self, Py_BuildValue("(O)", other)));
        Py_DECREF(tmp0);

        if (q_dot < 0) {
            tmp0 = Quaternion_neg(self);
            Py_DECREF(self);
            self = (Quaternion *)tmp0;
            q_dot *= -1.;
        }

        if (q_dot > thres) {
            tmp0 = PyFloat_FromDouble(t);
            tmp1 = Quaternion_sub(other, (PyObject *)self);
            tmp2 = Quaternion_mult(tmp0, tmp1);
            Py_DECREF(tmp0);
            Py_DECREF(tmp1);
            tmp0 = Quaternion_add(self, tmp2);
            Py_DECREF(tmp2);
            result = (Quaternion *)Quaternion_normalize((Quaternion *)tmp0);
            Py_DECREF(tmp0);
        } else {
            q_dot = fmax(fmin(q_dot, 1.0), -1.0);
            omega = acos(q_dot);
            so = sin(omega);
            tmp0 = PyFloat_FromDouble(sin((1.0 - t)*omega)/so);
            tmp1 = Quaternion_mult(tmp0, (PyObject *)self);
            Py_DECREF(tmp0);
            tmp0 = PyFloat_FromDouble(sin(t*omega)/so);
            tmp2 = Quaternion_mult(tmp0, (PyObject *)other);
            Py_DECREF(tmp0);
            result = (Quaternion *)Quaternion_add((Quaternion *)tmp1, tmp2);
            Py_DECREF(tmp1);
            Py_DECREF(tmp2);
        } 
        Py_DECREF(self);
        Py_DECREF(other);
    }

    return (PyObject *)result;
}


/**
 *
 * distance implementation
 *
 **/
static PyObject * Quaternion_distance(Quaternion *self, PyObject *args) {
    PyObject  *tmp0, *tmp1, *result = NULL;
    Quaternion *other;

    if (PyArg_ParseTuple(args, "O!", &QuaternionType, &other)) {
        tmp0 = Py_BuildValue("(O)", other);
        tmp1 = Quaternion_dot(self, tmp0);
        Py_DECREF(tmp0);
        result = PyFloat_FromDouble(1.0 - pow(PyFloat_AsDouble(tmp1), 2));
        Py_DECREF(tmp1);
    }

    return result;
}


/**
 *
 * random implementation
 *
 **/
static PyObject * Quaternion_random(Quaternion *self, PyObject *args) {
    PyObject *rng = NULL;
    Quaternion *result = NULL;
    double u1, u2, u3;

    if (PyArg_ParseTuple(args, "|O", &rng)) {
        if (rng) {
            u1 = PyFloat_AsDouble(PyObject_CallObject(rng, NULL));
            u2 = PyFloat_AsDouble(PyObject_CallObject(rng, NULL));
            u3 = PyFloat_AsDouble(PyObject_CallObject(rng, NULL));
        } else {
            srand((unsigned)time(NULL));
            u1 = ((double)rand()/(double)RAND_MAX);
            u2 = ((double)rand()/(double)RAND_MAX);
            u3 = ((double)rand()/(double)RAND_MAX);
        }
        result = Quaternion_create(sqrt(1-u1)*sin(2*M_PI*u2), sqrt(1-u1)*cos(2*M_PI*u2), sqrt(u1)*sin(2*M_PI*u3), sqrt(u1)*cos(2*M_PI*u3));
    }

    return (PyObject *)result;
}


/**
 *
 * Quaternion type methods
 *
 **/
static PyMethodDef Quaternion_methods[] = {
    {"clone", (PyCFunction)Quaternion_clone, METH_NOARGS, "Creates a new instance of this quaternion"},
    {"conjugate", (PyCFunction)Quaternion_conjugate, METH_NOARGS, "Calculates the conjugate of this quaternion"},
    {"dot", (PyCFunction)Quaternion_dot, METH_VARARGS, "Calculates the dot product of two quaternions"},
    {"norm2", (PyCFunction)Quaternion_norm2, METH_NOARGS, "Calculates quaternion norm^2"},
    {"norm", (PyCFunction)Quaternion_norm, METH_NOARGS, "Calculates quaternion norm"},
    {"normalize", (PyCFunction)Quaternion_normalize, METH_NOARGS, "Normalizes a given quaternion"},
    {"inverse", (PyCFunction)Quaternion_inverse, METH_NOARGS, "Calculates the inverse of this quaternion"},
    {"rotate", (PyCFunction)Quaternion_rotate, METH_VARARGS | METH_KEYWORDS, "Rotates vec3 using quaternion"},
    {"lerp", (PyCFunction)Quaternion_lerp, METH_VARARGS, "Calculates the linear interpolation between two quaternions"},
    {"slerp", (PyCFunction)Quaternion_slerp, METH_VARARGS, "Calculates the spherical linear interpolation of two quaternions given a t step"},
    {"distance", (PyCFunction)Quaternion_distance, METH_VARARGS | METH_KEYWORDS, "Calculates the closeness of two orientations represented in quaternions space.\n\
     \n\
     Quaternions must be normalized. Distance is 0 when quaternions are equal and 1 when\n\
     the orientations are 180 degrees apart.\n\
     See http://math.stackexchange.com/questions/90081/quaternion-distance"},
    {"random", (PyCFunction)Quaternion_random, METH_VARARGS | METH_STATIC, "Generates a random quaternion uniformly distributed:\n\
     http://planning.cs.uiuc.edu/node198.html"},
    {NULL}  /* Sentinel */
};


/**
 *
 * Quaternion type number methods
 *
 **/
static PyNumberMethods quaternion_as_number = {
    (binaryfunc)Quaternion_add,     /* binaryfunc nb_add;          __add__      */
    (binaryfunc)Quaternion_sub,     /* binaryfunc nb_subtract;     __sub__      */
    (binaryfunc)Quaternion_mult,    /* binaryfunc nb_multiply;     __mul__      */
    (binaryfunc)Quaternion_div,     /* binaryfunc nb_divide;       __div__      */
    0,                              /* binaryfunc nb_remainder;    __mod__      */
    0,                              /* binaryfunc nb_divmod;       __divmod__   */
    0,                              /* ternaryfunc nb_power;       __pow__      */
    (unaryfunc)Quaternion_neg,      /* unaryfunc nb_negative;      __neg__      */
    0,                              /* unaryfunc nb_positive;      __pos__      */
    0,                              /* unaryfunc nb_absolute;      __abs__      */
    0,                              /* inquiry nb_nonzero;         __nonzero__  */
    0,                              /* unaryfunc nb_invert;        __invert__   */
    0,                              /* binaryfunc nb_lshift;       __lshift__   */
    0,                              /* binaryfunc nb_rshift;       __rshift__   */
    0,                              /* binaryfunc nb_and;          __and__      */
    0,                              /* binaryfunc nb_xor;          __xor__      */
    0,                              /* binaryfunc nb_or;           __or__       */
    0,                              /* coercion nb_coerce;         __coerce__   */
    0,                              /* unaryfunc nb_int;           __int__      */
    0,                              /* unaryfunc nb_long;          __long__     */
    0,                              /* unaryfunc nb_float;         __float__    */
    0,                              /* unaryfunc nb_oct;           __oct__      */
    0,                              /* unaryfunc nb_hex;           __hex__      */
};


/**
 *
 * Quaternion type definition
 *
 **/
static PyTypeObject QuaternionType = {
    PyObject_HEAD_INIT(NULL)
    0,                              /*ob_size*/
    "Quaternion",                   /*tp_name*/
    sizeof(Quaternion),             /*tp_basicsize*/
    0,                              /*tp_itemsize*/
    (destructor)Quaternion_dealloc, /*tp_dealloc*/
    0,                              /*tp_print*/
    0,                              /*tp_getattr*/
    0,                              /*tp_setattr*/
    0,                              /*tp_compare*/
    (reprfunc)Quaternion_repr,      /*tp_repr*/
    &quaternion_as_number,          /*tp_as_number*/
    0,                              /*tp_as_sequence*/
    0,                              /*tp_as_mapping*/
    0,                              /*tp_hash */
    0,                              /*tp_call*/
    0,                              /*tp_str*/
    0,                              /*tp_getattro*/
    0,                              /*tp_setattro*/
    0,                              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_RICHCOMPARE, /*tp_flags*/
    "Quaternion object",            /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    (richcmpfunc)&Quaternion_richcompare, /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    Quaternion_methods,             /* tp_methods */
    Quaternion_members,             /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc)Quaternion_init,      /* tp_init */
    0,                              /* tp_alloc */
    0,                              /* tp_new */
};


/**
 *
 * Module methods table
 *
 **/
static PyMethodDef module_methods[] = {
    {NULL}
};


/**
 *
 * Initialization function
 *
 **/
PyMODINIT_FUNC initquaternion(void) {
    QuaternionType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&QuaternionType) < 0) {
        return;
    }

    PyObject *m = Py_InitModule3("quaternion", module_methods, "Quaternion object");

    PyModule_AddObject(m, "Quaternion", (PyObject *)&QuaternionType);

    QuaternionError = PyErr_NewException("quaternion.error", NULL, NULL);
    PyModule_AddObject(m, "error", QuaternionError);
}

