#include "Python.h"

static PyObject *ErrorObject;

static PyObject *find_vertices(PyObject *self, PyObject *args){
  PyObject *result;
  PyObject *matrix;
  int length, rows, columns;

  if ( !PyArg_ParseTuple(args, "iiO:find_vertices", &rows, &columns, &matrix) )
    return NULL;

  if ( !PySequence_Check(matrix) ){
    PyErr_SetString( ErrorObject, 
     "Argument 3 to find_vertices must support the sequence protocol.");
    return NULL;
  }

  length = PySequence_Length(matrix);

  if ( rows < 0 || columns < 0 || length != rows*columns ){
    PyErr_SetString( ErrorObject,
		     "Bad arguments to find_vertices: rows*columns != size of matrix.");
    return NULL;
  }

  result = PyList_New(0);
  {
    int i;
    PyObject *Item;
    int A[length];
    for (i=0; i< length; i++) {
      Item = PySequence_GetItem(matrix, i);
      A[i] = PyInt_AsLong(Item);
      Py_DECREF(Item);
      }

  }
  return result;
}


/* Documentation */
static char find_vertices_doc[]=
"Someone should document this.\n"
"\n";

/* List of functions defined in the module */

static PyMethodDef vertex_methods[] = {
	{"find_vertices",  find_vertices,  METH_VARARGS, find_vertices_doc },
	{NULL,		NULL}		/* sentinel */
};


/* Initialization function for the module */

DL_EXPORT(void)
init_normal(void)
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule("_normal", vertex_methods);

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	ErrorObject = PyErr_NewException("_normal.error", NULL, NULL);
	PyDict_SetItemString(d, "error", ErrorObject);
}
