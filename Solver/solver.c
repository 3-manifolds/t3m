#include "Python.h"
#include "vertex.h"

static PyObject *ErrorObject;

static PyObject *t3m_find_vertices(PyObject *self, PyObject *args, PyObject *keywds){
  PyObject *result;
  PyObject *pymatrix;
  int rows, columns, length;
  int modp = 0;
  static char *kwlist[] = {"rows", "columns", "matrix", "modp"};

  if ( !PyArg_ParseTupleAndKeywords(args, keywds, "iiO|i:find_vertices", kwlist,
			 &rows, &columns, &pymatrix, &modp) )
    return NULL;

  if ( !PySequence_Check(pymatrix) ){
    PyErr_SetString( ErrorObject, 
     "Argument 3 to find_vertices must support the sequence protocol.");
    return NULL;
  }

  length = PySequence_Length(pymatrix);

  if ( rows < 0 || columns < 0 || length != rows*columns ){
    PyErr_SetString( ErrorObject,
		     "Bad arguments to find_vertices: rows*columns != size of matrix.");
    return NULL;
  }

  {
    int i;
    PyObject *Item;
    matrix_t *matrix = new_matrix(rows, columns);
    filter_list_t *filter = embedded_filter(columns/3);
    for (i=0; i< length; i++) {
      Item = PySequence_GetItem(pymatrix, i);
      matrix->matrix[i] = PyInt_AsLong(Item);
      Py_DECREF(Item);
      }

    if (modp)
      result = find_vertices_mod_p(matrix, filter);
    else
      result = find_vertices(matrix, filter);
  }
  return result;
}

/* Documentation */
static char find_vertices_doc[]=
"Someone should document this.\n"
"\n";

/* List of functions defined in the module */

static PyMethodDef vertex_methods[] = {
  {"find_vertices", (PyCFunction)t3m_find_vertices, METH_VARARGS|METH_KEYWORDS,
   find_vertices_doc},
  {NULL,		NULL}		/* sentinel */
};


/* Initialization function for the module */

DL_EXPORT(void)
initsolver(void)
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule("solver", vertex_methods);

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	ErrorObject = PyErr_NewException("solver.error", NULL, NULL);
	PyDict_SetItemString(d, "error", ErrorObject);
}
