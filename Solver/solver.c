#include "Python.h"
#include "vertex.h"

static PyObject *ErrorObject;

static void *build_vertex_list(vertex_stack_t *stack, int dimension){
  PyObject *result, *coeff;
  vertex_t *V = *stack;
  int i;

  result = PyList_New(0);
  if (result != NULL){
    for (; V != NULL; V = V->next ) {
      coeff = PyTuple_New(dimension);
      for (i=0; i<dimension ; i++)
	PyTuple_SetItem(coeff, i, PyInt_FromLong((long)V->vector[i]));
      PyList_Append(result, coeff);
    }
  }
  return result;
}

static PyObject *t3m_find_vertices(PyObject *self, PyObject *args, PyObject *keywds){
  PyObject *result;
  PyObject *pymatrix;
  int rows, columns, length;
  int modp = 0;
  int filtering = 1;
  static char *kwlist[] = {"rows", "columns", "matrix", "modp", "filtering"};

  if ( !PyArg_ParseTupleAndKeywords(args, keywds, "iiO|ii:find_vertices", kwlist,
			 &rows, &columns, &pymatrix, &modp, &filtering) )
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
    filter_list_t *filter = NULL;
    if (filtering)
      filter = embedded_filter(columns/3);
    for (i=0; i< length; i++) {
      Item = PySequence_GetItem(pymatrix, i);
      matrix->matrix[i] = PyInt_AsLong(Item);
      Py_DECREF(Item);
      }

    if (modp)
      result = find_vertices_mod_p(matrix, filter, build_vertex_list);
    else
      result = find_vertices(matrix, filter, build_vertex_list);
    if (filter)
      destroy_filter_list(filter);
    destroy_matrix(matrix);
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
