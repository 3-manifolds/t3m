//  !!!! Modified for t3m by MC
//
//  SnapPeaC.c
//
//  The Python interface to the SnapPea kernel consists of two parts.
//
//      SnapPea.py (in a different file)(this file) defines a set
//          of objects (Triangulation, AbelianGroup, etc.) purely in Python.
//
//      SnapPeaC.c (this file) implements SnapPea.py's methods as a set
//          of wrappers to the standard SnapPea kernel functions,
//          which are written in C.
//
//  Technical comment:  There's the awkward question of how the Python
//  interpreter is to keep copies of pointers to SnapPea data structures,
//  such as Triangulations, which SnapPea.h exports as "opaque typedefs".
//  The present solution is to pass them to the interpreter as long ints.

//  If your compiler can't find Python.h, please go to the makefile
//  and adjust the line
// 
//      PYTHON_DIRECTORY = /usr/include/python1.5
//
//  If you don't know where Python.h is, try "locate Python.h".
#include <Python.h>

#include <ctype.h>
#include "assert.h"
#include "SnapPea.h"

//  from SnapPea's unix kit
#include "unix_cusped_census.h"
#include "unix_file_io.h"

// start t3m
// declare this here for now
Triangulation* easy_triangulate_punctured_torus_bundle(
               Boolean negative_determinant, 
	       Boolean negative_trace, char* LR_factors);
// end t3m

#define FALSE   0
#define TRUE    1

       void initSnapPeaC(void);
static char *RelationToString(int *aSnapPeaRelation);


#define DECLARE_WRAPPER(name)   \
    static PyObject *name(PyObject *self, PyObject *args);

DECLARE_WRAPPER(wrap_verify_my_malloc_usage);
DECLARE_WRAPPER(wrap_get_triangulation);
DECLARE_WRAPPER(wrap_save_triangulation);
DECLARE_WRAPPER(wrap_copy_triangulation);
DECLARE_WRAPPER(wrap_get_cusped_census_manifold);
DECLARE_WRAPPER(wrap_free_triangulation);
DECLARE_WRAPPER(wrap_get_triangulation_name);
DECLARE_WRAPPER(wrap_set_triangulation_name);
DECLARE_WRAPPER(wrap_get_triangulation_is_orientable);
DECLARE_WRAPPER(wrap_get_solution_type);
DECLARE_WRAPPER(wrap_get_num_tetrahedra);
DECLARE_WRAPPER(wrap_get_num_cusps);
DECLARE_WRAPPER(wrap_get_cusp_is_complete);
DECLARE_WRAPPER(wrap_get_cusp_is_orientable);
DECLARE_WRAPPER(wrap_get_cusp_m);
DECLARE_WRAPPER(wrap_get_cusp_l);
DECLARE_WRAPPER(wrap_set_cusp_info);
DECLARE_WRAPPER(wrap_cusp_is_fillable);
DECLARE_WRAPPER(wrap_remove_Dehn_fillings);
DECLARE_WRAPPER(wrap_volume);
DECLARE_WRAPPER(wrap_homology);
DECLARE_WRAPPER(wrap_fill_cusp);
DECLARE_WRAPPER(wrap_get_drillable_curves);
DECLARE_WRAPPER(wrap_drill_curve);
DECLARE_WRAPPER(wrap_get_normal_surfaces);
DECLARE_WRAPPER(wrap_split_along_normal_surface);
DECLARE_WRAPPER(wrap_fundamental_group);
DECLARE_WRAPPER(wrap_free_group_presentation);
DECLARE_WRAPPER(wrap_fg_get_num_generators);
DECLARE_WRAPPER(wrap_fg_get_relations);
DECLARE_WRAPPER(wrap_fg_representation);
DECLARE_WRAPPER(wrap_fg_peripheral_curves);
DECLARE_WRAPPER(wrap_core_geodesic);
DECLARE_WRAPPER(wrap_shortest_curves_become_meridians);
DECLARE_WRAPPER(wrap_current_fillings_become_meridians);
DECLARE_WRAPPER(wrap_basic_simplification);
DECLARE_WRAPPER(wrap_randomize_triangulation);
DECLARE_WRAPPER(wrap_reorient);
DECLARE_WRAPPER(wrap_proto_canonize);
DECLARE_WRAPPER(wrap_is_canonical_triangulation);
DECLARE_WRAPPER(wrap_symmetry_group);
DECLARE_WRAPPER(wrap_free_symmetry_group);
DECLARE_WRAPPER(wrap_symmetry_group_order);
DECLARE_WRAPPER(wrap_symmetry_group_is_abelian);
DECLARE_WRAPPER(wrap_symmetry_group_abelian_description);
DECLARE_WRAPPER(wrap_symmetry_group_is_dihedral);
DECLARE_WRAPPER(wrap_symmetry_group_is_polyhedral);
DECLARE_WRAPPER(wrap_symmetry_group_polyhedral_description);
DECLARE_WRAPPER(wrap_symmetry_group_is_S5);
DECLARE_WRAPPER(wrap_symmetry_group_is_direct_product);
DECLARE_WRAPPER(wrap_symmetry_group_factor);
DECLARE_WRAPPER(wrap_symmetry_group_is_amphicheiral);
DECLARE_WRAPPER(wrap_symmetry_group_invertible_knot);
DECLARE_WRAPPER(wrap_symmetry_group_commutator_subgroup);
DECLARE_WRAPPER(wrap_symmetry_group_abelianization);
DECLARE_WRAPPER(wrap_symmetry_group_center);
DECLARE_WRAPPER(wrap_symmetry_group_presentation);
DECLARE_WRAPPER(wrap_tet_shapes);
DECLARE_WRAPPER(wrap_Dirichlet);
DECLARE_WRAPPER(wrap_free_Dirichlet_domain);
DECLARE_WRAPPER(wrap_Dirichlet_num_vertices);
DECLARE_WRAPPER(wrap_Dirichlet_num_edges);
DECLARE_WRAPPER(wrap_Dirichlet_num_faces);
DECLARE_WRAPPER(wrap_Dirichlet_vertices);
DECLARE_WRAPPER(wrap_Dirichlet_faces);
DECLARE_WRAPPER(wrap_Dirichlet_face_colors);
DECLARE_WRAPPER(wrap_Dirichlet_face_pairings);
// start t3m
DECLARE_WRAPPER(t3m_get_gluing_data);
DECLARE_WRAPPER(t3m_get_triangulation_from_DT);
DECLARE_WRAPPER(t3m_easy_triangulate_punctured_torus_bundle);
// end t3m

static struct PyMethodDef SnapPeaCMethods[] =
{
    {"verify_my_malloc_usage",                  wrap_verify_my_malloc_usage,                METH_VARARGS},
    {"get_triangulation",                       wrap_get_triangulation,                     METH_VARARGS},
    {"save_triangulation",                      wrap_save_triangulation,                    METH_VARARGS},
    {"copy_triangulation",                      wrap_copy_triangulation,                    METH_VARARGS},
    {"get_cusped_census_manifold",              wrap_get_cusped_census_manifold,            METH_VARARGS},
    {"free_triangulation",                      wrap_free_triangulation,                    METH_VARARGS},
    {"get_triangulation_name",                  wrap_get_triangulation_name,                METH_VARARGS},
    {"set_triangulation_name",                  wrap_set_triangulation_name,                METH_VARARGS},
    {"get_triangulation_is_orientable",         wrap_get_triangulation_is_orientable,       METH_VARARGS},
    {"get_solution_type",                       wrap_get_solution_type,                     METH_VARARGS},
    {"get_num_tetrahedra",                      wrap_get_num_tetrahedra,                    METH_VARARGS},
    {"get_num_cusps",                           wrap_get_num_cusps,                         METH_VARARGS},
    {"get_cusp_is_complete",                    wrap_get_cusp_is_complete,                  METH_VARARGS},
    {"get_cusp_is_orientable",                  wrap_get_cusp_is_orientable,                METH_VARARGS},
    {"get_cusp_m",                              wrap_get_cusp_m,                            METH_VARARGS},
    {"get_cusp_l",                              wrap_get_cusp_l,                            METH_VARARGS},
    {"set_cusp_info",                           wrap_set_cusp_info,                         METH_VARARGS},
    {"cusp_is_fillable",                        wrap_cusp_is_fillable,                      METH_VARARGS},
    {"remove_Dehn_fillings",                    wrap_remove_Dehn_fillings,                  METH_VARARGS},
    {"volume",                                  wrap_volume,                                METH_VARARGS},
    {"homology",                                wrap_homology,                              METH_VARARGS},
    {"fill_cusp",                               wrap_fill_cusp,                             METH_VARARGS},
    {"get_drillable_curves",                    wrap_get_drillable_curves,                  METH_VARARGS},
    {"drill_curve",                             wrap_drill_curve,                           METH_VARARGS},
    {"get_normal_surfaces",                     wrap_get_normal_surfaces,                   METH_VARARGS},
    {"split_along_normal_surface",              wrap_split_along_normal_surface,            METH_VARARGS},
    {"fundamental_group",                       wrap_fundamental_group,                     METH_VARARGS},
    {"free_group_presentation",                 wrap_free_group_presentation,               METH_VARARGS},
    {"fg_get_num_generators",                   wrap_fg_get_num_generators,                 METH_VARARGS},
    {"fg_get_relations",                        wrap_fg_get_relations,                      METH_VARARGS},
    {"fg_representation",                       wrap_fg_representation,                     METH_VARARGS},
    {"fg_peripheral_curves",                    wrap_fg_peripheral_curves,                  METH_VARARGS},
    {"core_geodesic",                           wrap_core_geodesic,                         METH_VARARGS},
    {"shortest_curves_become_meridians",        wrap_shortest_curves_become_meridians,      METH_VARARGS},
    {"current_fillings_become_meridians",       wrap_current_fillings_become_meridians,     METH_VARARGS},
    {"basic_simplification",                    wrap_basic_simplification,                  METH_VARARGS},
    {"randomize_triangulation",                 wrap_randomize_triangulation,               METH_VARARGS},
    {"reorient",                                wrap_reorient,                              METH_VARARGS},
    {"proto_canonize",                          wrap_proto_canonize,                        METH_VARARGS},
    {"is_canonical_triangulation",              wrap_is_canonical_triangulation,            METH_VARARGS},
    {"symmetry_group",                          wrap_symmetry_group,                        METH_VARARGS},
    {"free_symmetry_group",                     wrap_free_symmetry_group,                   METH_VARARGS},
    {"symmetry_group_order",                    wrap_symmetry_group_order,                  METH_VARARGS},
    {"symmetry_group_is_abelian",               wrap_symmetry_group_is_abelian,             METH_VARARGS},
    {"symmetry_group_abelian_description",      wrap_symmetry_group_abelian_description,    METH_VARARGS},
    {"symmetry_group_is_dihedral",              wrap_symmetry_group_is_dihedral,            METH_VARARGS},
    {"symmetry_group_is_polyhedral",            wrap_symmetry_group_is_polyhedral,          METH_VARARGS},
    {"symmetry_group_polyhedral_description",   wrap_symmetry_group_polyhedral_description, METH_VARARGS},
    {"symmetry_group_is_S5",                    wrap_symmetry_group_is_S5,                  METH_VARARGS},
    {"symmetry_group_is_direct_product",        wrap_symmetry_group_is_direct_product,      METH_VARARGS},
    {"symmetry_group_factor",                   wrap_symmetry_group_factor,                 METH_VARARGS},
    {"symmetry_group_is_amphicheiral",          wrap_symmetry_group_is_amphicheiral,        METH_VARARGS},
    {"symmetry_group_invertible_knot",          wrap_symmetry_group_invertible_knot,        METH_VARARGS},
    {"symmetry_group_commutator_subgroup",      wrap_symmetry_group_commutator_subgroup,    METH_VARARGS},
    {"symmetry_group_abelianization",           wrap_symmetry_group_abelianization,         METH_VARARGS},
    {"symmetry_group_center",                   wrap_symmetry_group_center,                 METH_VARARGS},
    {"symmetry_group_presentation",             wrap_symmetry_group_presentation,           METH_VARARGS},
    {"tet_shapes",                              wrap_tet_shapes,                            METH_VARARGS},
    {"Dirichlet",                               wrap_Dirichlet,                             METH_VARARGS},
    {"free_Dirichlet_domain",                   wrap_free_Dirichlet_domain,                 METH_VARARGS},
    {"Dirichlet_num_vertices",                  wrap_Dirichlet_num_vertices,                METH_VARARGS},
    {"Dirichlet_num_edges",                     wrap_Dirichlet_num_edges,                   METH_VARARGS},
    {"Dirichlet_num_faces",                     wrap_Dirichlet_num_faces,                   METH_VARARGS},
    {"Dirichlet_vertices",                      wrap_Dirichlet_vertices,                    METH_VARARGS},
    {"Dirichlet_faces",                         wrap_Dirichlet_faces,                       METH_VARARGS},
    {"Dirichlet_face_colors",                   wrap_Dirichlet_face_colors,                 METH_VARARGS},
    {"Dirichlet_face_pairings",                 wrap_Dirichlet_face_pairings,               METH_VARARGS},
// start t3m 
    {"get_gluing_data",                         t3m_get_gluing_data,                        METH_VARARGS},
    {"get_triangulation_from_DT",               t3m_get_triangulation_from_DT,              METH_VARARGS},
    {"easy_triangulate_punctured_torus_bundle", t3m_easy_triangulate_punctured_torus_bundle,METH_VARARGS},
// end t3m 
    {NULL,                                      NULL                                                    }
};


void initSnapPeaC(void)
{
    (void) Py_InitModule("SnapPeaC", SnapPeaCMethods);
}

// start t3m
static PyObject *t3m_get_gluing_data(PyObject *self, PyObject *args)
{
    long theTriangulation;
    int fill_first, fill_how_many = 0, num_cusps;
    Triangulation *filled_Triangulation = NULL;
    TriangulationData *data = NULL;
    int i,j,k; 
    PyObject *gluing_data, *neighbors, *gluings, *permutation, *tetrahedron_datum;

    if (!PyArg_ParseTuple(args, "li:t3m_get_gluing_data",
			  &theTriangulation, &fill_first))
    {
        PyErr_SetString(PyExc_TypeError, 
			"Invalid arguments passed to SnapPeaC.get_gluing_data.\n"
			"           Expected ((long)triangulation, (int)fill_first).");
        return NULL;
    }

    if (fill_first) {
      num_cusps = get_num_cusps((Triangulation *)theTriangulation);
      {
	Boolean fill_list[num_cusps];

	for (i = 0; i<num_cusps; i++)
	  if ( (fill_list[i] = cusp_is_fillable((Triangulation *)theTriangulation, i)) )
	    ++fill_how_many;

	if (fill_how_many > 0)
	  {
	    if (fill_how_many < num_cusps)
	      filled_Triangulation = fill_cusps((Triangulation *)theTriangulation,
						fill_list, "", FALSE);
	    else 
	      filled_Triangulation = fill_cusps((Triangulation *)theTriangulation,
						 NULL, "", TRUE);
	  }
      }
      triangulation_to_data(filled_Triangulation, &data);
    }
    else
      triangulation_to_data((Triangulation *)theTriangulation, &data);

    gluing_data = PyList_New(0);
       for (i = 0; i < data->num_tetrahedra; i++)
       {
           neighbors = PyTuple_New(4);
           gluings = PyTuple_New(4);
           for (j = 0; j < 4; j++)
           {
               permutation = PyTuple_New(4);
               for (k=0;k<4;k++)
                  PyTuple_SET_ITEM(permutation, k, PyInt_FromLong(data->tetrahedron_data[i].gluing[j][k]));
               PyTuple_SET_ITEM(neighbors, j, PyInt_FromLong(data->tetrahedron_data[i].neighbor_index[j]));
               PyTuple_SET_ITEM(gluings, j, permutation);
           }
           tetrahedron_datum = PyTuple_New(2);
           PyTuple_SET_ITEM(tetrahedron_datum, 0, neighbors);
           PyTuple_SET_ITEM(tetrahedron_datum, 1, gluings);
           PyList_Append(gluing_data, tetrahedron_datum);
		   Py_DECREF(tetrahedron_datum);
      }
    if (filled_Triangulation) free_triangulation(filled_Triangulation);
    if (data) free_triangulation_data(data);
    return gluing_data;
}

/*
   This function takes as argument a python list which is a
   Dowker-Thistlethwaite description of a knot, and returns a snapPea
   triangulation.  No consistency checking is done on the DT
   description.  That should be done in a python wrapper.
*/
extern Triangulation   *DT_int_to_triangulation(int aNumCrossings, int *aDTCode);

static PyObject *t3m_get_triangulation_from_DT(PyObject *self, PyObject *args)
{
    PyObject        *theDTlist, *aPyInt;
    Triangulation   *theTriangulation;
    int aNumCrossings, *aDTCode, i;

    if (!PyArg_ParseTuple(args, "O", &theDTlist)  || !PyList_Check(theDTlist) )
    {
        PyErr_SetString(PyExc_TypeError, "t3m_get_triangulation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    aNumCrossings = PyList_Size(theDTlist);
    aDTCode = malloc(aNumCrossings * sizeof(int));
    for (i=0; i < aNumCrossings; i++)
    {
        aPyInt = PyList_GetItem(theDTlist, i);
        if (aPyInt == NULL )
        {
            PyErr_SetString(PyExc_RuntimeError, "Invalid DT list.");
            return NULL;
	}
        aDTCode[i] = (int)PyInt_AsLong(aPyInt);
    }
    
    theTriangulation = DT_int_to_triangulation(aNumCrossings, aDTCode);
    
    if (theTriangulation == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not create triangulation.");
        return NULL;
    }
    
    return Py_BuildValue("l", (long int) theTriangulation);
}

extern Triangulation *triangulate_punctured_torus_bundle(
                      LRFactorization *anLRFactorization);

/*
 *      If the manifold is hyperbolic (i.e. if the number of LR factors
 *      is at least two for an orientable bundle, or at least one for a
 *      nonorientable bundle), triangulates the complement and returns
 *      a pointer to it.  Otherwise returns NULL.
 */

Triangulation* easy_triangulate_punctured_torus_bundle(
               Boolean negative_determinant, 
	       Boolean negative_trace, char* LR_factors){
    Triangulation* manifold;
    LRFactorization*  glueing;

    glueing = alloc_LR_factorization(strlen(LR_factors));
    glueing->is_available = TRUE;
    glueing->negative_determinant = negative_determinant;
    glueing->negative_trace = negative_trace;
    strcpy(glueing->LR_factors, LR_factors);
    manifold =  triangulate_punctured_torus_bundle(glueing);
    free_LR_factorization(glueing);
    return manifold;
}

static PyObject *t3m_easy_triangulate_punctured_torus_bundle(PyObject *self, PyObject *args) {
    Triangulation *theTriangulation;
    Boolean  _arg0;
    Boolean  _arg1;
    char * _arg2;

    if(!PyArg_ParseTuple(args,"bbs:easy_triangulate_punctured_torus_bundle",&_arg0,&_arg1,&_arg2)) 
        return NULL;
    theTriangulation = (Triangulation *)easy_triangulate_punctured_torus_bundle(_arg0,_arg1,_arg2);
    return Py_BuildValue("l", (long int) theTriangulation);
}

// end t3m 

static PyObject *wrap_verify_my_malloc_usage(PyObject *self, PyObject *args)
{
    verify_my_malloc_usage();
    return Py_BuildValue("");
}


static PyObject *wrap_get_triangulation(PyObject *self, PyObject *args)
{
    char            *theName;
    Triangulation   *theTriangulation;

    if (!PyArg_ParseTuple(args, "s", &theName))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_triangulation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    theTriangulation = get_triangulation(theName);
    
    if (theTriangulation == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not read named triangulation.");
        return NULL;
    }
    
    return Py_BuildValue("l", (long int) theTriangulation);
}


static PyObject *wrap_save_triangulation(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    char    *theName;

    if (!PyArg_ParseTuple(args, "ls", &theTriangulation, &theName))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_save_triangulation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    save_triangulation((Triangulation *)theTriangulation, theName);
    
    return Py_BuildValue("");
}


static PyObject *wrap_copy_triangulation(PyObject *self, PyObject *args)
{
    long            theTriangulation;
    Triangulation   *theCopy;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_copy_triangulation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    copy_triangulation((Triangulation *)theTriangulation, &theCopy);
    
    if (theCopy == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not copy the triangulation.");
        return NULL;
    }
        
    return Py_BuildValue("l", (long int) theCopy);
}


static PyObject *wrap_get_cusped_census_manifold(PyObject *self, PyObject *args)
{
    int             theCensus;          //  5, 6, or 7
    int             theOrientability;   //  0 (false) or 1 (true)
    int             theIndex;           //  0 through #manifolds - 1
    Triangulation*  theTriangulation;

    if (!PyArg_ParseTuple(args, "iii", &theCensus, &theOrientability, &theIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_cusped_census_manifold() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    theTriangulation = GetCuspedCensusManifold(
        theCensus,
        theOrientability ? oriented_manifold : nonorientable_manifold,
        theIndex);
    
    if (theTriangulation == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not load requested census manifold.");
        return NULL;
    }
    
    return Py_BuildValue("l", (long int) theTriangulation);
}


static PyObject *wrap_free_triangulation(PyObject *self, PyObject *args)
{
    long    theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_free_triangulation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    free_triangulation((Triangulation *) theTriangulation);
    
    return Py_BuildValue("");
}


static PyObject *wrap_get_triangulation_name(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    char    *theName;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_triangulation_name() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    //  We get a pointer to the actual name string (not a copy),
    //  so we won't need to free it.
    theName = get_triangulation_name((Triangulation *) theTriangulation);
    
    return Py_BuildValue("s", theName);
}


static PyObject *wrap_set_triangulation_name(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    char    *theName;

    if (!PyArg_ParseTuple(args, "ls", &theTriangulation, &theName))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_set_triangulation_name() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    set_triangulation_name((Triangulation *)theTriangulation, theName);
    
    return Py_BuildValue("");
}


static PyObject *wrap_get_triangulation_is_orientable(PyObject *self, PyObject *args)
{
    long            theTriangulation;
    Orientability   theOrientability;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_triangulation_is_orientable() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    theOrientability = get_orientability((Triangulation *) theTriangulation);
    if (theOrientability != oriented_manifold
     && theOrientability != nonorientable_manifold)
    {
        PyErr_SetString(PyExc_RuntimeError, "Triangulation's orientability is unknown.");
        return NULL;
    }
    
    return Py_BuildValue("i", theOrientability == oriented_manifold);
}


static PyObject *wrap_get_solution_type(PyObject *self, PyObject *args)
{
    long    theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_solution_type() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    switch (get_filled_solution_type((Triangulation *) theTriangulation))
    {
        case not_attempted:         return Py_BuildValue("s", "not attempted");
        case geometric_solution:    return Py_BuildValue("s", "all tetrahedra positively oriented");
        case nongeometric_solution: return Py_BuildValue("s", "contains negatively oriented tetrahedra");
        case flat_solution:         return Py_BuildValue("s", "contains flat tetrahedra");
        case degenerate_solution:   return Py_BuildValue("s", "contains degenerate tetrahedra");
        case other_solution:        return Py_BuildValue("s", "unrecognized solution type");
        case no_solution:           return Py_BuildValue("s", "no solution found");
        default:
            PyErr_SetString(PyExc_RuntimeError, "Bad solution type.");
            return NULL;
    }
}


static PyObject *wrap_get_num_tetrahedra(PyObject *self, PyObject *args)
{
    long    theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_num_tetrahedra() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    return Py_BuildValue("i", get_num_tetrahedra((Triangulation *) theTriangulation));
}


static PyObject *wrap_get_num_cusps(PyObject *self, PyObject *args)
{
    long    theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_num_cusps() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    return Py_BuildValue("i", get_num_cusps((Triangulation *) theTriangulation));
}


static PyObject *wrap_get_cusp_is_complete(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    int     theCuspIndex;
    Boolean theCuspIsComplete;

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theCuspIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_cusp_is_complete() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    get_cusp_info(  (Triangulation *) theTriangulation, theCuspIndex,
        NULL, &theCuspIsComplete, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    return Py_BuildValue("i", theCuspIsComplete);
}


static PyObject *wrap_get_cusp_is_orientable(PyObject *self, PyObject *args)
{
    long            theTriangulation;
    int             theCuspIndex;
    CuspTopology    theCuspTopology;

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theCuspIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_cusp_is_orientable() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    get_cusp_info(  (Triangulation *) theTriangulation, theCuspIndex,
        &theCuspTopology, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    //  Python scripts should never encounter Cusps of unknown CuspTopology.
    if (theCuspTopology != torus_cusp  &&  theCuspTopology != Klein_cusp)
    {
        PyErr_SetString(PyExc_RuntimeError, "Unknown cusp topology.");
        return NULL;
    }

    return Py_BuildValue("i", theCuspTopology == torus_cusp);
}


static PyObject *wrap_get_cusp_m(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    int     theCuspIndex;
    double  m;

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theCuspIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_cusp_m() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    get_cusp_info(  (Triangulation *) theTriangulation, theCuspIndex,
                    NULL, NULL, &m, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    return Py_BuildValue("d", m);
}


static PyObject *wrap_get_cusp_l(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    int     theCuspIndex;
    double  l;

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theCuspIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_cusp_l() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    get_cusp_info(  (Triangulation *) theTriangulation, theCuspIndex,
                    NULL, NULL, NULL, &l, NULL, NULL, NULL, NULL, NULL, NULL);

    return Py_BuildValue("d", l);
}


static PyObject *wrap_set_cusp_info(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    int     theCuspIndex;
    double  m,
            l;
    int     theRecomputeFlag;

    if (!PyArg_ParseTuple(args, "liddi", &theTriangulation, &theCuspIndex,
                &m, &l, &theRecomputeFlag))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_set_cusp_info() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    set_cusp_info(  (Triangulation *) theTriangulation, theCuspIndex,
                    (m == 0.0 && l == 0.0), m, l);

    if (theRecomputeFlag)
        (void) do_Dehn_filling((Triangulation *) theTriangulation);

    return Py_BuildValue("");
}


static PyObject *wrap_cusp_is_fillable(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    int     theCuspIndex;

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theCuspIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_cusp_is_fillable() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    return Py_BuildValue("i",
        cusp_is_fillable((Triangulation *) theTriangulation, theCuspIndex));
}

static PyObject *wrap_remove_Dehn_fillings(PyObject *self, PyObject *args)
{
    long    theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_remove_Dehn_fillings() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    (void) remove_Dehn_fillings((Triangulation *) theTriangulation);

    return Py_BuildValue("");
}

static PyObject *wrap_volume(PyObject *self, PyObject *args)
{
    long    theTriangulation;
    double  theVolume;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_volume() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    theVolume = volume((Triangulation *) theTriangulation, NULL);
    
    return Py_BuildValue("d", theVolume);
}


static PyObject *wrap_homology(PyObject *self, PyObject *args)
{
    long            theTriangulation;
    AbelianGroup    *theHomology;
    PyObject        *theList;
    int             i;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_homology() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    theHomology = homology((Triangulation *) theTriangulation);
    
    if (theHomology == NULL)
        return Py_BuildValue("");

    compress_abelian_group(theHomology);

    theList = PyList_New(theHomology->num_torsion_coefficients);
    for (i = 0; i < theHomology->num_torsion_coefficients; i++)
        PyList_SetItem(theList, i, Py_BuildValue("l", theHomology->torsion_coefficients[i]));
    
    free_abelian_group(theHomology);

    return theList;
}


static PyObject *wrap_fill_cusp(PyObject *self, PyObject *args)
{
    long            theTri;
    Triangulation   *theTriangulation,
                    *theResult;
    int             theIndex,
                    n,
                    i;
    Boolean         *theFillArray;

    if (!PyArg_ParseTuple(args, "li", &theTri, &theIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_fill_cusp() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    theTriangulation = (Triangulation *) theTri;
    
    if (theIndex < 0 || theIndex >= get_num_cusps(theTriangulation))
    {
        PyErr_SetString(PyExc_RuntimeError, "Cusp index out of bounds.");
        return NULL;
    }
    
    //  If the cusp isn't fillable, return the original Triangulation.
    if (cusp_is_fillable(theTriangulation, theIndex) == FALSE)
        return Py_BuildValue("l", theTri);

    //  Otherwise do the filling, and return the result.
    n = get_num_cusps(theTriangulation);
    theFillArray = (Boolean *) malloc(n * sizeof(Boolean));
    if (theFillArray == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't allocate theFillArray.");
        return NULL;
    }
    for (i = 0; i < n; i++)
        theFillArray[i] = FALSE;
    theFillArray[theIndex] = TRUE;
    theResult = fill_cusps( theTriangulation,
                            theFillArray,
                            get_triangulation_name(theTriangulation),
                            FALSE);
    free(theFillArray);
    
    if (theResult == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't fill the requested cusps.");
        return NULL;
    }
    else
        return Py_BuildValue("l", (long int) theResult);
}


static PyObject *wrap_get_drillable_curves(PyObject *self, PyObject *args)
{
    long                    theTriangulation;
    int                     theMaxSegments,
                            theNumCurves,
                            i;
    DualOneSkeletonCurve    **theCurves;
    MatrixParity            theParity;
    Complex                 theCompleteLength,
                            theFilledLength;
    PyObject                *theReturnData,
                            *theCurveData;

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theMaxSegments))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_drillable_curves() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    dual_curves((Triangulation *) theTriangulation,
                theMaxSegments,
                &theNumCurves,
                &theCurves);
    
    theReturnData = PyList_New(theNumCurves);
    
    for (i = 0; i < theNumCurves; i++)
    {
        get_dual_curve_info(theCurves[i],
                            &theCompleteLength,
                            &theFilledLength,
                            &theParity);
        
        theCurveData = Py_BuildValue("[i,d,d,d,d]",
                                    theParity == orientation_reversing,
                                    theFilledLength.real,
                                    theFilledLength.imag,
                                    theCompleteLength.real,
                                    theCompleteLength.imag);

        PyList_SetItem(theReturnData, i, theCurveData);
        //  PyList_SetItem() takes ownership of the reference to theCurveData,
        //  so we do NOT decrement theCurveData's reference count.
    }
    
    free_dual_curves(theNumCurves, theCurves);
    
    return theReturnData;
}


static PyObject *wrap_drill_curve(PyObject *self, PyObject *args)
{
    long                    theTriangulation;
    int                     theMaxSegments,
                            theIndex,
                            theNumCurves;
    Triangulation           *theResult;
    DualOneSkeletonCurve    **theCurves;

    if (!PyArg_ParseTuple(args, "lii", &theTriangulation, &theIndex, &theMaxSegments))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_drill_curve() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    dual_curves((Triangulation *) theTriangulation,
                theMaxSegments,
                &theNumCurves,
                &theCurves);
    
    if (theIndex < 0 || theIndex >= theNumCurves)
    {
        free_dual_curves(theNumCurves, theCurves);
        PyErr_SetString(PyExc_RuntimeError, "Curve index out of bounds.");
        return NULL;
    }
    
    theResult = drill_cusp( (Triangulation *) theTriangulation,
                            theCurves[theIndex],
                            get_triangulation_name((Triangulation *) theTriangulation) );

    free_dual_curves(theNumCurves, theCurves);
    
    if (theResult == NULL)
    {
        //  theResult will be NULL if the curve is boundary parallel.
        //  At present such curves aren't offered to the user,
        //  so NULL results shouldn't occur.  If we want to allow
        //  them in the future, replace "return NULL" with
        //
        //      return Py_BuildValue("l", (long int) 0);
        //
        PyErr_SetString(PyExc_RuntimeError, "Requested curve cannot be drilled because it's parallel to the boundary.");
        return NULL;
    }
    else
        return Py_BuildValue("l", (long int) theResult);
}


static PyObject *wrap_get_normal_surfaces(PyObject *self, PyObject *args)
{
    long                theTriangulation;
    NormalSurfaceList   *theSurfaceList;
    int                 theNumSurfaces,
                        i;
    PyObject            *theReturnData,
                        *theSurfaceDescription;
    char                theDescription[128];

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_get_normal_surfaces() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    switch (find_normal_surfaces((Triangulation *) theTriangulation, &theSurfaceList))
    {
        case func_OK:

            theNumSurfaces = number_of_normal_surfaces_on_list(theSurfaceList);
            
            theReturnData = PyList_New(theNumSurfaces);

            for (i = 0; i < theNumSurfaces; i++)
            {
                if (normal_surface_is_two_sided(theSurfaceList, i))
                    strcpy(theDescription, "2-sided ");
                else
                    strcpy(theDescription, "1-sided ");

                switch (normal_surface_Euler_characteristic(theSurfaceList, i))
                {
                    case 2:
                        assert(normal_surface_is_orientable(theSurfaceList, i) == TRUE);
                        strcat(theDescription, "sphere");
                        break;

                    case 1:
                        assert(normal_surface_is_orientable(theSurfaceList, i) == FALSE);
                        strcat(theDescription, "projective plane");
                        break;

                    case 0:
                        if (normal_surface_is_orientable(theSurfaceList, i))
                            strcat(theDescription, "torus");
                        else
                            strcat(theDescription, "Klein bottle");
                        break;
                    
                    default:    //  shouldn't ever happen
                        Py_DECREF(theReturnData);
                        free_normal_surfaces(theSurfaceList);
                        PyErr_SetString(PyExc_RuntimeError, "'Impossible' situation:  splitting surface has Euler characteristic not equal to 0, 1, or 2.");
                        return NULL;
                }

                theSurfaceDescription = Py_BuildValue("s", theDescription);

                PyList_SetItem(theReturnData, i, theSurfaceDescription);
                //  PyList_SetItem() takes ownership of the reference
                //  to theSurfaceDescription, so we do NOT decrement
                //  theSurfaceDescription's reference count.
            }

            free_normal_surfaces(theSurfaceList);

            return theReturnData;
        
        case func_bad_input:
            return Py_BuildValue("[s]", "splittings not yet available for closed manifolds");
        
        default:
            return Py_BuildValue("[s]", "???"); // should never occur
    }
}


static PyObject *wrap_split_along_normal_surface(PyObject *self, PyObject *args)
{
    long                theTriangulation;
    int                 theIndex;
    NormalSurfaceList   *theSurfaceList;
    int                 theNumSurfaces;
    Triangulation       *thePieces[2];
    PyObject            *thePieceList;
    char                *theOldName,
                        *theNewName;
    int                 i;

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_split_along_normal_surface() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    if (func_OK != find_normal_surfaces((Triangulation *) theTriangulation,
                                        &theSurfaceList))
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't split along normal surface.  2-sided projective plane?");
        return NULL;
    }

    theNumSurfaces = number_of_normal_surfaces_on_list(theSurfaceList);
    assert(theIndex >= 0  &&  theIndex < theNumSurfaces);
    
    //  Even if split_along_normal_surface() returns func_bad_input,
    //  we can proceed normally.  thePieces[] will both be NULL.
    (void) split_along_normal_surface(theSurfaceList, theIndex, thePieces);
    
    free_normal_surfaces(theSurfaceList);
    
    thePieceList = PyList_New(0);

    theOldName = get_triangulation_name((Triangulation *) theTriangulation);
    for (i = 0; i < 2; i++)
        if (thePieces[i] != NULL)
        {
            theNewName = (char *) malloc((strlen(theOldName) + 3) * sizeof(char));
            strcpy(theNewName, theOldName);
            strcat(theNewName, i == 0 ? ".a" : ".b");
            set_triangulation_name(thePieces[i], theNewName);
            free(theNewName);

            PyList_Append(thePieceList, Py_BuildValue("l", (long) thePieces[i]));
        }

    return thePieceList;
}


static char *RelationToString(int *aSnapPeaRelation)
{
    int     *theLetter,
            theRelationLength,
            i;
    char    *theCString;
    
    //  How long is aSnapPeaRelation?
    theRelationLength = 0;
    for (   theLetter = aSnapPeaRelation;
            *theLetter != 0;
            theLetter++)
        theRelationLength++;
    
    //  Allocate memory for theCString.
    theCString = malloc(theRelationLength + 1);
    assert(theCString != NULL);

    //  Write the relation.
    for (   theLetter = aSnapPeaRelation, i = 0;
            *theLetter != 0;
            theLetter++, i++)
    {
        if (*theLetter > 0)
            theCString[i] = 'a' - 1 + *theLetter;
        else
            theCString[i] = 'A' - 1 - *theLetter;
    }
    assert(i == theRelationLength);
    theCString[i] = 0;
    
    return theCString;
}

static PyObject *wrap_fundamental_group(PyObject *self, PyObject *args)
{
    long                theTriangulation;
    int                 theSimplifyFlag,
                        theFillingsAffectGeneratorsFlag,
                        theMinimizeNumGeneratorsFlag;
    GroupPresentation   *theFundamentalGroup;

    if (!PyArg_ParseTuple(  args,
                            "liii",
                            &theTriangulation,
                            &theSimplifyFlag,
                            &theFillingsAffectGeneratorsFlag,
                            &theMinimizeNumGeneratorsFlag))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_fundamental_group() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    theFundamentalGroup = fundamental_group(
                        (Triangulation *) theTriangulation,
                        theSimplifyFlag,
                        theFillingsAffectGeneratorsFlag,
                        theMinimizeNumGeneratorsFlag);

    return Py_BuildValue("l", (long) theFundamentalGroup);
}

static PyObject *wrap_free_group_presentation(PyObject *self, PyObject *args)
{
    long                    theFundamentalGroup;

    if (!PyArg_ParseTuple(args, "l", &theFundamentalGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_free_group_presentation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    free_group_presentation((GroupPresentation *) theFundamentalGroup);
    
    return Py_BuildValue("");
}

static PyObject *wrap_fg_get_num_generators(PyObject *self, PyObject *args)
{
    long    theFundamentalGroup;

    if (!PyArg_ParseTuple(args, "l", &theFundamentalGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_fg_get_num_generators() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    return Py_BuildValue("i",
        (fg_integer_fillings((GroupPresentation *) theFundamentalGroup) == TRUE) ?
        fg_get_num_generators((GroupPresentation *) theFundamentalGroup) :
        0);
}

static PyObject *wrap_fg_get_relations(PyObject *self, PyObject *args)
{
    long        theFundamentalGroup;
    int         theNumRelations,
                *theRawRelation,
                i;
    char        *theString;
    PyObject    *theRelationList;

    if (!PyArg_ParseTuple(args, "l", &theFundamentalGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_fg_get_relations() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    if (fg_integer_fillings((GroupPresentation *) theFundamentalGroup) == TRUE)
    {
        theNumRelations = fg_get_num_relations((GroupPresentation *) theFundamentalGroup);
        
        theRelationList = PyList_New(theNumRelations);

        for (i = 0; i < theNumRelations; i++)
        {
            theRawRelation = fg_get_relation((GroupPresentation *) theFundamentalGroup, i);
            theString = RelationToString(theRawRelation);
            fg_free_relation(theRawRelation);
            PyList_SetItem(theRelationList, i, Py_BuildValue("s", theString));
            free(theString);
        }
        
        return theRelationList;
    }
    else    
        return Py_BuildValue("[]");
}

static PyObject *wrap_fg_representation(PyObject *self, PyObject *args)
{
    long                    theFundamentalGroup;
    int                     theWordLength;
    char                    *theWord;
    int                     *theSnapPeaWord;
    int                     i;
    FuncResult              theError;
    O31Matrix               theO31Matrix;
    MoebiusTransformation   theMoebiusTransformation;
    double                  (*o)[4];
    Complex                 (*s)[2];

    if (!PyArg_ParseTuple(  args,
                            "ls",
                            &theFundamentalGroup,
                            &theWord))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_fg_representation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    theWordLength   = strlen(theWord);
    theSnapPeaWord  = malloc((theWordLength + 1) * sizeof(int));
    
    for (i = 0; i < theWordLength; i++)
    {
        if (islower(theWord[i]))
            theSnapPeaWord[i] =  (theWord[i] - 'a' + 1);
        else if (isupper(theWord[i]))
            theSnapPeaWord[i] = -(theWord[i] - 'A' + 1);
    }
    theSnapPeaWord[theWordLength] = 0;

    theError = fg_word_to_matrix(
                        (GroupPresentation *) theFundamentalGroup,
                        theSnapPeaWord,
                        theO31Matrix,
                        &theMoebiusTransformation);
    
    free(theSnapPeaWord);
    
    if (theError == func_bad_input)
    {
        PyErr_SetString(PyExc_ValueError, "letter in word is not a valid generator");
        return NULL;
    }

    o = theO31Matrix;
    s = theMoebiusTransformation.matrix;
    return Py_BuildValue
    (   "(((dddd)(dddd)(dddd)(dddd)) (i(((dd)(dd))((dd)(dd)))))",

        o[0][0], o[0][1], o[0][2], o[0][3], 
        o[1][0], o[1][1], o[1][2], o[1][3], 
        o[2][0], o[2][1], o[2][2], o[2][3], 
        o[3][0], o[3][1], o[3][2], o[3][3], 
        
        theMoebiusTransformation.parity == orientation_preserving ? 0 : 1,
        s[0][0].real, s[0][0].imag,    s[0][1].real, s[0][1].imag,
        s[1][0].real, s[1][0].imag,    s[1][1].real, s[1][1].imag
    );
}

static PyObject *wrap_fg_peripheral_curves(PyObject *self, PyObject *args)
{
    long                    theFundamentalGroup;
    int                     theNumCusps,
                            i;
    int                     *theMeridian,
                            *theLongitude;
    char                    *theMString,
                            *theLString;
    PyObject                *thePeripheralCurves;

    if (!PyArg_ParseTuple(args, "l", &theFundamentalGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_fg_peripheral_curves() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    theNumCusps = fg_get_num_cusps((GroupPresentation *) theFundamentalGroup);

    thePeripheralCurves = PyList_New(theNumCusps);
    
    for (i = 0; i < theNumCusps; i++)
    {
        theMeridian  = fg_get_meridian ((GroupPresentation *) theFundamentalGroup, i);
        theLongitude = fg_get_longitude((GroupPresentation *) theFundamentalGroup, i);
        
        theMString = RelationToString(theMeridian );
        theLString = RelationToString(theLongitude);
        
        fg_free_relation(theMeridian);
        fg_free_relation(theLongitude);
        
        PyList_SetItem(thePeripheralCurves, i,
            Py_BuildValue("[ss]", theMString, theLString));
        
        free(theMString);
        free(theLString);
    }
    
    return thePeripheralCurves;
}


static PyObject *wrap_core_geodesic(PyObject *self, PyObject *args)
{
    long        theTriangulation;
    int         theCuspIndex;
    int         theSingularityIndex;
    Complex     theComplexLength,
                theHolonomy,
                theTraceSquared,
                theTrace,
                theEigenvalue;
    int         thePrecision;
    PyObject    *theDictionary,
                *theObject;

    static Complex  zero = { 0.0, 0.0},
                    one  = { 1.0, 0.0},
                    two  = { 2.0, 0.0};

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theCuspIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_core_geodesic() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    core_geodesic(  (Triangulation *) theTriangulation,
                    theCuspIndex,
                    &theSingularityIndex,
                    &theComplexLength,
                    &thePrecision);
    
    if (theSingularityIndex > 0)
    {
        //  complex_length.c explains the relationship between
        //  holonomy, complex length and trace.  It uses the
        //  variable 'k' for theHolonomy.
        theHolonomy     = complex_exp(theComplexLength);
        theTraceSquared = complex_plus (
                            complex_plus(
                                theHolonomy,
                                complex_div(one, theHolonomy)),
                            two);
        theTrace        = complex_sqrt(theTraceSquared);
        theEigenvalue   = complex_exp(complex_real_mult(0.5, theComplexLength));
    }
    else
    {
        theComplexLength    = zero;
        theHolonomy         = zero;
        theTraceSquared     = zero;
        theTrace            = zero;
        theEigenvalue       = zero;
    }

    theDictionary = PyDict_New();

    theObject = Py_BuildValue("i", theSingularityIndex);
    PyDict_SetItemString(theDictionary, "singularity index", theObject);
    Py_DECREF(theObject);

    theObject = Py_BuildValue("i", thePrecision);
    PyDict_SetItemString(theDictionary, "precision", theObject);
    Py_DECREF(theObject);

    theObject = PyComplex_FromDoubles(theComplexLength.real, theComplexLength.imag  );
    PyDict_SetItemString(theDictionary, "complex length", theObject);
    Py_DECREF(theObject);

    theObject = PyComplex_FromDoubles(theHolonomy.real,      theHolonomy.imag    );
    PyDict_SetItemString(theDictionary, "holonomy", theObject);
    Py_DECREF(theObject);

    theObject = PyComplex_FromDoubles(theTraceSquared.real,  theTraceSquared.imag);
    PyDict_SetItemString(theDictionary, "trace squared", theObject);
    Py_DECREF(theObject);

    theObject = PyComplex_FromDoubles(theTrace.real,         theTrace.imag       );
    PyDict_SetItemString(theDictionary, "trace", theObject);
    Py_DECREF(theObject);

    theObject = PyComplex_FromDoubles(theEigenvalue.real,    theEigenvalue.imag  );
    PyDict_SetItemString(theDictionary, "eigenvalue", theObject);
    Py_DECREF(theObject);
    
    return theDictionary;
}


static PyObject *wrap_shortest_curves_become_meridians(PyObject *self, PyObject *args)
{
    long        theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_shortest_curves_become_meridians() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    install_shortest_bases((Triangulation *)theTriangulation);

    return Py_BuildValue("");
}

static PyObject *wrap_current_fillings_become_meridians(PyObject *self, PyObject *args)
{
    long        theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_current_fillings_become_meridians() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    install_current_curve_bases((Triangulation *)theTriangulation);

    return Py_BuildValue("");
}


static PyObject *wrap_basic_simplification(PyObject *self, PyObject *args)
{
    long        theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_basic_simplification() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    basic_simplification((Triangulation *)theTriangulation);

    return Py_BuildValue("");
}

static PyObject *wrap_randomize_triangulation(PyObject *self, PyObject *args)
{
    long        theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_randomize_triangulation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    randomize_triangulation((Triangulation *)theTriangulation);

    return Py_BuildValue("");
}

static PyObject *wrap_reorient(PyObject *self, PyObject *args)
{
    long        theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_reorient() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    reorient((Triangulation *)theTriangulation);

    return Py_BuildValue("");
}

static PyObject *wrap_proto_canonize(PyObject *self, PyObject *args)
{
    long        theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_proto_canonize() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    switch (proto_canonize((Triangulation *)theTriangulation))
    {
        case func_OK:
            if (is_canonical_triangulation((Triangulation *)theTriangulation) == FALSE)
                uAcknowledge("The canonical cell decomposition contains cells other than tetrahedra.  Such cells have been arbitrarily subdivided into tetrahedra.");
            break;
        
        case func_failed:
            uAcknowledge("Only hyperbolic manifolds have canonical decompositions.");
            break;
        
        default:
            //  Should not occur.
            PyErr_SetString(PyExc_RuntimeError, "Unknown error in proto_canonize().");
            return NULL;
    }

    return Py_BuildValue("");
}

static PyObject *wrap_is_canonical_triangulation(PyObject *self, PyObject *args)
{
    long        theTriangulation;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_is_canonical_triangulation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    return Py_BuildValue("i", is_canonical_triangulation((Triangulation *)theTriangulation));
}

static PyObject *wrap_symmetry_group(PyObject *self, PyObject *args)
{
    long            theTriangulation;
    SymmetryGroup   *symmetry_group_of_manifold,
                    *symmetry_group_of_link;
    Triangulation   *symmetric_triangulation;
    Boolean         is_full_group;

    if (!PyArg_ParseTuple(args, "l", &theTriangulation))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    symmetry_group_of_manifold  = NULL;
    symmetry_group_of_link      = NULL;
    symmetric_triangulation     = NULL;
    is_full_group               = FALSE;

    //  In the rare event that compute_symmetry_group() fails,
    //  the various pointers will be left as NULL.
    
    (void) compute_symmetry_group(  (Triangulation *)theTriangulation,
                                    &symmetry_group_of_manifold,
                                    &symmetry_group_of_link,
                                    &symmetric_triangulation,
                                    &is_full_group);

    return Py_BuildValue("[llli]",
                            (long) symmetry_group_of_manifold,
                            (long) symmetry_group_of_link,
                            (long) symmetric_triangulation,
                            is_full_group);
}

static PyObject *wrap_free_symmetry_group(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_free_symmetry_group() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    free_symmetry_group((SymmetryGroup *) theSymmetryGroup);
    
    return Py_BuildValue("");
}

static PyObject *wrap_symmetry_group_order(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_order() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    return Py_BuildValue("i", symmetry_group_order((SymmetryGroup *) theSymmetryGroup));
}

static PyObject *wrap_symmetry_group_is_abelian(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_is_abelian() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    return Py_BuildValue("i", (int) symmetry_group_is_abelian((SymmetryGroup *) theSymmetryGroup, NULL));
}

static PyObject *wrap_symmetry_group_abelian_description(PyObject *self, PyObject *args)
{
    long            theSymmetryGroup;
    AbelianGroup    *theAbelianDescription;
    PyObject        *theList;
    int             i;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_abelian_description() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    if (symmetry_group_is_abelian(  (SymmetryGroup *) theSymmetryGroup,
                                    &theAbelianDescription)
        != TRUE)
    {
        PyErr_SetString(PyExc_TypeError, "Abelian desciptions aren't available for nonabelian groups.  Use SymmetryGroup.is_abelian() to test.");
        return NULL;
    }

    theList = PyList_New(theAbelianDescription->num_torsion_coefficients);
    for (i = 0; i < theAbelianDescription->num_torsion_coefficients; i++)
        PyList_SetItem(theList, i, Py_BuildValue("l", theAbelianDescription->torsion_coefficients[i]));
    
    return theList;
}

static PyObject *wrap_symmetry_group_is_dihedral(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_is_dihedral() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    return Py_BuildValue("i", (int) symmetry_group_is_dihedral((SymmetryGroup *) theSymmetryGroup));
}

static PyObject *wrap_symmetry_group_is_polyhedral(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_is_polyhedral() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    return Py_BuildValue("i", (int) symmetry_group_is_polyhedral(
        (SymmetryGroup *) theSymmetryGroup, NULL, NULL, NULL, NULL));
}

static PyObject *wrap_symmetry_group_polyhedral_description(PyObject *self, PyObject *args)
{
    long        theSymmetryGroup;
    Boolean     is_binary_group;
    int         p,
                q,
                r;
    PyObject    *theDictionary,
                *theObject;
    char        theName[64];

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_polyhedral_description() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    if (symmetry_group_is_polyhedral(   (SymmetryGroup *) theSymmetryGroup,
                                        &is_binary_group,
                                        &p,
                                        &q,
                                        &r
        ) != TRUE)
    {
        PyErr_SetString(PyExc_ValueError, "Symmetry group is not a polyhedral group, and therefore does not have a polyhedral description.  Use SymmetryGroup.is_polyhedral() to test.");
        return NULL;
    }

    //  Use the fact that the SnapPea kernel generates the (p,q,r)
    //  in ascending order.
    
    assert(p == 2);
    
    switch(q)
    {
        case 2:
            //  Plain dihedral groups are best handled
            //  in the specialized dihedral group code,
            //  which puts the elements into a natural order.
            assert(is_binary_group == TRUE);
            
            sprintf(theName, "binary dihedral group <2,2,%d>", r);

            break;
    
        case 3:

            strcpy(theName, is_binary_group ? "binary " : "");
    
            switch (r)
            {
                case 3:  strcat(theName, "tetrahedral group");  break;
                case 4:  strcat(theName,  "octahedral group");  break;
                case 5:  strcat(theName, "icosahedral group");  break;
                default:  assert(FALSE);
            }
            break;
        
        default:  assert(FALSE);
    }

    theDictionary = PyDict_New();
    
    theObject = Py_BuildValue("i", is_binary_group ? 1 : 0);
    PyDict_SetItemString(   theDictionary,  "is binary",    theObject);
    Py_DECREF(theObject);

    theObject = Py_BuildValue("(iii)", p, q, r);
    PyDict_SetItemString(   theDictionary,  "pqr",          theObject);
    Py_DECREF(theObject);

    theObject = Py_BuildValue("s", theName);
    PyDict_SetItemString(   theDictionary,  "name",         theObject);
    Py_DECREF(theObject);
    
    return theDictionary;
}

static PyObject *wrap_symmetry_group_is_S5(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_is_S5() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    return Py_BuildValue("i", (int) symmetry_group_is_S5((SymmetryGroup *) theSymmetryGroup));
}

static PyObject *wrap_symmetry_group_is_direct_product(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_is_direct_product() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    return Py_BuildValue("i", (int) symmetry_group_is_direct_product((SymmetryGroup *) theSymmetryGroup));
}

static PyObject *wrap_symmetry_group_factor(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;
    int     theFactorIndex;

    if (!PyArg_ParseTuple(args, "li", &theSymmetryGroup, &theFactorIndex))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_factor() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    if (theFactorIndex != 0 && theFactorIndex != 1)
    {
        PyErr_SetString(PyExc_ValueError, "wrap_symmetry_group_factor() provides only two factors, indexed 0 and 1, which may be factored recursively if necessary.");
        return NULL;
    }
        
    return Py_BuildValue("l", (long int) get_symmetry_group_factor((SymmetryGroup *) theSymmetryGroup, theFactorIndex));
}

static PyObject *wrap_symmetry_group_is_amphicheiral(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_is_amphicheiral() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    return Py_BuildValue("i", symmetry_group_is_amphicheiral((SymmetryGroup *) theSymmetryGroup));
}

static PyObject *wrap_symmetry_group_invertible_knot(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_invertible_knot() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    return Py_BuildValue("i", symmetry_group_invertible_knot((SymmetryGroup *) theSymmetryGroup));
}


static PyObject *wrap_symmetry_group_commutator_subgroup(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_commutator_subgroup() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
        
    return Py_BuildValue("l", (long int) get_commutator_subgroup((SymmetryGroup *) theSymmetryGroup));
}

static PyObject *wrap_symmetry_group_abelianization(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_abelianization() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    return Py_BuildValue("l", (long int) get_abelianization((SymmetryGroup *) theSymmetryGroup));
}

static PyObject *wrap_symmetry_group_center(PyObject *self, PyObject *args)
{
    long    theSymmetryGroup;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_center() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    return Py_BuildValue("l", (long int) get_center((SymmetryGroup *) theSymmetryGroup));
}

static PyObject *wrap_symmetry_group_presentation(PyObject *self, PyObject *args)
{
    long                        theSymmetryGroup;
    SymmetryGroupPresentation   *thePresentation;
    PyObject                    *theRelations,
                                *theRelation,
                                *theDictionary,
                                *theObject;
    int                         theNumGenerators,
                                theNumRelations,
                                theNumFactors,
                                theGenerator,
                                thePower,
                                i,
                                j;

    if (!PyArg_ParseTuple(args, "l", &theSymmetryGroup))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_symmetry_group_presentation() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    thePresentation = get_symmetry_group_presentation((SymmetryGroup *) theSymmetryGroup);

    theNumGenerators = sg_get_num_generators(thePresentation);
    theNumRelations  = sg_get_num_relations(thePresentation);
    
    theRelations = PyList_New(theNumRelations);
    for (i = 0; i < theNumRelations; i++)
    {
        theNumFactors = sg_get_num_factors(thePresentation, i);
        theRelation = PyList_New(theNumFactors);
        for (j = 0; j < theNumFactors; j++)
        {
            sg_get_factor(thePresentation, i, j, &theGenerator, &thePower);
            PyList_SetItem(theRelation, j, Py_BuildValue("(ii)", theGenerator, thePower));
        }
        PyList_SetItem(theRelations, i, theRelation);
    }

    free_symmetry_group_presentation(thePresentation);

    theDictionary = PyDict_New();

    theObject = Py_BuildValue("i", theNumGenerators);
    PyDict_SetItemString(theDictionary, "number of generators", theObject);
    Py_DECREF(theObject);

    theObject = Py_BuildValue("i", theNumRelations);
    PyDict_SetItemString(theDictionary, "number of relations",  theObject);
    Py_DECREF(theObject);

    PyDict_SetItemString(theDictionary, "relations",            theRelations);
    Py_DECREF(theRelations);
    
    return theDictionary;
}


static PyObject *wrap_tet_shapes(PyObject *self, PyObject *args)
{
    long        theTriangulation;
    int         theFixedAlignment;
    int         theNumTetrahedra,
                i;
    PyObject    *theList,
                *theDictionary,
                *theObject;
    double      theShapeRectReal,
                theShapeRectImag,
                theShapeLogReal,
                theShapeLogImag;
    int         thePrecisionRectReal,
                thePrecisionRectImag,
                thePrecisionLogReal,
                thePrecisionLogImag;
    Boolean     theGeometricFlag;

    if (!PyArg_ParseTuple(args, "li", &theTriangulation, &theFixedAlignment))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_tet_shapes() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    theNumTetrahedra = get_num_tetrahedra((Triangulation *)theTriangulation);

    theList = PyList_New(theNumTetrahedra);

    for (i = 0; i < theNumTetrahedra; i++)
    {
        get_tet_shape(  (Triangulation *)theTriangulation,
                        i,
                        theFixedAlignment,
                        &theShapeRectReal,
                        &theShapeRectImag,
                        &theShapeLogReal,
                        &theShapeLogImag,
                        &thePrecisionRectReal,
                        &thePrecisionRectImag,
                        &thePrecisionLogReal,
                        &thePrecisionLogImag,
                        &theGeometricFlag);
        
        theDictionary = PyDict_New();

        theObject = PyComplex_FromDoubles(theShapeRectReal, theShapeRectImag);
        PyDict_SetItemString(theDictionary, "shape rect", theObject);
        Py_DECREF(theObject);

        theObject = PyComplex_FromDoubles(theShapeLogReal,  theShapeLogImag );
        PyDict_SetItemString(theDictionary, "shape log",  theObject);
        Py_DECREF(theObject);

        theObject = Py_BuildValue("i", thePrecisionRectReal);
        PyDict_SetItemString(theDictionary, "precision rect real", theObject);
        Py_DECREF(theObject);

        theObject = Py_BuildValue("i", thePrecisionRectImag);
        PyDict_SetItemString(theDictionary, "precision rect imag", theObject);
        Py_DECREF(theObject);

        theObject = Py_BuildValue("i", thePrecisionLogReal);
        PyDict_SetItemString(theDictionary, "precision log real",  theObject);
        Py_DECREF(theObject);

        theObject = Py_BuildValue("i", thePrecisionLogImag);
        PyDict_SetItemString(theDictionary, "precision log imag",  theObject);
        Py_DECREF(theObject);

        theObject = Py_BuildValue("i", theGeometricFlag);
        PyDict_SetItemString(theDictionary, "is geometric",  theObject);
        Py_DECREF(theObject);
        
        PyList_SetItem(theList, i, theDictionary);
    }
    
    return theList;
}


static PyObject *wrap_Dirichlet(PyObject *self, PyObject *args)
{
    long            theTriangulation;
    WEPolyhedron    *theDirichletDomain;
    int             theCentroidAtOriginFlag,
                    theMaximizeInjRadiusFlag;
    double          displacement[3];

    if (!PyArg_ParseTuple(args, "liiddd",
            &theTriangulation,
            &theCentroidAtOriginFlag,
            &theMaximizeInjRadiusFlag,
            &displacement[0],
            &displacement[1],
            &displacement[2]))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_Dirichlet() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }

    theDirichletDomain = Dirichlet_with_displacement(
        (Triangulation *) theTriangulation, //  manifold
        displacement,                       //  displacement
        1e-8,                               //  vertex epsilon
        theCentroidAtOriginFlag,            //  centroid_at_origin
        Dirichlet_keep_going,               //  DirichletInteractivity
        theMaximizeInjRadiusFlag);          //  maximize_injectivity_radius
    
    return Py_BuildValue("l", (long) theDirichletDomain);
}

static PyObject *wrap_free_Dirichlet_domain(PyObject *self, PyObject *args)
{
    long    theDirichletDomain;

    if (!PyArg_ParseTuple(args, "l", &theDirichletDomain))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_free_Dirichlet_domain() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    
    if (theDirichletDomain != 0)
        free_Dirichlet_domain((WEPolyhedron *) theDirichletDomain);
    
    return Py_BuildValue("");
}

static PyObject *wrap_Dirichlet_num_vertices(PyObject *self, PyObject *args)
{
    long            theDirichletDomain;
    WEPolyhedron    *theDD;

    if (!PyArg_ParseTuple(args, "l", &theDirichletDomain))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_Dirichlet_num_vertices() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    theDD = (WEPolyhedron *) theDirichletDomain;
    
    if (theDD != NULL)
        return Py_BuildValue("i", theDD->num_vertices);
    else
        return Py_BuildValue("i", 0);
}

static PyObject *wrap_Dirichlet_num_edges(PyObject *self, PyObject *args)
{
    long            theDirichletDomain;
    WEPolyhedron    *theDD;

    if (!PyArg_ParseTuple(args, "l", &theDirichletDomain))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_Dirichlet_num_edges() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    theDD = (WEPolyhedron *) theDirichletDomain;
    
    if (theDD != NULL)
        return Py_BuildValue("i", theDD->num_edges);
    else
        return Py_BuildValue("i", 0);
}

static PyObject *wrap_Dirichlet_num_faces(PyObject *self, PyObject *args)
{
    long            theDirichletDomain;
    WEPolyhedron    *theDD;

    if (!PyArg_ParseTuple(args, "l", &theDirichletDomain))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_Dirichlet_num_faces() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    theDD = (WEPolyhedron *) theDirichletDomain;
    
    if (theDD != NULL)
        return Py_BuildValue("i", theDD->num_faces);
    else
        return Py_BuildValue("i", 0);
}

static PyObject *wrap_Dirichlet_vertices(PyObject *self, PyObject *args)
{
    long            theDirichletDomain;
    WEPolyhedron    *theDD;
    PyObject        *theList;
    WEVertex        *theVertex;
    int             theCount;

    if (!PyArg_ParseTuple(args, "l", &theDirichletDomain))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_Dirichlet_vertices() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    theDD = (WEPolyhedron *) theDirichletDomain;

    if (theDD == NULL)
        return Py_BuildValue("[]");

    theList = PyList_New(theDD->num_vertices);

    for (theVertex = theDD->vertex_list_begin.next, theCount = 0;
         theVertex != &theDD->vertex_list_end;
         theVertex = theVertex->next, theCount++)
    {
        PyList_SetItem(theList, theCount, 
            Py_BuildValue("(fff)", theVertex->x[1], theVertex->x[2], theVertex->x[3]));
    }
    assert(theCount == theDD->num_vertices);
    
    return theList;
}

static PyObject *wrap_Dirichlet_faces(PyObject *self, PyObject *args)
{
    long            theDirichletDomain;
    WEPolyhedron    *theDD;
    WEVertex        *theVertex,
                    *theCandidateVertex;
    WEEdge          *theEdge;
    WEFace          *theFace;
    PyObject        *theList,
                    *theVertexList;
    int             theFaceCount,
                    theVertexCount,
                    theVertexIndex;

    if (!PyArg_ParseTuple(args, "l", &theDirichletDomain))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_Dirichlet_faces() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    theDD = (WEPolyhedron *) theDirichletDomain;

    if (theDD == NULL)
        return Py_BuildValue("[]");

    theList = PyList_New(theDD->num_faces);

    for (theFace = theDD->face_list_begin.next, theFaceCount = 0;
         theFace != &theDD->face_list_end;
         theFace = theFace->next, theFaceCount++)
    {
        theVertexList   = PyList_New(theFace->num_sides);
        theVertexCount  = 0;

        //  Enumerate the indices of the vertices, travelling counterclockwise.
        theEdge = theFace->some_edge;
        do
        {
            //  Find the more clockwise of theEdge's two endpoints.
            if (theEdge->f[left] == theFace)
                theVertex = theEdge->v[tip];
            else
                theVertex = theEdge->v[tail];

            //  The vertices don't have explicit indices, so we have
            //  to figure out on the fly where each vertex is on the
            //  list.  This isn't very efficient, but this isn't
            //  a time-critical routine.
            theVertexIndex = 0;
            theCandidateVertex = theDD->vertex_list_begin.next;
            while (theCandidateVertex != theVertex)
            {
                theVertexIndex++;
                theCandidateVertex = theCandidateVertex->next;
            }
            PyList_SetItem( theVertexList,
                            theVertexCount++,
                            Py_BuildValue("i", theVertexIndex));
    
            //  Move on to the next edge.
            if (theEdge->f[left] == theFace)
                theEdge = theEdge->e[tip][left];
            else
                theEdge = theEdge->e[tail][right];
    
        } while (theEdge != theFace->some_edge);
        assert(theVertexCount == theFace->num_sides);
        
        PyList_SetItem(theList, theFaceCount, theVertexList);
    }
    assert(theFaceCount == theDD->num_faces);

    return theList;
}

static PyObject *wrap_Dirichlet_face_colors(PyObject *self, PyObject *args)
{
    long            theDirichletDomain;
    WEPolyhedron    *theDD;
    PyObject        *theList;
    WEFace          *theFace;
    int             theCount;
    double          h,
                    r,
                    g,
                    b;

    if (!PyArg_ParseTuple(args, "l", &theDirichletDomain))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_Dirichlet_face_colors() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    theDD = (WEPolyhedron *) theDirichletDomain;

    if (theDD == NULL)
        return Py_BuildValue("[]");

    theList = PyList_New(theDD->num_faces);

    for (theFace = theDD->face_list_begin.next, theCount = 0;
         theFace != &theDD->face_list_end;
         theFace = theFace->next, theCount++)
    {
        //  Determine the face color.
        
        //  Uh-oh.  I know I want
        //
        //      hue         = theFace->f_class->hue;
        //      saturation  = 1.0;
        //      value       = 0.75;
        //
        //  but I have no idea how to convert HSV to RGB on unix.
        //  For now, let's just fudge it.
        
        h = 3.0 * theFace->f_class->hue;
        if (h < 1.0)
        {
            r = 0.25 + 0.75*(1.0 - h);
            g = 0.25 + 0.75*(h - 0.0);
            b = 0.25;
        }
        else if (h < 2.0)
        {
            r = 0.25;
            g = 0.25 + 0.75*(2.0 - h);
            b = 0.25 + 0.75*(h - 1.0);
        }
        else
        {
            r = 0.25 + 0.75*(h - 2.0);
            g = 0.25;
            b = 0.25 + 0.75*(3.0 - h);
        }

        PyList_SetItem(theList, theCount,
            Py_BuildValue("(fff)", r, g, b));
    }
    assert(theCount == theDD->num_faces);
    
    return theList;
}

static PyObject *wrap_Dirichlet_face_pairings(PyObject *self, PyObject *args)
{
    long            theDirichletDomain;
    WEPolyhedron    *theDD;
    PyObject        *theList;
    WEFace          *theFace;
    int             theCount;

    if (!PyArg_ParseTuple(args, "l", &theDirichletDomain))
    {
        PyErr_SetString(PyExc_TypeError, "wrap_Dirichlet_face_pairings() in SnapPeaC.c received data of the wrong type.");
        return NULL;
    }
    theDD = (WEPolyhedron *) theDirichletDomain;

    if (theDD == NULL)
        return Py_BuildValue("[]");

    theList = PyList_New(theDD->num_faces);

    for (theFace = theDD->face_list_begin.next, theCount = 0;
         theFace != &theDD->face_list_end;
         theFace = theFace->next, theCount++)
    {
        PyList_SetItem(theList, theCount,
            Py_BuildValue("(ffff)(ffff)(ffff)(ffff)",
            (*(theFace->group_element))[0][0],
            (*(theFace->group_element))[0][1],
            (*(theFace->group_element))[0][2],
            (*(theFace->group_element))[0][3],
            (*(theFace->group_element))[1][0],
            (*(theFace->group_element))[1][1],
            (*(theFace->group_element))[1][2],
            (*(theFace->group_element))[1][3],
            (*(theFace->group_element))[2][0],
            (*(theFace->group_element))[2][1],
            (*(theFace->group_element))[2][2],
            (*(theFace->group_element))[2][3],
            (*(theFace->group_element))[3][0],
            (*(theFace->group_element))[3][1],
            (*(theFace->group_element))[3][2],
            (*(theFace->group_element))[3][3]
            ));
    }
    assert(theCount == theDD->num_faces);
    
    return theList;
}
