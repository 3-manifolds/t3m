/*
 *  unix_closed_census_main.c
 *
 *  Reads input from stdin.
 *  Each line of the input describes one closed manifold.
 *  Here's a sample input file:

0.9427 5 003 -3 1
0.9813 5 003 -2 3
1.0149 5 007  3 1

 *  The first entry in each row is the manifold's volume
 *          (for the benefit of human readers, not for the program itself).
 *  The second entry says which SnapPea census of cusped manifolds
 *          ( <=5 tetrahedra, 6 tetrahedra, or 7 tetrahedra ).
 *  The third entry is the index of the given cusped manifold within the census.
 *  The last two entries are the Dehn filling coefficients (m,l).
 *
 *  Note:  I have a huge file giving the data for vast numbers of
 *  low-volume closed orientable hyperbolic 3-manifolds, and a more
 *  modest collection of nonorientable manifolds.
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "SnapPea.h"
#include "unix_cusped_census.h"

#define FALSE   0
#define TRUE    1

static void do_something(Triangulation *aTriangulation);


int main(void)
{
    /*
     *  main() just iterates through the manifolds.
     *  The real action is in do_something() below.
     */
    
    double          theVolume;
    int             theCensus,
                    theIndex,
                    m,
                    l;
    Triangulation   *theTriangulation;

    while (5 == scanf(  "%lf%d%d%d%d",
                        &theVolume,
                        &theCensus,
                        &theIndex,
                        &m,
                        &l))
    {
        theTriangulation = GetCuspedCensusManifold(
                            theCensus,
                            oriented_manifold /* ignored for 5-tet census */,
                            theIndex);

        if (theTriangulation != NULL)
        {
            if (m != 0 || l != 0)
                set_cusp_info(theTriangulation, 0, FALSE, m, l);
            else
                set_cusp_info(theTriangulation, 0, TRUE, m, l);

            switch (do_Dehn_filling(theTriangulation))
            {
                case geometric_solution:
                case nongeometric_solution:
                    do_something(theTriangulation);
                    break;

                default:
                    printf("Couldn't find hyperbolic structure for manifold %d %d(%d,%d).\n",
                            theCensus, theIndex, m, l);
                    printf("Expected volume was %lf.\n", theVolume);
                    break;
            }
            free_triangulation(theTriangulation);
        }
        else
                printf("Couldn't read census manifold %d %d.\n", theCensus, theIndex);

        /* check for memory leaks */
        verify_my_malloc_usage();
    }

    return 0;
}


static void do_something(
    Triangulation   *aTriangulation)
{
    /*
     *  Your code goes here, to do whatever
     *  you want with each manifold.
     */

    printf( "volume is %14.10lf\n", volume(aTriangulation, NULL));
}
