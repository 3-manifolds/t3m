/*
 *  unix_cusped_census_main.c
 *
 *  The main() in this file illustrates the use of the unix-style
 *  GetCuspedCensusManifold().  It reads all the manifolds of
 *  5 or fewer tetrahedra and prints their volumes.
 */

#include "SnapPea.h"
#include "unix_cusped_census.h"
#include <stdio.h>
#include <stdlib.h>


int main(void)
{
    int             theIndex;
    Triangulation   *theTriangulation;

    for (theIndex = 0; theIndex < gNumOrientableCuspedCensusMflds[5]; theIndex++)
    {
        theTriangulation = GetCuspedCensusManifold(
                    5,
                    oriented_manifold /* ignored for 5-tet census */,
                    theIndex);

        if (theTriangulation != NULL)
        {
            printf("%s  %lf\n", get_triangulation_name(theTriangulation), volume(theTriangulation, NULL));
            free_triangulation(theTriangulation);
        }
        else
            printf("Couldn't read census manifold.\n");
    }
    
    return 0;
}
