/*
 *  unix_sample_main.c
 *
 *  This file provides a sample main() illustrating how a unix-style
 *  program can read individual manifolds from files and do computations
 *  with them.  This sample main just computes the volume, but you could
 *  of course compute anything else available in the SnapPea kernel.
 *
 *  unix_cusped_census_main.c unix_closed_census_main.c provide
 *  unix-style code for looping through the census manifolds,
 *  so your program can do computations in batch mode.
 *  The new Python SnapPea is even more convenient for such purposes.
 */

#include <string.h>
#include <stdio.h>
#include "SnapPea.h"
#include "unix_file_io.h"

#define MANIFOLD_DIRECTORY          ""
#define DIRECTORY_NAME_LENGTH       0
#define MAX_MANIFOLD_NAME_LENGTH    50

static Triangulation    *get_manifold(void);


int main(void)
{
    Triangulation   *manifold;
    double          vol;
    int             digits;

    while ((manifold = get_manifold()) != NULL)
    {
        /*
         *  Here's where you put your code to do whatever
         *  you want with the manifold.
         */
        vol = volume(manifold, &digits);
        printf("volume is %.*lf\n", digits, vol);

        /*
         *  Free the manifold, and check for memory leaks.
         */
        free_triangulation(manifold);
        verify_my_malloc_usage();
    }
    
    return 0;
}


static Triangulation *get_manifold()
{
    char            manifold_name[MAX_MANIFOLD_NAME_LENGTH],
                    path_name[DIRECTORY_NAME_LENGTH + MAX_MANIFOLD_NAME_LENGTH];
    Triangulation   *manifold;

    do
    {
        printf("\nmanifold name: ");
        fgets(manifold_name, MAX_MANIFOLD_NAME_LENGTH, stdin);
        manifold_name[strlen(manifold_name) - 1] = 0;   /* overwrite the newline */
        if (manifold_name[0] == 0)
            return NULL;
    
        strcpy(path_name, MANIFOLD_DIRECTORY);
        strcat(path_name, manifold_name);
    
        manifold = get_triangulation(path_name);
        
    } while (manifold == NULL);
    
    return manifold;
}
