//  detab.c
//
//  SnapPea's source code looks good when tab stops are
//  set to 4 spaces.  Unfortunately the unix tradition
//  is 8 space tab stops.  So to view SnapPea code on
//  a unix machine, you should either
//
//  (1) set your editor to use 4 space tab stops, or,
//      if that's not convenient,
//
//  (2) use this program to replace the tab stops with
//      spaces.  To do so, just type
//
//          gcc -o detab detab.c
//          chmod 644 */*.c */*.h
//          detab */*.c */*.h
//          chmod 444 */*.c */*.h

#include <stdio.h>

#define TAB_SIZE        4
#define SCRATCH_FILE    "/tmp/detab.temp"


int main(int argc, char *argv[])
{
    int     i;
    FILE    *fp1,
            *fp2;
    int     pos,    // position in line, modulo TAB_SIZE
            j;
    int     ch;
    char    scratch[128];
    
    for (i = 1; i < argc; i++)
    {
        fp1 = fopen(argv[i], "r");
        if (fp1 == NULL)
        {
            printf("couldn't open %s\n", argv[i]);
            continue;
        }
        fp2 = fopen(SCRATCH_FILE, "w");
        if (fp2 == NULL)
        {
            printf("couldn't open %s\n", SCRATCH_FILE);
            fclose(fp1);
            continue;
        }
        
        pos = 0;
        while ( (ch = getc(fp1)) != EOF )
            switch (ch)
            {
                case '\t':
                    for (j = TAB_SIZE - pos%TAB_SIZE; --j >= 0; )
                        putc(' ', fp2);
                    pos = 0;
                    break;

                case '\n':
                    putc(ch, fp2);
                    pos = 0;
                    break;

                default:
                    putc(ch, fp2);
                    pos++;
                    break;
            }
        
        fclose(fp1);
        fclose(fp2);

        fp1 = fopen(argv[i], "w");
        if (fp1 == NULL)
        {
            printf("couldn't open %s\n", argv[i]);
            continue;
        }
        fp2 = fopen(SCRATCH_FILE, "r");
        if (fp2 == NULL)
        {
            printf("couldn't open %s\n", SCRATCH_FILE);
            fclose(fp1);
            continue;
        }
        
        while ((ch = getc(fp2)) != EOF)
            putc(ch, fp1);
        
        fclose(fp1);
        fclose(fp2);
    }

    strcpy(scratch, "rm -f ");  
    strcat(scratch, SCRATCH_FILE);  
    system(scratch);
    return 0;
}
