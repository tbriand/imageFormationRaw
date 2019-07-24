#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "xfopen.c"
#include "parsenumbers.c"

int main(int c, char *v[])
{
    if (c < 3) {
            fprintf(stderr,"usage:\n\t%s filename opt opt2\n", *v);
            //                         0 1        2   3
            return EXIT_FAILURE;
    }
    char *filename = c > 1 ? v[1] : "-";
    int opt  = c > 2 ? atoi(v[2]) : 0;
    int opt2 = c > 3 ? atoi(v[3]) : 1;

    // read file
    FILE *f = xfopen(filename, "r");
    int nvalues;
    double *values = read_ascii_doubles(f, &nvalues);
    xfclose(f);

    // first computation of the mean
    double mean = 0;
    for (int i = 0; i < nvalues; i++) mean += values[i];
    mean /= nvalues;

    if ( opt != 1 ) {
        if ( opt2 )
            printf("%.5lf\n", mean);
        else
            printf("%lf\n", mean);
        if ( opt == 0)
            return EXIT_SUCCESS;
    }

    // first computation of the std
    double std = 0;
    for (int i = 0; i < nvalues; i++) std += (values[i] - mean)*(values[i] - mean);
    std = sqrt(std/nvalues);

    if ( opt > 0 ) {
        if ( opt2 )
            printf("%.5lf\n", std);
        else
            printf("%lf\n", std);
        if ( opt == 1)
            return EXIT_SUCCESS;
    }

    // free memory
    free(values);

    return EXIT_SUCCESS;
}
