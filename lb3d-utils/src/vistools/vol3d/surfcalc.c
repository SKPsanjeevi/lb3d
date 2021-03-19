#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vol3df.h"
#include "isosurface.h"

int main(int argc, char *argv[])
{
    char *volFilename=NULL,*dataFilename=NULL;
    struct vol3df_file_hints volHints = vol3df_default_file_hints;
    struct vol3df_file_hints dataHints = vol3df_default_file_hints;
    struct vol3d *surfVol=NULL,*dataVol=NULL;
    struct mcubes_state *mc=NULL;
    float sumScalarArea,sumScalar2Area,sumArea;
    float meanScalar,mean2Scalar;
    float variance;
    float isoVal=0.0;

    if (4!=argc) {
            fprintf(stderr,"%s volfile isoval datafile\n",argv[0]);
            return -1;
    }
    volHints.filename = volFilename=argv[1]; 
    if (NULL==(surfVol=vol3df_new_from_file_heuristic(volFilename,&volHints)))
    {
            fprintf(stderr,"Failed to read input volume file\n");
            return -1;
    }

    isoVal = atof(argv[2]);

    dataHints.filename = dataFilename=argv[3]; 
    if (NULL==(dataVol=vol3df_new_from_file_heuristic(dataFilename,&dataHints)))
    {
            fprintf(stderr,"Failed to read input volume file\n");
            return -1;
    }

    if (NULL==(mc=mcubes_new(surfVol,isoVal,dataVol,1))) {
            fprintf(stderr,"mcubes_new failed.\n");
            return -1;
    }

    mcubes_count_scalar(mc);

    if (0!=mcubes_calc_vertices(mc)) {
            fprintf(stderr,"mcubes_calc_vertices failed.\n");
            return -1;
    }


    mcubes_lewiner(mc);
    sumArea = mcubes_area(mc);
    sumScalarArea = mcubes_sum_scalar_area(mc);
    sumScalar2Area = mcubes_sum_scalar2_area(mc);

    meanScalar = sumScalarArea/sumArea;
    mean2Scalar= sumScalar2Area/sumArea;

    variance = mean2Scalar-meanScalar*meanScalar;
    

    printf("(sumScalarArea,sumArea,mean,sd)= %.12f %.12f %.12f %.12f\n",
            sumScalarArea,sumArea,meanScalar,sqrt(variance)
            );
    mcubes_destroy(mc);
    return 0;

}
