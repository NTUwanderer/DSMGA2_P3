/***************************************************************************
 *   Copyright (C) 2015 by Tian-Li Yu                                      *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/


#include <math.h>
#include <iostream>
#include <cstdlib>

#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#include "chromosome.h"

using namespace std;


int
main (int argc, char *argv[]) {


    if (argc < 8) {
        printf ("DSMGA2 ell function maxGen maxFe repeat display rand_seed [pNum]\n");
        printf ("function: \n");
        printf ("     ONEMAX:  0\n");
        printf ("     MK    :  1\n");
        printf ("     FTRAP :  2\n");
        printf ("     CYC   :  3\n");
        printf ("     NK    :  4\n");
        printf ("     SPIN  :  5\n");
        printf ("     SAT   :  6\n");

        return -1;
    }

    int ell = atoi (argv[1]); // problem size
    //int nInitial = 0; // initial population size
    int nInitial = (int)(10*(log(ell)/log(2.71828))+2); // initial population size
    //int nInitial = (int) 1; // initial population size
    //int nInitial = (int) 2; // initial population size
    int fffff = atoi (argv[2]); // function
    int maxGen = atoi (argv[3]); // max generation
    int maxFe = atoi (argv[4]); // max fe
    int repeat = atoi (argv[5]); // how many time to repeat
    int display = atoi (argv[6]); // display each generation or not
    int rand_seed = atoi (argv[7]);  // rand seed
    int pNum = (argc==9)? atoi (argv[8]):0;


    if (fffff == 4) {

        char filename[200];
        sprintf(filename, "../../NK_Instance/pnk%d_%d_%d_%d", ell, 4, 3, pNum);

        printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");
        loadNKWAProblem(fp, &nkwa);
        fclose(fp);
    }
    else if (fffff == 5) {

        int side = (int)sqrt(ell+0.001);
        char filename[200];
        FILE *fp;
        if (side <= 17) {
            sprintf(filename, "../../SPIN/%dx%d/%d_%d.%d", side, side, side, ell, pNum);
            printf("Loading: %s\n", filename);
            fp = fopen(filename, "r");
            loadSimonsSpinGlassInstance(fp, &mySpinGlassParams);
            fclose(fp);
        }
        else {
            sprintf(filename, "../../SPIN/%dx%d/s_%dx%d.%d", side, side, side, side, pNum);
            printf("Loading: %s\n", filename);
            fp = fopen(filename, "r");
            loadSpinGlassInstance(fp, &mySpinGlassParams);
            fclose(fp);
        }
    }
    else if (fffff == 6) {

        char filename[200];
        sprintf(filename, "../../SAT/uf%d/uf%d-0%d.cnf", ell, ell, pNum);

        printf("Loading: %s\n", filename);
        loadSAT(filename, &mySAT);
    }

    if (rand_seed != -1)  // time
        myRand.seed((unsigned long)rand_seed);

    int i;

    Statistics stGen, stFE, stLSFE, stN;
    int usedGen;

    int failNum = 0;



    for (i = 0; i < repeat; i++) {

        DSMGA2 ga (ell, nInitial, maxGen, maxFe, fffff);

        if (display == 1)
            usedGen = ga.doIt (true);
        else
            usedGen = ga.doIt (false);


        if (!ga.foundOptima()) {
            failNum++;
            printf ("-");
        } else {
            stN.record (ga.nCurrent);
            stFE.record (Chromosome::hitnfe);
            stLSFE.record (Chromosome::lsnfe);
            stGen.record (usedGen);
            printf ("+");
        }

        fflush (NULL);

    }

    if (fffff == 4)
        freeNKWAProblem(&nkwa);
    else if (fffff == 5)
        freeSpinGlassInstance(&mySpinGlassParams);


    cout << endl;
    printf ("\n");
    printf ("Gen: %f\n", stGen.getMean ());
    printf ("FailNum: %d\n", failNum);
    printf ("LSNFE: %f\n", stLSFE.getMean());
    printf ("NFE: %f\n", stFE.getMean());
    printf ("Popu: %f\n", stN.getMean());


    return EXIT_SUCCESS;
}
