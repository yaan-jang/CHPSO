#include "chpso.h"
#include "chpso_search_space.h"
#include "chpso_objective_func.h"
#include "chpso_loop.cpp"

void about();
int main()
{
    about();              // Information of CHPSO in NiFoPS

    long start, end;
    start = clock();

    int seed;             // Random seed
    int fnum = 0;         // Number of objective functions
    int run, run_max;     // Counter and numbers of trials
    double gopt;          // Best-so-far
    double avg;           // Average of the best-so-fars 
    int dim;              // Dimension of an objective function
    double **ss = NULL;   // Search space
    run_max = 20;         // Numbers of trials
    for(int fidx = fnum; fidx<= fnum; fidx++) // Loop on problems
    {
        srand((unsigned)time(NULL));          // Initialize seed
        printf("Function No.: F%i\n", fidx+1);
            
        dim = 30;
        
        ss = create2darray(2, dim);
        chpso_search_space(ss, dim);
        avg = 0.0;
        seed = rand() % 100000;               // In the range 0 to 99999;
        for (run = 0; run < run_max; run++)
        {
            seed += run;
            chpso_loop(&gopt, ss, dim, seed);
            avg += gopt;
            if(run == 0)printf("The %2ist trial, best-so-far: %e\n", run+1, gopt);
            else if(run == 1)printf("The %2ind trial, best-so-far: %e\n", run+1, gopt);
            else if(run == 2)printf("The %2ird trial, best-so-far: %e\n", run+1, gopt);
            else printf("The %2ith trial, best-so-far: %e\n", run+1, gopt);
        }
        avg /= run_max;
        printf("\nThe average best-so-far of %i trials is: %e\n", run_max, avg);
        
    }
    end = clock();
    printf("Total time used: %f seconds.\n\n", ((double)(end - start)) / CLOCKS_PER_SEC);
    delete2darray(ss, 2, dim);
    return 0; 
}// End of main program

/*================== About =============*/
void about()
{
    // system("clear");
    printf( "======================================================================\n");
    printf( "|           Heterogeneous Particle Swarm Optimizer (CHPSO)           |\n");
    printf( "======================================================================\n");
    printf( "|           Author:        N. J. Cheung                              |\n");
    printf( "|           Email:         nyaam.ch@gmail.com                        |\n");
    printf( "|           Release:       Version 1.0                               |\n");
    printf( "|           Release Date:  Sep. 15, 2013.                            |\n");
    printf( "|           Last Modified: Aug. 15, 2014.                            |\n");
    printf( "|           http://godzilla.uchicago.edu/pages/ngaam/index.html      |\n");
    printf( "----------------------------------------------------------------------\n");
    printf( "   CHPSO Incentive Product - Sibe Executable Build                    \n");
    printf( "   Copyright (C) 2013-   by NJC (nyaam.ch@gmail.com)                  \n");
    printf( "   James Franck Institute, The University of Chicago                  \n");
    printf( "   929 East 57 Street, Chicago, IL 60637.                             \n");
    printf( "\n");
    printf( "   Permission is granted for all academic users and not-for-profit    \n");
    printf( "   institutions, who can copy, use, or modify this software for       \n");
    printf( "   any non-commercial purposes. Commercial users wishing an           \n");
    printf( "   evaluation copy should contact the OWNER. Commercial users MUST    \n");
    printf( "   license this software product after completing the license         \n");
    printf( "   agreement and sending it to nyaam.ch@gmail.com. Any other usage    \n");
    printf( "   is specifically prohibited and may constitute a violation of       \n");
    printf( "   United States and international copyright laws.                    \n");
    printf( "\n");
    printf( " References\n");
    printf( " [1] Ngaam J. Cheung, Xue-Ming Ding, Hong-Bin Shen.                   \n");
    printf( "     IEEE Transactions on Fuzzy Systems, 22(4): 919-933, 2014.        \n");
    printf( " [2] Ngaam J. Cheung, Zhen-Kai Xu, Xue-Ming Ding, Hong-Bin Shen.      \n");
    printf( "     European Journal of Operational Research, 247(2): 349-358, 2015. \n");
    printf( " [3] Ngaam J. Cheung, Xue-Ming Ding, Hong-Bin Shen.                   \n");
    printf( "     Journal of Computational Chemistry, 37(4): 426-436, 2016.        \n");
    printf( "----------------------------------------------------------------------\n");
    printf("\n");
}

