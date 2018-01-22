// Standard libraries
#include <stdlib.h>     // standard C library function
#include <stdio.h>
#include <cmath>        // math
#include <fstream>      // file i/o
#include <iostream>     // i/o stream
//#include <strstream>    // char *stream
#include <iomanip>      // i/o format
#include <string>       // Standard string library
#include <map>          // STL associative container
#include <limits.h>     // Noeric limit
#include <algorithm>    // STL algorithm
#include <vector>       // STL vector
#include <functional>   // STL functional
#include <time.h>
using namespace std;


// General settings
#define    NGen 2000     // Max number of generations

// Default coefficients. 
#define SNP 8            // Size of sub-population
#define acc 1.49445      // Acceleration factor
#define iwt_max  0.9     // Maximum inertia weight
#define iwt_min  0.35    // Minimum inertia weight
#define vel_limt 0.75    // Limit of velocity

// For random numbers
#define IM1 2147483563
#define IM2 2147483399
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

#define zero  1.0e-400  // To avoid numerical instabilities on some machines
#define infinity 1.e+100

// Create 2D array
double** create2darray(int rows, int cols)
{
    double** array = new double*[rows];
    for (int row=0; row < rows; row++)
        {
                array[row] = new double[cols];
        }
    return array;
}
// Delete 2D array
void delete2darray(double **ar, int rows, int cols)
{
    for (int row=0; row < rows; row++)
    {
        delete [] ar[row];
    }
    delete [] ar;
}

// Selection of better solution
void selection(double *a, double *b, int d)
{
  if (a[d] > b[d])
  {
    for(int k = 0; k < d; k++)
    {
      a[k] = b[k];
    }
    a[d] = b[d];
  }
}

// SRC-FUNCTION   :rnd_uni()                                        
// LONG_NAME      :random_uniform                                   
// AUTHOR         :(see below)                                      
//                                                                 
// DESCRIPTION    :rnd_uni() generates an equally distributed ran-  
//                dom number in the interval [0,1]. For further    
//                reference see Press, W.H. et al, Numerical     
//                Recipes in C, Cambridge University Press, 1992.                                                                                                          
// PARAMETERS     :*idum    serves as a seed value                                                                             
// PRECONDITIONS  :*idum must be negative on the first call.                                                                   
// POSTCONDITIONS :*idum will be changed 

float rnd_uni(long *idum)
{
    long j, k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) 
    {
        if (-(*idum) < 1) { *idum = 1; } 
        else { *idum = -(*idum); }
        idum2 = (*idum);
        for (j=NTAB+7;j >= 0;j--) 
        {
            k = (*idum)/IQ1;
            *idum = IA1*(*idum - k*IQ1) - k*IR1;
            if (*idum < 0) { *idum += IM1; }
            if (j < NTAB) { iv[j] = *idum; }
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum = IA1*(*idum - k*IQ1) - k*IR1;
    if (*idum < 0) { *idum += IM1; }
    k = idum2/IQ2;
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
    if (idum2 < 0) { idum2 += IM2; }
    j = iy/NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) { iy += IMM1; }
    double AM = (1.0/IM1);
    temp = (float) AM*iy;
    if (temp > RNMX) { return (float) RNMX; } 
    else { return temp; }
} 


