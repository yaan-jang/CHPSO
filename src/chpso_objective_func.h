void test(double *eval, double xx[], int dim, int function) 
{    // Evaluate the fitness of the particle
    *eval = 0.0; 
    if (function == 0) // Sphere
    {
        for (int j = 0; j < dim; j++) 
        {
            *eval = *eval + xx[j] * xx[j]; 
        }
    }
}
