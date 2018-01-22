void chpso_loop(double *opt, double **ss, int dim, int seed) 
{
    long rnd_uni_init = -(long)seed;  // Initialization of rnd_uni()

    double iwt = 0.0;                       // Inertia weight
    double acc1 = acc,                // Acceleration factors
           acc2 = acc;

    bool stop = true;                 // Flag for stop loop
    int cnt = 0;                      // Iteration number (time step)
    double gamma, gamma_1, gamma_2;   // Fitness trackers
  
    int rfc = 0, RFC = 5;             // Negative operator

    double rmin, rmax;                // Lower and upper boundaries
    double mv[dim];                   // Limits of velocity
    // for(int j = 0; j < dim; j++) mv[j] = vel_limt * (rmax - rmin);
    
    // Definitions of the four populations
    double **pos_a  = NULL, **velocity_a = NULL, **localpos_a = NULL, **global_a = NULL;
    double **pos_b  = NULL, **velocity_b = NULL, **localpos_b = NULL, **global_b = NULL;
    double **pos_c  = NULL, **velocity_c = NULL, **localpos_c = NULL, **global_c = NULL;
    double **pos_d  = NULL, **velocity_d = NULL, **localpos_d = NULL, **global_d = NULL;
    double **global = NULL;

    pos_a = create2darray(SNP, dim+1);
    pos_b = create2darray(SNP, dim+1);
    pos_c = create2darray(SNP, dim+1);
    pos_d = create2darray(SNP, dim+1);

    velocity_a = create2darray(SNP, dim);
    velocity_b = create2darray(SNP, dim);
    velocity_c = create2darray(SNP, dim);
    velocity_d = create2darray(SNP, dim);

    localpos_a = create2darray(SNP, dim+1);
    localpos_b = create2darray(SNP, dim+1);
    localpos_c = create2darray(SNP, dim+1);
    localpos_d = create2darray(SNP, dim+1);

    global_a = create2darray(1, dim+1);
    global_b = create2darray(1, dim+1);
    global_c = create2darray(1, dim+1);
    global_d = create2darray(1, dim+1);

    global = create2darray(1, dim+1);

    // Initialization of global information
    for(int j = 0; j < dim; j++)
    {
        rmin = ss[0][j];
        rmax = ss[1][j];
        global_a[0][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
        global_b[0][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
        global_c[0][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
        global_d[0][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
        global[0][j]   = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
    }
    // global_a[0][dim] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
    // global_b[0][dim] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
    // global_c[0][dim] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
    // global_d[0][dim] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
    // global[0][dim]   = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);

    test(&global_a[0][dim], global_a[0], dim, 0); 
    test(&global_b[0][dim], global_b[0], dim, 0);
    test(&global_c[0][dim], global_c[0], dim, 0);
    test(&global_d[0][dim], global_d[0], dim, 0);
    test(&global[0][dim],   global[0],   dim, 0);

    // Initializations of positions and velocities    
    for (int i = 0; i < SNP; i++)   
    {
        for(int j = 0; j < dim; j++)
        {
            rmin = ss[0][j];
            rmax = ss[1][j];
            mv[j] = vel_limt * (rmax - rmin);
            
            pos_a[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            pos_b[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            pos_c[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            pos_d[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);

            velocity_a[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            velocity_b[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            velocity_c[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            velocity_d[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);

            localpos_a[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            localpos_b[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            localpos_c[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
            localpos_d[i][j] = rmin + (rmax - rmin) * rnd_uni(&rnd_uni_init);
        } 

        // Evaluations of populations
        test(&pos_a[i][dim], pos_a[i], dim, 0);
        test(&pos_b[i][dim], pos_b[i], dim, 0);
        test(&pos_c[i][dim], pos_c[i], dim, 0);
        test(&pos_d[i][dim], pos_d[i], dim, 0);
        // Evaluations of local best records
        test(&localpos_a[i][dim], localpos_a[i], dim, 0); 
        test(&localpos_b[i][dim], localpos_b[i], dim, 0); 
        test(&localpos_c[i][dim], localpos_c[i], dim, 0); 
        test(&localpos_d[i][dim], localpos_d[i], dim, 0); 
        // Local best selections
        selection(localpos_a[i], pos_a[i], dim);
        selection(localpos_b[i], pos_b[i], dim);
        selection(localpos_c[i], pos_c[i], dim);
        selection(localpos_d[i], pos_d[i], dim);
        // Global best selections
        selection(global_a[0], pos_a[i], dim);
        selection(global_b[0], pos_b[i], dim);
        selection(global_c[0], pos_c[i], dim);
        selection(global_d[0], pos_d[i], dim);
    }
    selection(global[0], global_a[0], dim);
    selection(global[0], global_b[0], dim);
    selection(global[0], global_c[0], dim);
    selection(global[0], global_d[0], dim);

    // Main loop
    while (stop == true) 
    {    
        cnt = cnt + 1;
        iwt = iwt_max - (iwt_max - iwt_min) * cnt/NGen;
        for(int i = 0; i < SNP; i++)
        {
            for(int j = 0; j < dim; j++)
            {
                velocity_a[i][j] =  iwt  * velocity_a[i][j]
                                   + acc1 * (localpos_a[i][j] - pos_a[i][j]) * rnd_uni(&rnd_uni_init)
                                   + acc2 * (global[0][j]     - pos_a[i][j]) * rnd_uni(&rnd_uni_init);
                velocity_a[i][j] =  (velocity_a[i][j] >  mv[j]) * mv[j]
                                  + (velocity_a[i][j] <= mv[j]) * velocity_a[i][j];
                velocity_a[i][j] =  (velocity_a[i][j] < -mv[j]) * (-mv[j])
                                  + (velocity_a[i][j] >=-mv[j]) * velocity_a[i][j];
                pos_a[i][j] = pos_a[i][j] + velocity_a[i][j];
                // pos_a[i][j] =  (pos_a[i][j] >  rmax) * rmax
                //              + (pos_a[i][j] <= rmax) * pos_a[i][j];
                // pos_a[i][j] =  (pos_a[i][j] <  rmin) * rmin
                //              + (pos_a[i][j] >= rmin) * pos_a[i][j];

                velocity_b[i][j] =  iwt  * velocity_b[i][j]
                                  + acc1 * (localpos_b[i][j] - pos_b[i][j]) * rnd_uni(&rnd_uni_init)
                                  + acc2 * (global[0][j]     - pos_b[i][j]) * rnd_uni(&rnd_uni_init);
                velocity_b[i][j] =  (velocity_b[i][j] >  mv[j]) * mv[j]
                                  + (velocity_b[i][j] <= mv[j]) * velocity_b[i][j];
                velocity_b[i][j] =  (velocity_b[i][j] < -mv[j]) * (-mv[j])
                                  + (velocity_b[i][j] >=-mv[j]) * velocity_b[i][j];
                pos_b[i][j] = pos_b[i][j] + velocity_b[i][j];
                // pos_b[i][j] =  (pos_b[i][j] >  rmax) * rmax
                //              + (pos_b[i][j] <= rmax) * pos_b[i][j];
                // pos_b[i][j] =  (pos_b[i][j] <  rmin) * rmin
                //              + (pos_b[i][j] >= rmin) * pos_b[i][j];
            }
            test(&pos_a[i][dim], pos_a[i], dim, 0);
            test(&pos_b[i][dim], pos_b[i], dim, 0);
            gamma_1 = pos_a[i][dim];
            gamma_2 = pos_b[i][dim];
            gamma   = gamma_1 + gamma_2;

            for(int j = 0; j < dim; j++)
            {
                velocity_c[i][j] =  iwt  * (( gamma * velocity_a[i][j] / gamma_1 
                                            + gamma * velocity_b[i][j] / gamma_2) + velocity_c[i][j])
                                   + acc1 * (localpos_c[i][j] - pos_c[i][j]) * rnd_uni(&rnd_uni_init)
                                  + acc2 * (global[0][j]     - pos_c[i][j]) * rnd_uni(&rnd_uni_init);
                velocity_c[i][j] =  (velocity_c[i][j] >  mv[j]) * mv[j]
                                  + (velocity_c[i][j] <= mv[j]) * velocity_c[i][j];
                velocity_c[i][j] =  (velocity_c[i][j] < -mv[j]) * (-mv[j])
                                  + (velocity_c[i][j] >=-mv[j]) * velocity_c[i][j];
                pos_c[i][j] = pos_c[i][j] + velocity_c[i][j];
                // pos_c[i][j] =  (pos_c[i][j] >  rmax) * rmax
                //              + (pos_c[i][j] <= rmax) * pos_c[i][j];
                // pos_c[i][j] =  (pos_c[i][j] <  rmin) * rmin
                //              + (pos_c[i][j] >= rmin) * pos_c[i][j];
            // }
            // test(&pos_c[i][dim], pos_c[i], dim, 0);

            // for(int j = 0; j < dim; j++)
            // {
                velocity_d[i][j] =  velocity_a[i][j] + velocity_b[i][j] - velocity_c[i][j];
                velocity_d[i][j] =  (velocity_d[i][j] >  mv[j]) * mv[j]
                                  + (velocity_d[i][j] <= mv[j]) * velocity_d[i][j];
                velocity_d[i][j] =  (velocity_d[i][j] < -mv[j]) * (-mv[j])
                                  + (velocity_d[i][j] >=-mv[j]) * velocity_d[i][j];
                pos_d[i][j] = pos_d[i][j]/6.0 + localpos_d[i][j]/3.0 + global[0][j]/2.0 + velocity_d[i][j];
                // pos_d[i][j] =  (pos_d[i][j] >  rmax) * rmax
                //              + (pos_d[i][j] <= rmax) * pos_d[i][j];
                // pos_d[i][j] =  (pos_d[i][j] <  rmin) * rmin
                //              + (pos_d[i][j] >= rmin) * pos_d[i][j];
            }
            test(&pos_c[i][dim], pos_c[i], dim, 0);
            test(&pos_d[i][dim], pos_d[i], dim, 0);
            // Local best selections
            selection(localpos_a[i], pos_a[i], dim);
            selection(localpos_b[i], pos_b[i], dim);
            selection(localpos_c[i], pos_c[i], dim);
            selection(localpos_d[i], pos_d[i], dim);
            // Global best selections
            selection(global_a[0], pos_a[i], dim);
            selection(global_b[0], pos_b[i], dim);
            selection(global_c[0], pos_c[i], dim);
            selection(global_d[0], pos_d[i], dim);
        }
        selection(global[0], global_a[0], dim);
        selection(global[0], global_b[0], dim);
        selection(global[0], global_c[0], dim);
        selection(global[0], global_d[0], dim);

        rfc = rfc + 1;
        if(rfc >= RFC)
        {
            int widx_a = 0, widx_b = 0, widx_c = 0, widx_d = 0;
            double wval_a = pos_a[widx_a][dim];
            double wval_b = pos_b[widx_b][dim];
            double wval_c = pos_c[widx_c][dim];
            double wval_d = pos_d[widx_d][dim];
            rfc = 0;
            for(int i = 1; i < SNP; i++)
            {
                if(pos_a[i][dim] > wval_a)
                {
                    widx_a = i;
                    wval_a = pos_a[i][dim];
                }
                if(pos_b[i][dim] > wval_b)
                {
                    widx_b = i;
                    wval_b = pos_b[i][dim];
                }
                if(pos_c[i][dim] > wval_c)
                {
                    widx_c = i;
                    wval_c = pos_c[i][dim];
                }
                if(pos_d[i][dim] > wval_d)
                {
                    widx_d = i;
                    wval_d = pos_d[i][dim];
                }
            }
            for(int j = 0; j< dim; j++) 
            {
                pos_a[widx_a][j] = -global[0][j];
                pos_b[widx_b][j] = -global[0][j];
                pos_c[widx_c][j] = -global[0][j];
                pos_d[widx_d][j] = -global[0][j];
            }
        }
        if (cnt > NGen) stop = false;
    }
    *opt = global[0][dim];

    delete2darray(pos_a, SNP, dim+1);
    delete2darray(pos_b, SNP, dim+1);
    delete2darray(pos_c, SNP, dim+1);
    delete2darray(pos_d, SNP, dim+1);

    delete2darray(velocity_a, SNP, dim+1);
    delete2darray(velocity_b, SNP, dim+1);
    delete2darray(velocity_c, SNP, dim+1);
    delete2darray(velocity_d, SNP, dim+1);

    delete2darray(localpos_a, SNP, dim+1);
    delete2darray(localpos_b, SNP, dim+1);
    delete2darray(localpos_c, SNP, dim+1);
    delete2darray(localpos_d, SNP, dim+1);

    delete2darray(global_a, 1, dim+1);
    delete2darray(global_b, 1, dim+1);
    delete2darray(global_c, 1, dim+1);
    delete2darray(global_d, 1, dim+1);

    delete2darray(global, 1, dim+1);   
 
}
