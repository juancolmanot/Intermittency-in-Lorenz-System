#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../../../../../cursos/modulosc/lorenz.h"
// #include "../../../../../cursos/modulosc/gsl_utilities.h"

int main(void) {

    FILE *poincare = fopen("../datafiles/poincare_section_lorenz.dat", "w");
    FILE *evol = fopen("../datafiles/lorenz_evol.dat", "w");

    unsigned long int seed = (unsigned long int)time(NULL);

    Parameters *p = malloc(sizeof(Parameters));

    p->a = 10.0;
    p->b = 166.2;
    p->c = 8.0 / 3.0;

    double rel_err = 1e-14, abs_err = 1e-12;
    double h = 1e-6;
    const gsl_odeiv2_step_type *integrator = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(integrator, 3);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(abs_err, rel_err);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(3);

    gsl_odeiv2_system sys = {lorenz, jaclorenz, 3, p};

    double t = 0.0, t1 = 10000000.0;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, seed);
    double x[3];

    x[0] = 2.0;
    x[1] = -1.0;
    x[2] = 150.0;

    /* Passes counter */
    unsigned int counter = 0;
    unsigned int pass_target = 100000;

    /* Arrays for x, y, z */
    double xi[3], yi[3], zi[3];

    /* Arrays for adjustment */
    double xfit[3], yfit[3], zfit[3];
    double *coefficients = calloc(3, sizeof(double));
    double a_c, b_c, c_c;
    double *yreg, *zreg;
    yreg = calloc(pass_target, sizeof(double));
    zreg = calloc(pass_target, sizeof(double));
    gsl_matrix *ty = gsl_matrix_alloc(2, pass_target);
    gsl_matrix *tz = gsl_matrix_alloc(2, pass_target);

    /* Poincaré section */
    double xp = 0;
    xi[0] = x[0];
    yi[0] = x[1];
    zi[0] = x[2];

    /* Initialize the first three points */
    for (unsigned int i = 1; i < 3; i++) {
        int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h, x);
        xi[i] = x[0];
        yi[i] = x[1];
        zi[i] = x[2];
        if (status != GSL_SUCCESS) {
            printf("error initializing, return value=%d\n", status);
            break;
        }
    }

    while (counter < pass_target) {
        int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h, x);

        /* Update values */
        xi[0] = xi[1];
        xi[1] = xi[2];
        xi[2] = x[0];
        
        yi[0] = yi[1];
        yi[1] = yi[2];
        yi[2] = x[1];

        zi[0] = zi[1];
        zi[1] = zi[2];
        zi[2] = x[2];

        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n", status);
            break;
        }

        if (t > 1000.0) {
            if (xi[1] > xp) {
                if (xi[0] < xp) {
                    xfit[0] = xi[0];
                    xfit[1] = xi[1];
                    xfit[2] = xi[2];
                    yfit[0] = yi[0];
                    yfit[1] = yi[1];
                    yfit[2] = yi[2];
                    zfit[0] = zi[0];
                    zfit[1] = zi[1];
                    zfit[2] = zi[2];
                    gsl_regression_quadratic(xfit, yfit, 3, coefficients);
                    a_c = coefficients[2];
                    b_c = coefficients[1];
                    c_c = coefficients[0];
                    // yreg[counter] = yi[1];
                    yreg[counter] = a_c * xp * xp + b_c * xp + c_c;
                    gsl_regression_quadratic(xfit, zfit, 3, coefficients);
                    a_c = coefficients[2];
                    b_c = coefficients[1];
                    c_c = coefficients[0];
                    zreg[counter] = a_c * xp * xp + b_c * xp + c_c;
                    fprintf(evol, "%d  %5.5f  %5.5f\n", counter, yreg[counter], zreg[counter]);
                    // printf("%d  %5.5f  %5.5f\n", counter, yreg[counter], zreg[counter]);
                    // zreg[counter] = zi[1];
                    gsl_matrix_set(ty, 0, counter, t);
                    gsl_matrix_set(ty, 1, counter, yreg[counter]);
                    gsl_matrix_set(tz, 0, counter, t);
                    gsl_matrix_set(tz, 1, counter, zreg[counter]);
                    counter++;
                    if (counter % 5000 == 0) {
                        printf("counter: %d\n", counter);
                    }
                }
            }
        }
    }

    /* Laminar region */
    FILE *reinjection = fopen("../datafiles/reinject_lorenz.dat", "w");
    double yf = 41.2861;
    double clam = 1.85;
    unsigned int count_reinject = 0;

    for (unsigned int i = 1; i < pass_target; i++) {
        if (yreg[i] >= yf - clam && yreg[i] <= yf + clam) {
            if (yreg[i - 1] > yf + clam || yreg[i - 1] < yf - clam) {
                count_reinject++;
                fprintf(reinjection, "%d %5.15f  %5.15f\n", count_reinject, yreg[i - 1], yreg[i]);
            }
        }
    }

    printf("Count reinject: %d\n", count_reinject);

    FILE *RPD = fopen("../datafiles/rpd_lorenz.dat", "w");
    unsigned int n_bins = 750;
    double ymin, ymax, yj, dy;
    ymin = yf - clam;
    ymax = yf + clam;
    dy = (ymax - ymin) / (double)(n_bins - 1);
    yj = ymin;
    double *yreinjected = calloc(n_bins, sizeof(double));
    double knorm = 0;
    
    /* Increment y from min to max to generate bins */
    for (unsigned int i = 0; i < n_bins; i++) {

        /* Loop through whole array of y fitted */
        for (unsigned int j = 1; j < pass_target; j++) {
            
            /* Ask for belonging to bin and reinjection */
            if (yreg[j] > yj && yreg[j] <= yj + dy) {
                if (yreg[j] >= yf - clam && yreg[j] <= yf + clam) {
                    if (yreg[j - 1] < yf - clam || yreg[j - 1] > yf + clam) {
                        yreinjected[i]++;
                    }
                }
            }
        }
        yreinjected[i] = yreinjected[i] / (double)count_reinject;
        fprintf(RPD, "%5.15f  %5.15f\n", yj, yreinjected[i]);
        yj += dy;
    }

    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    gsl_rng_free(rng);
    fclose(poincare);
    fclose(reinjection);
    fclose(RPD);
    fclose(evol);
    return 0;
}