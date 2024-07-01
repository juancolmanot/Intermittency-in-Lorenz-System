#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "lorenz.h"

void lorenz(
    long double t,
    long double x[],
    long double f[],
    long double p[]
) {
    long double sigma, rho, beta;
    sigma = p[0];
    rho = p[1];
    beta = p[2];

    f[0] = sigma * (x[1] - x[0]);
    f[1] = ((rho - x[2]) * x[0]) - x[1];
    f[2] = (x[0] * x[1]) - beta * x[2];
}

void jaclorenz(
    long double t,
    long double x[],
    long double dfdx[][3],
    long double p[]
) {
    long double sigma = p[0];
    long double rho = p[1];
    long double beta = p[2];

    dfdx[0][0] = -sigma;
    dfdx[0][1] = sigma;
    dfdx[0][2] = 0.0;
    dfdx[1][0] = rho - x[2];
    dfdx[1][1] = -1.0;
    dfdx[1][2] = -x[0];
    dfdx[2][0] = x[1];
    dfdx[2][1] = x[0];
    dfdx[2][2] = -beta;   
}

int lorenz_gsl(
    double t,
    const double x[],
    double f[],
    void *par
) {
    double *p = (double *)par;
    double sigma, rho, beta;
    sigma = p[0];
    rho = p[1];
    beta = p[2];

    f[0] = sigma * (x[1] - x[0]);
    f[1] = ((rho - x[2]) * x[0]) - x[1];
    f[2] = (x[0] * x[1]) - beta * x[2];

    return GSL_SUCCESS;
}

int jaclorenz_gsl(
    double t,
    const double x[],
    double *dfdx,
    double dfdt[],
    void *par
) {
    double *p = (double *)par;
    double sigma = p[0];
    double rho = p[1];
    double beta = p[2];

    dfdx[0 * 3 + 0] = -sigma;
    dfdx[0 * 3 + 1] = sigma;
    dfdx[0 * 3 + 2] = 0.0;
    dfdx[1 * 3 + 0] = rho - x[2];
    dfdx[1 * 3 + 1] = -1.0;
    dfdx[1 * 3 + 2] = -x[0];
    dfdx[2 * 3 + 0] = x[1];
    dfdx[2 * 3 + 1] = x[0];
    dfdx[2 * 3 + 2] = -beta;

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;

    return GSL_SUCCESS;
}