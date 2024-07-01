#ifndef LORENZ_H
#define LORENZ_H

typedef struct {
    double sigma;
    double rho;
    double beta;
} lorenz_par;

void lorenz(
    long double t,
    long double x[],
    long double f[],
    long double p[]
);

void jaclorenz(
    long double t,
    long double x[],
    long double dfdx[][3],
    long double p[]
);

int lorenz_gsl(
    double t,
    const double x[],
    double f[],
    void *par
);

int jaclorenz_gsl(
    double t,
    const double x[],
    double *dfdx,
    double dfdt[],
    void *par
);

#endif