#ifndef PHYSPROP_H
#define PHYSPROP_H

#define NMIN 39
#define NMODE 19
#define PI 3.14159265358979
#define ISQRTPI 1.12837916709551 
#define gravacc 10

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Function Prototypes */
double ref_density(int i, double Tref);
void adj_mf(double T, double mf[], double *cmf);
double adj_alpha(double t, double T);
void computedensity(int aflag, double t, double Tref, double P, double T,
  double mf[], double *rho, double *alpha, double *vf);
void computeheatcap(int aflag, double t, double Tref,
  double P, double T, double vf[], double *Cp);
void computethermcond(double P, double T, double vf[], double *lambda, int opt);
void densityPT(int i, double Tref, double P, double T, double *rhoPT,
  double *alphaPT, double R);
double densityP(int i, double P, double rhor);
double pressure(int i, double rho0, double rhoP);
double dpressure(int i, double rho0, double rhoP);
double thermexp(int i, double T);
double ithermexp(int i, double T, double Tref);
double dthermexp(int i, double T);
double heatcap(int i, double Tref, double P, double T, double R);
void phasefrac(double vf[], double *X);
void bvalue(double vf[], double *b);
void condconst(double vf[], double T, double *lambda0);
void thermcond(double P, double T, double vf[],
  double *lamda_lat, double *lambda_rad);
//void ithermcond(double P, double T, double vf[], double *G);
void dthermcond(double P, double T, double vf[], double *dlambda);
double adiabat(double Tp, double alphaPT, double Cp);
double sp2gt(double T);
void modevalue(double vmin[], double vf[], double *vmode);
void modefrac(double vf[], double *modefrac);
int checkwts(double vf[]);
double wtmean(double w[], double vmode[]);
double geomean(double w[], double vmode[]);
void vf2mf(double vf[], double mf[], double rho[]);
void mf2vf(double mf[], double vf[], double rho[]);

#endif
