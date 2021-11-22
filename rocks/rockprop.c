#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "physprop.h"

void usage(void);
char *getrecord(FILE *fptr);
int loadrock(FILE *fin, double *mf, int *type, int *kopt);
void outputmineralogy(int N, double mf[]);
void readstep(FILE *fin, int *nT, double *T0, double *dT,
    int *nP, double *P0, double *dP);
void numrocks(FILE *fin, int *N, int *kopt);
void stepcompute(FILE *fin, double mf[], int kopt);
void PTcompute(FILE *fin, double mf[], int kopt);
void outputheader(int opt);
void outputresults(double P, double T, double alpha, double rho, double Cp, double lambda, double kappa);

/* Mineral index key                                                                            *
 *  0: alpha-quartz      10: Fe-epidote          20: wollastonite        30: grossular          *
 *  1: beta-quartz       11: clinochlore         21: enstatite           31: spessartine        *
 *  2: orthoclase        12: daphnite            22: ferrosillite        32: aragonite          *
 *  3: sanadine          13: hornblende          23: forsterite          33: calcite            *
 *  4: albite            14: tremolite           24: fayalite            34: dolomite           *
 *  5: anorthite         15: ferro-actinolite    25: spinel              35: andalucite         *
 *  6: muscovite         16: pargasite           26: magnetite           36: kyanite            *
 *  7: paragonite        17: diopside            27: hercynite           37: sillimanite        *
 *  8: phlogopite        18: hedenbergite        28: pyrope              38: rutile             *
 *  9: annite            19: jadite              29: almandine                                  */
char minname[NMIN][30] = {
    "(qz) alpha-quartz",      "(qz) beta-quartz",     "(or) orthoclase",
    "(or) sanidine",          "(plag) albite",        "(plag) anorthite",
    "(musc) muscovite",       "(musc) paragonite",    "(bt) phlogopite",
    "(bt) annite",            "(ep) epidote",         "(chl) clinochore",
    "(chl) daphnite",         "(amph) hornblende",    "(amph) tremolite",
    "(amph) ferroactinolite", "(amph) pargasite",     "(cpx) diopside",
    "(cpx) hedenbergite",     "(cpx) jadite",         "(wo) wollastonite",
    "(opx) enstatite",        "(opx) ferrosilite",    "(ol) fosterite",
    "(ol) fayalite",          "(sp) spinel",          "(sp) magnetite",
    "(sp) hercynite",         "(gt) pyrope",          "(gt) almandine",
    "(gt) grossular",         "(gt) spessartine",     "(carb) aragonite",
    "(carb) calcite",         "(carb) dolomite",      "(alsi) andalucite",
    "(alsi) kyanite",         "(alsi) sillimanite",   "(ru) rutile"};


void usage(void) {
    fprintf(stderr, "\nusage: rockprop infile [outfile]\n\n");
    fprintf(stderr, " Mineral index key:\n\n");
    fprintf(stderr, "    0: alpha-quartz        20: wollastonite\n");
    fprintf(stderr, "    1: beta-quartz         21: enstatite\n");
    fprintf(stderr, "    2: orthoclase          22: ferrosillite\n");
    fprintf(stderr, "    3: sanadine            23: forsterite\n");
    fprintf(stderr, "    4: albite              24: fayalite\n");
    fprintf(stderr, "    5: anorthite           25: spinel\n");
    fprintf(stderr, "    6: muscovite           26: magnetite\n");
    fprintf(stderr, "    7: paragonite          27: hercynite\n");
    fprintf(stderr, "    8: phlogopite          28: pyrope\n");
    fprintf(stderr, "    9: annite              29: almandine\n");
    fprintf(stderr, "   10: Fe-epidote          30: grossular\n");
    fprintf(stderr, "   11: clinochlore         31: spessartine\n");
    fprintf(stderr, "   12: daphnite            32: aragonite\n");
    fprintf(stderr, "   13: hornblende          33: calcite\n");
    fprintf(stderr, "   14: tremolite           34: dolomite\n");
    fprintf(stderr, "   15: ferro-actinolite    35: andalucite\n");
    fprintf(stderr, "   16: pargasite           36: kyanite\n");
    fprintf(stderr, "   17: diopside            37: sillimanite\n");
    fprintf(stderr, "   18: hedenbergite        38: rutile\n");
    fprintf(stderr, "   19: jadite\n\n");
}


/* getrecord.c                                                  *
 *                                                              *
 * Get a single record from a file.                             *
 * Returns NULL on EOF, or for comment/blank lines.             *
 *                                                              *
 * Copyright (c) 2000, 2002 by                                  *
 * Paul Gettings,                                               *
 * Department of Geology & Geophysics,                          *
 * University of Utah.                                          *
 *                                                              *
 * All Rights Reserved.                                         *
 *                                                              *
 * This is free software; you can redistribute it and/or        *
 * modify it under the terms of the GNU General Public          *
 * License, version 2, as published by the Free Software        *
 * Foundation.                                                  */
char *getrecord(FILE *fptr) {
    static char line[1025];
    char str[1025];
    int i;

    fgets(line, 1025, fptr);

    /* EOF */
    if(feof(fptr)) { return(NULL); }

    /* strip leading whitespace */
    i = sscanf(line, "%1024s", str); 
    //printf("i: %i str: '%s' | line: '%s'\n",i,str,line);

    /* empty line/EOF */
    if (i <= 0) { return(NULL); }
    /* comments */
    if ((str[0] == '#') || (str[0] == '%') || (strcmp(str,"\n") == 0)) { return(NULL); }

    return(line);
}


/* LOADROCK */
int loadrock(FILE *fin, double *mf, int *type, int *kopt) {

    int i,j;
    int N;
    char *line;
    double temp;
    double s = 0;

    /* Number of Minerals */
    while ((line = getrecord(fin)) == NULL) { continue; } 
    sscanf(line, "%i", &N);         

    /* Initialize Mineral Fractions */
    for (j = 0; j < NMIN; j++) {
        mf[j] = 0;
    }

    /* Read Mineral ID and Mass Fraction */
    for (j = 0; j < N; j++) {
        while ((line = getrecord(fin)) == NULL) { continue; } 

        sscanf(line, "%i %lg", &i, &temp);

        mf[i] = temp/100;
        s = s + mf[i];
    }
    /* Normalize mineral fractions if necessary */
    if (s != 1) {
        fprintf(stdout,"%% Warning: Mineral fractions have been normalized.\n");
        fprintf(stdout,"%%   Normalization factor: %4.2lf.\n",s);
        for (j = 0; j < NMIN; j++) {
            mf[j] = mf[j]/s;
        }
    }

    /* Read PT - type (0: discrete steps/ 1: exact values) and load     *
     * conductivity option (0: effective / 1: lattice / 2: radiative)   */
    while ((line = getrecord(fin)) == NULL) { continue; } 
    if (*kopt == -1) {
        sscanf(line, "%i %i", type, kopt);
    } else {
        sscanf(line, "%i", type);
    }

    return(0);
}


void readstep(FILE *fin, int *nT, double *T0, double *dT,
    int *nP, double *P0, double *dP) {

    char *line;

    /* Temperature (in kelvins) */
    while ((line = getrecord(fin)) == NULL) { continue; } 
    sscanf(line, "%i %lf %lf", nT,dT,T0);

    /* Pressure (in GPa) */
    while ((line = getrecord(fin)) == NULL) { continue; } 
    sscanf(line, "%i %lf %lf", nP,dP,P0);
}



void numrocks(FILE *fin, int *N, int *kopt) {

    char *line;

    while ((line = getrecord(fin)) == NULL) { continue; } 

    sscanf(line,"%i %i",N,kopt);
}


void outputmineralogy(int N, double mf[]) {

    int i;

    fprintf(stdout,"%% Rock: %i\n",N+1);

    for (i = 0; i < NMIN; i++) {
        if (mf[i] != 0) {
            fprintf(stdout,"%%  %-30s : %7.3lf\n",minname[i],mf[i]*100);
        }
    }
}

void outputheader(int opt) {
    switch (opt) {
        case 0 :
            fprintf(stdout,"%%P       T   alpha      rho     Cp      k_eff  kappa\n");
            break;
        case 1 :
            fprintf(stdout,"%%P        T   alpha      rho     Cp      k_lat  kappa\n");
            break;
        case 2 :
            fprintf(stdout,"%%P        T   alpha      rho     Cp      k_rad  kappa\n");
            break;
        default :
            fprintf(stderr,"Warning (outputheader): unknown option.\n");
            break;
    }
}

void outputresults(double P, double T, double alpha, double rho, double Cp, double lambda, double kappa) {
    fprintf(stdout,"%7.4lf %4.0lf %6.4le %7.2lf %7.2lf %6.3lf %6.4le \n",P,T,alpha,rho,Cp,lambda,kappa);
}

void stepcompute(FILE *fin, double mf[], int kopt) {

    int j,k;

    int nP,nT;
    double dP,P0,P;
    double dT,T0,T;

    double alpha,rho,Cp,lambda,kappa;
    static double dummy[NMIN];

    readstep(fin, &nT,&T0,&dT, &nP,&P0,&dP);

    for (j = 0; j < nP; j++) {
        P = P0 + j*dP;
        for (k = 0; k < nT; k++) {
            T = T0 + k*dT;
            if (T > 845 && mf[0] != 0) {
                mf[1] = mf[0];
                mf[0] = 0;
            } else if (T <= 845 && mf[1] != 0) {
                mf[0] = mf[1];
                mf[1] = 0;
            }

            computedensity(0,0, 1573, P,T, mf, &rho,&alpha, dummy);
            computeheatcap(0,0, 1573, P,T, mf, &Cp);
            computethermcond(P,T, mf, &lambda,kopt);
            kappa = lambda / (rho * Cp);

            outputresults(P,T, alpha,rho,Cp,lambda,kappa);
        }
    }
}


void PTcompute(FILE *fin, double mf[], int kopt) {

    char *line;
    int i,N;

    double P,T;
    double alpha,rho,Cp,lambda,kappa;
    static double dummy[NMIN];

    while ((line = getrecord(fin)) == NULL) { continue; } 
    sscanf(line, "%i", &N);


    for (i = 0; i < N; i++) {
        while ((line = getrecord(fin)) == NULL) { continue; } 
        sscanf(line, "%lg %lg", &P, &T);

        computedensity(0,0, 1573, P,T, mf, &rho,&alpha, dummy);
        computeheatcap(0,0, 1573, P,T, mf, &Cp);
        computethermcond(P,T, mf, &lambda, kopt);
        kappa = lambda / (rho * Cp);

        outputresults(P,T, alpha,rho,Cp,lambda,kappa);
    }
}


int main(int argc, char *argv[]) {

    char infile[512];                   // input file of parameters
    FILE *fin;

    int i;
    int N,ktype,kopt;

    double mf[NMIN];
    int type;

    if(argc < 2 || argc > 3) {
        usage();
        exit(1);
    }
    strcpy(infile,argv[1]);
    fprintf(stdout,"%% Infile: %s\n",infile);

    if((fin = fopen(infile, "r")) == NULL) {
        fprintf(stderr, "Cannot open %s for read!!!! DIE!!!!!!\n", infile);
        return(-1);
    }

    numrocks(fin,&N,&ktype);
    for (i = 0; i < N; i++) {
        kopt = ktype;
        if (loadrock(fin, mf, &type, &kopt) != 0) {
            fprintf(stderr,"Error loading rock %i.\n",i+1);
            return(-1);
        }
        outputmineralogy(i,mf);
        outputheader(kopt);

        if (type == 0) {
            stepcompute(fin, mf, kopt);
        } else {
            PTcompute(fin, mf, kopt);
        }
        fprintf(stdout,"NaN     NaN  NaN        NaN     NaN     NaN    NaN\n");
        fprintf(stdout,"\n");
    }

    fclose(fin);

    return(0);
}
