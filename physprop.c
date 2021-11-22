#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "physprop.h"

/* Mineral references and notes:                                        *
 *                                                                      *
 *  1 from Jackson et al. [Eur. J. Mineral., 2003]                      *
 *  2 average model of Hugh-Jones [Am. Mineral., 1997] and Yang &       *
 *    Ghose [Phys. Chem. Mineral., 1994]                                *
 *  3 from Sueda et al.[PEPI, 2008]                                     *
 *  4 KT and K' for Enstatite from Angel & Jackson [Am. Mineral.,       *
 *    2002], which shows results fairly consistent with other           *
 *    measurements.  The second order quadratic fit to the data by      *
 *    Flesch et al. [Am. Mineral., 1998] yields a K'' of -0.0064        *
 *    GPa^-1.                                                           *
 *  5 K' and gammaT from Yang & Ghose [Phys. Chem. Mineral., 1994].     *
 *    change is also made there.                                        *
 *  6 Computed from adiabatic bulk modulii in Schutt & Lescher [JGR,    *
 *    2006].  The Gruenesion parameter is estimated from pyrope and     *
 *    the Anderson-Gruenesion parameter also from pyrope.               *
 *  7 All values from Sueda et al.[PEPI, 2008].  Although other         *
 *    published Gruenesion parameters for spinel are ~1.3.              *
 *  8 from Angel [Contrib. Mineral. Petrol., 2004], Downs [Am.          *
 *    Mineralogist, 1994], the temperature derivative of the bulk       *
 *    modulus actually comes from that of Nepheline (Na3KAl4Si4O16),    *
 *    a feldspathoid of similar composition to Albite.  This value is   *
 *    an order of magnitude lower than most common mantle silicates     *
 *    [Bass, AGU, Handbook of Physical Constants, 1995].                *
 *    I would guess that Anorthite has a temperature derivative more    *
 *    similar to other silicates with higher bulk modulii hence the     *
 *    estimate at -0.01, but still lower in keeping with feldspar(?).   *
 *    Basically pulled the number out of my ass, but couldn't find a    *
 *    published value anywhere... maybe someday I'll learn to compute   *
 *    a more reasonable value.                                          *
 *    Grueneison parameters for plagioclase phases from Hacker et al.   *
 *    [JGR, 2003].                                                      *
 *  9 from Hacker et al. [JGR, 2003]                                    *
 * 10 reformulated from Holland and Powell [J. Metamorph. Petrol.,      *
 *    1998]                                                             *
 * 11 Berman & Aranovich [Contrib. Mineral.  Petrol. 1996]             *
 * 12 Berman & Brown [Contrib. Mineral. Petrol. 1985]                   *
 * 13 reformulated consants from Holland & Powell [1998]                *
 * 14 [AGU Reference Desk, Handbook of Physical Properties, 1995]       *
 * 15 from Stixrude & Lithgow-Bertelloni [Geophys. J. Int., 2005] the   *
 *    delta_T values are estimated by delta_T = gammaT + K' (see        *
 *    Hacker & Abers, [JGR, 2003] for approximation)                    *
 * 16 values from Pawley et al. [Am. Mineral., 2002] and values for     *
 *    clinochlore mean from Welch & Crichton [Eur. J. Mineral., 2002]   *
 *    and Pawley et al. [Am. Mineral., 2002]                            *
 * 17 hornblende is most closely related to pargasite and in the        *
 *    absence of additional specific heat capacity data, the values for *
 *    paragasite are used.                                              * 
 * 18 Muscovite values from Comodi and Zanazzi [Phys. Chem. Minerals,   *
 *    1995].                                                            *
      Annite values for KT and K' and alpha from Pavese et al. [Am.     *
 *    Min., 2007].  Phlogopite values are average of Hazen and Finger   *
 *    [Am. Min. 1978], Comodi et al. [Eur. J. Min. 2004], Pavese et al. *
 *    [Eur. J. Min. 2003] and therefore given the same values.          *
 * 19 Holland and Powell [J. Metamorphic Geol., 1998]                   *
 * 20 Values from Robertson & Hemingway [USGS OFR 95-622, 1995] switch  *
 *    between aqz and bqz at 844 K                                      *
 * 21 Values from Hacker et al. [JGR, 2003]                             *
 * 22 Values for bulk modulus from Usui et al. [Japan Geoscience Union, *
 *    Mtg., 2003].  K' is assumed.                                      *
 * 23 Values from Ross and Reeder [Am. Mineralogist 1992] pressure      *
 *    derivative is assumed to be 4.                                    *
 * 24 Isaak et al. [Phys. Chem. Minerals, 1998].                        *
 * 25 Yang et al. [Phys. Chem. Minerals, 1997].                         *
 * 26 Estimated                                                         *
 * 27 Estimated form Comodi & Zanazzi [Phys. Chem. Minerals, 1997]      *
 *                                                                      * 
 * Mineral abbreviations and index number                               *
 *                                                                      *
 *  0:  1: qtz-aqz      alpha-quartz      SiO2                          *
 *  1:  2: qtz-bqz      beta-quartz       SiO2                          *
 *  2:  3: or-or        orthoclase        KAlSi3O8                      *
 *  3:  3: or-san       sanadine          KAlSi3O8                      *
 *  4:  4: plag-ab      albite            NaAl3Si3O8                    *
 *  5:  4: plag-an      anorthite         CaAl2Si2O8                    *
 *  6:  5: musc-ms      muscovite         KAl3Si3O10(OH)2               *
 *  7:  5: musc-pg      paragonite        NaAl3Si3O10(OH)2              *
 *  8:  6: bt-phlog     phlogopite        KMg3AlSi3O10(FOH)2            *
 *  9:  6: bt-ann       annite            KFe3AlSi3O10(FOH)2            *
 * 10:  7: ep-ep        Fe-epidote        Ca2FeAl2Si3O12(OH)            *
 * 11:  8: chl-clin     clinochlore       Mg5Al(AlSi3)O10(OH)8          *
 * 12:  8: chl-daph     daphnite          Fe5Al(AlSi3)O10(OH)8          *
 * 13:  9: amph-hbl     hornblende        Ca2(MgFeAl)5(AlSi)8O22(OH)2   *
 * 14:  9: amph-tr      tremolite         Ca2Mg5Si8O22(OH)2             *
 * 15:  9: amph-fact    ferro-actinolite  Ca2Fe5Si8O22(OH)2             *
 * 16:  9: amph-parg    pargasite         NaCa2Mg4Al(Al2Si6)O22(OH)2    *
 * 17: 10: cpx-di       diopside          CaMgSi2O6                     *
 * 18: 10: cpx-hd       hedenbergite      CaFeSi2O6                     *
 * 19: 10: cpx-jd       jadite            NaAlSi2O6                     *
 * 20: 11: wo-wo        wollastonite      CaSiO3                        *
 * 21: 12: opx-en       enstatite         MgSiO3                        *
 * 22: 12: opx-fs       ferrosillite      FeSiO3                        *
 * 23: 13: ol-fo        forsterite        MgSiO4                        *
 * 24: 13: ol-fa        fayalite          FeSiO4                        *
 * 25: 14: sp-sp        spinel            MgAl2O4                       *
 * 26: 14: sp-mt        magnetite         Fe3O4                         *
 * 27: 14: sp-herc      hercynite         FeAl2O4                       *
 * 28: 15: gt-py        pyrope            Mg3Al2Si3O12                  *
 * 29: 15: gt-al        almandine         Fe3Al2Si3O12                  *
 * 30: 15: gt-gr        grossular         Ca3Al2Si3O12                  *
 * 31: 15: gt-spess     spessartine       Mn3Al2Si3O12                  *
 * 32: 16: carb-arag    aragonite         CaCO3                         *
 * 33: 16: carb-cc      calcite           CaCO3                         *
 * 34: 16: carb-do      dolomite          CaMg(CO3)2                    *
 * 35: 17: als-and      andalucite        Al2SiO5                       *
 * 36: 17: als-ky       kyanite           Al2SiO5                       *
 * 37: 17: als-sil      sillimanite       Al2SiO5                       *
 * 38: 18: rut-rut      rutile            TiO2                          */

/* indicies for each mineral mode [start,end]                           */
static const int ind[NMODE][2] = {
    { 0,  0},  // a-quartz
    { 1,  1},  // b-quartz
    { 2,  2},  // orthoclase
    { 4,  5},  // plagioclase
    { 6,  7},  // muscovite
    { 8,  9},  // biotite
    {10, 10},  // epidote
    {11, 12},  // chlorite
    {13, 16},  // amphibole
    {17, 19},  // clinopyroxene
    {20, 20},  // wollastonite
    {21, 22},  // orthopyroxene
    {23, 24},  // olivine
    {25, 27},  // spinel
    {28, 31},  // garnet
    {32, 34},  // carbonate
    {35, 37},  // aluminosilicate
    {38, 38},  // rutile
    { 3,  3}}; // sanidine

/* formula weights for minerals in table                                */
static const double min_wt[NMIN] = {
     60.08,  // alpha-quartz
     60.08,  // beta-quartz
    278.34,  // orthoclase
    278.34,  // sanadine
    262.22,  // albite
    278.21,  // anorthite
    398.31,  // muscovite
    382.20,  // paragonite
    417.29,  // phlogopite
    511.89,  // annite
    512.09,  // Fe-epidote
    555.83,  // clinochlore
    713.51,  // daphnite
    864.70,  // hornblende
    812.41,  // tremolite
    970.08,  // ferro-actinolite
    835.86,  // pargasite
    216.56,  // diopside
    248.10,  // hedenbergite
    202.14,  // jadite
    116.16,  // wollastonite
    100.40,  // enstatite
    131.93,  // ferrosillite
    140.71,  // forsterite
    203.78,  // fayalite
    142.27,  // spinel
    231.54,  // magnetite
    173.81,  // hercynite
    403.15,  // pyrope
    497.75,  // almandine
    450.45,  // grossular
    495.03,  // spessartine
    100.09,  // aragonite
    100.09,  // calcite
    184.41,  // dolomite
    162.05,  // andalucite
    162.05,  // kyanite
    162.05,  // sillimanite
     79.90}; // rutile

/* Density at P = 0, T = 298 K                                          *
 * Values from Hacker et al. [JGR, 2004]                                */
static const double rho298[NMIN] = {
    2648,  // alpha-quartz
    2530,  // beta-quartz
    2555,  // orthoclase
    2553,  // sanadine
    2620,  // albite
    2760,  // anorthite
    2828,  // muscovite
    2780,  // paragonite
    2788,  // phlogopite
    3317,  // annite
    3498,  // Fe-epidote
    2635,  // clinochlore
    3343,  // daphnite
    3248,  // hornblende
    2979,  // tremolite
    3430,  // ferro-actinolite
    3074,  // pargasite
    3272,  // diopside
    3651,  // hedenbergite
    3346,  // jadite
    2910,  // wollastonite
    3206,  // enstatite
    4003,  // ferrosillite
    3222,  // forsterite
    4400,  // fayalite
    3575,  // spinel
    5201,  // magnetite
    4264,  // hercynite
    3565,  // pyrope
    4324,  // almandine
    3593,  // grossular
    4194,  // spessartine
    2931,  // aragonite
    2713,  // calcite
    2830,  // dolomite
    3150,  // andalucite
    3670,  // kyanite
    3249,  // sillimanite
    4242}; // rutile

/* Expansivity coefficients Fei [AGU, Handbook of Physical Constants,   *
 * 1995], functional form:                                              *
 *    alpha = a0 + a1*T + a2*T^-2                                       *
 *                                                                      *
 *    a0          a1          a2        Mineral               Note/Ref  */
static const double a[NMIN][3] = {
    { 1.4170e-5,  9.6581e-8, -1.6973},  // alpha-quartz
    {-0.4400e-5,  0.0000e-8,  0.0000},  // beta-quartz
    { 3.4000e-5,  0.0000e-8,  0.0000},  // orthoclase
    { 3.4000e-5,  0.0000e-8,  0.0000},  // sanadine
    { 1.9801e-5,  1.0065e-8, -0.9760},  // albite
    { 1.2491e-5, -0.0162e-8,  0.0161},  // anorthite
    { 3.5370e-5,  0.0000e-8,  0.0000},  // muscovite
    { 3.5370e-5,  0.0000e-8,  0.0000},  // paragonite
    { 5.8000e-5,  0.0000e-8,  0.0000},  // phlogopite
    { 3.7000e-5,  0.0000e-8,  0.0000},  // annite
    { 5.1000e-5,  0.0000e-8,  0.0000},  // Fe-epidote
    { 2.5000e-5,  0.0000e-8,  0.0000},  // clinochlore
    { 2.5000e-5,  0.0000e-8,  0.0000},  // daphnite
    { 2.0750e-5,  1.0270e-8,  0.0000},  // hornblende
    { 3.1310e-5,  0.0000e-8,  0.0000},  // tremolite
    { 3.1310e-5,  0.0000e-8,  0.0000},  // ferro-actinolite
    { 3.1310e-5,  0.0000e-8,  0.0000},  // pargasite
    { 3.3300e-5,  0.0000e-8,  0.0000},  // diopside
    { 2.9800e-5,  0.0000e-8,  0.0000},  // hedenbergite
    { 2.4700e-5,  0.0000e-8,  0.0000},  // jadite
    { 2.4740e-5,  1.0096e-8,  0.0000},  // wollastonite
    { 2.9720e-5,  0.5711e-8,  0.0000},  // enstatite
    { 2.8750e-5,  0.0000e-8,  0.0000},  // ferrosillite
    { 2.8540e-5,  1.0080e-8, -0.3842},  // forsterite
    { 2.3860e-5,  1.1530e-8, -0.0518},  // fayalite
    { 1.9600e-5,  1.6400e-8,  0.0000},  // spinel
    { 3.6775e-5,  1.2225e-8, -1.0268},  // magnetite
    { 0.9700e-5,  1.9392e-8,  0.0000},  // hercynite
    { 2.3110e-5,  0.5956e-8, -0.4538},  // pyrope
    { 1.7760e-5,  1.2140e-8, -0.5071},  // almandine
    { 1.9510e-5,  0.8089e-8, -0.4972},  // grossular
    { 2.9270e-5,  0.2726e-8, -1.1560},  // spessartine
    { 6.2210e-5,  0.0000e-8,  0.0000},  // aragonite
    { 0.7130e-5,  3.3941e-8, -1.2140},  // calcite
    { 1.9280e-5,  3.1703e-8, -0.5393},  // dolomite
    { 2.1810e-5,  0.3261e-8,  0.0000},  // andalucite
    { 2.5050e-5,  0.0000e-8,  0.0000},  // kyanite
    { 1.2600e-5,  0.2314e-8,  0.0000},  // sillimanite
    { 1.6230e-5,  1.6371e-8, -0.4799}}; // rutile

/* Elastic constants from Hacker et al. [JGR, 2003] and references      *
 * therein (unless noted otherwise)                                     *
 * functional form:                                                     *
 *    KT = KT,0 + K'*P + 0.5*K''*P^2 + dKs/dT*T                         *
 *                                                                      *
 *    KT    K'    K''   dKs/dT                                          */
static const double K[NMIN][4] = {
    { 37.1, 5.99,  0.0000,  0.0000},  // alpha-quartz
    { 57,   4,     0.0000,  0.0000},  // beta-quartz
    { 58.3, 4,     0.0000,  0.0000},  // orthoclase
    { 67,   4,     0.0000,  0.0000},  // sanadine
    { 53.8, 6,    -0.0600, -0.0037},  // albite
    { 82.5, 3.2,  -0.0500, -0.0100},  // anorthite
    { 50,   3,     0.0000,  0.0000},  // muscovite
    { 60,   2.9,   0.0000,  0.0000},  // paragonite
    { 54,   7.8,  -0.5953,  0.0000},  // phlogopite
    { 49,   8.1,  -0.5047,  0.0000},  // annite
    {162,   4,     0.0000,  0.0000},  // Fe-epidote
    { 85,   3.3,   0.0000, -0.0130},  // clinochlore
    { 86.9, 4,     0.0000, -0.0130},  // daphnite
    { 94,   4,     0.0000,  0.0000},  // hornblende
    { 85,   4,     0.0000,  0.0000},  // tremolite
    { 76,   4,     0.0000,  0.0000},  // ferro-actinolite
    { 91.2, 4,     0.0000,  0.0000},  // pargasite
    {113,   4.8,   0.0000, -0.0205},  // diopside
    {119,   4,     0.0000, -0.0205},  // hedenbergite
    {139,   5,     0.0000, -0.0165},  // jadite
    {107,   4,     0.0000,  0.0000},  // wollastonite
    {105.8, 8.5,  -1.6000, -0.0370},  // enstatite
    {100,   8.8,  -1.6000, -0.0370},  // ferrosillite
    {128,   4.2,   0.0000, -0.0170},  // forsterite
    {135,   4.2,   0.0000, -0.0204},  // fayalite
    {205,   4.1,   0.0000, -0.0300},  // spinel
    {181,   5.5,   0.0000,  0.0000},  // magnetite
    {209,   4,     0.0000,  0.0000},  // hercynite
    {173,   5,     0.0000, -0.0191},  // pyrope
    {174,   6,     0.0000, -0.0204},  // almandine
    {168,   5.5,   0.0000, -0.0148},  // grossular
    {182,   5.6,   0.0000, -0.0191},  // spessartine
    { 46,   4,     0.0000,  0.0000},  // aragonite
    { 73.5, 4,     0.0000,  0.0000},  // calcite
    { 94,   4,     0.0000,  0.0000},  // dolomite
    {151,   4,     0.0000,  0.0000},  // andalucite
    {156,   5.6,   0.0000,  0.0000},  // kyanite
    {171,   4,     0.0000,  0.0000},  // sillimanite
    {211.3, 6.8,   0.0000, -0.0500}}; // rutile

/* Gruenesion parameter                                                 *
 * From Hacker et al. [JGR, 2003]                                       */
static const double gammaT[NMIN] = {
    0.7,   // alpha-quartz
    0.1,   // beta-quartz
    0.4,   // orthoclase
    0.4,   // sanadine
    0.6,   // albite
    0.6,   // anorthite
    0.5,   // muscovite
    0.5,   // paragonite
    0.6,   // phlogopite
    0.6,   // annite
    1.1,   // Fe-epidote
    0.3,   // clinochlore
    0.3,   // daphnite
    1.1,   // hornblende
    0.74,  // tremolite
    0.73,  // ferro-actinolite
    0.84,  // pargasite
    1,     // diopside
    1.5,   // hedenbergite
    1,     // jadite
    1,     // wollastonite
    1.01,  // enstatite
    0.88,  // ferrosillite
    0.99,  // forsterite
    1.06,  // fayalite
    1.73,  // spinel
    2,     // magnetite
    1.2,   // hercynite
    1.29,  // pyrope
    1.29,  // almandine
    1.38,  // grossular
    1.29,  // spessartine
    1.4,   // aragonite
    0.7,   // calcite
    0.7,   // dolomite
    1,     // andalucite
    1,     // kyanite
    0.5,   // sillimanite
    1.62}; // rutile

/* Anderson-Gruenesion parameter                                        *
 * From Hacker et al. [JGR, 2003]                                       */
static const double deltaT[NMIN] = {
    8.42,  // alpha-quartz
    4.11,  // beta-quartz
    4.44,  // orthoclase
    4.44,  // sanadine
    6.57,  // albite
    3.47,  // anorthite
    7.42,  // muscovite
    3.4,   // paragonite
    4.55,  // phlogopite
    4.55,  // annite
    5.11,  // Fe-epidote
    4.3,   // clinochlore
    4.3,   // daphnite
    5.1,   // hornblende
    4.74,  // tremolite
    4.73,  // ferro-actinolite
    4.84,  // pargasite
    6.04,  // diopside
    5.21,  // hedenbergite
    4.99,  // jadite
    5,     // wollastonite
    9.39,  // enstatite
    9.05,  // ferrosillite
    5.19,  // forsterite
    5.26,  // fayalite
    6.5,   // spinel
    7.45,  // magnetite
    5.19,  // hercynite
    5.3,   // pyrope
    5.52,  // almandine
    4.57,  // grossular
    5.3,   // spessartine
    5.36,  // aragonite
    4.69,  // calcite
    4.7,   // dolomite
    5,     // andalucite
    6.56,  // kyanite
    4.5,   // sillimanite
    9.54}; // rutile

/* Fitting parameters                                                   *
 * functional form:                                                     *
 *    Cp = k0 + k1*T^-0.5 + k2*T^-2 + k3*T^-3                           */
static const double c[NMIN][4] = {
    {1.7699e3, -23.4900e3,  45.9470e6, -7.2006e9},  // alpha-quartz
    {0.9647e3,  -0.0004e3,  30.5460e6, -0.0018e9},  // beta-quartz
    {1.4445e3,  -9.5303e3, -26.3950e6,  3.6781e9},  // orthoclase
    {1.4445e3,  -9.5303e3, -26.3950e6,  3.6781e9},  // sanadine
    {1.5012e3,  -9.2116e3, -30.0990e6,  4.0829e9},  // albite
    {1.5793e3,  -4.7252e3,   0.0000e6, -1.1395e9},  // anorthite
    {1.6680e3, -10.9940e3, -36.7900e6,  5.0614e9},  // muscovite
    {1.7178e3, -10.4970e3, -51.4220e6,  8.3960e9},  // paragonite
    {1.4354e3,  -4.0508e3, -61.2800e6,  8.9950e9},  // phlogopite
    {1.2777e3,  -5.6610e3, -42.8130e6,  6.9201e9},  // annite
    {1.3009e3,  -7.6344e3,   8.5799e6, -6.2503e9},  // Fe-epidote
    {2.1726e3, -19.7390e3,  -2.9739e6, -1.7267e9},  // clinochlore
    {1.8204e3, -18.1750e3,   6.0436e6, -1.8045e9},  // daphnite
    {1.6559e3, -13.1260e3,   1.5163e6, -2.6059e9},  // hornblende
    {1.5725e3, -10.7350e3, -11.3070e6, -0.4465e9},  // tremolite
    {1.4695e3, -12.6660e3,   9.6112e6, -2.9282e9},  // ferro-actinolite
    {1.6559e3, -13.1260e3,   1.5163e6, -2.6059e9},  // pargasite
    {1.4103e3,  -7.4109e3, -33.0900e6,  4.2567e9},  // diopside
    {1.3712e3, -10.6680e3, -32.7290e6, -0.0042e9},  // hedenbergite
    {1.5400e3,  -9.9194e3, -26.4680e6,  3.2778e9},  // jadite
    {1.3714e3,  -9.2579e3,  -8.3273e6,  0.0709e9},  // wollastonite
    {1.6592e3, -11.9590e3, -22.6160e6,  2.7805e9},  // enstatite
    {1.3204e3, -10.5590e3,  -3.4440e6, -0.2858e9},  // ferrosillite
    {1.6572e3, -12.8040e3,   0.0000e6, -1.9042e9},  // forsterite
    {1.2366e3,  -9.8820e3,   0.0000e6, -0.3052e9},  // fayalite
    {1.7197e3, -14.0860e3,   0.0000e6,  0.0000e9},  // spinel
    {0.9894e3,  -3.0731e3, -27.8960e6,  3.1616e9},  // magnetite
    {1.4865e3, -11.5590e3, -15.9530e6,  3.1429e9},  // hercynite
    {1.4657e3,  -7.0120e3, -33.0420e6,  3.1262e9},  // pyrope
    {1.2485e3,  -6.6050e3, -30.2980e6,  4.4437e9},  // almandine
    {1.1531e3,  -0.1401e3, -62.1740e6,  7.7938e9},  // grossular
    {1.1660e3,  -5.1337e3, -17.2090e6, -0.3270e9},  // spessartine
    {1.6407e3, -15.1610e3,  22.1160e6, -5.1056e9},  // aragonite
    {1.7799e3, -17.1680e3,  -7.6997e6,  3.0985e9},  // calcite
    {1.8229e3, -15.2560e3, -16.7360e6,  2.7028e9},  // dolomite
    {1.5227e3,  -8.6527e3, -37.3920e6,  4.1310e9},  // andalucite
    {1.5202e3,  -8.3661e3, -40.3460e6,  4.4671e9},  // kyanite
    {1.5316e3,  -9.2260e3, -35.2810e6,  4.3266e9},  // sillimanite
    {1.3212e3, -12.9190e3,  19.2410e6, -2.6362e9}}; // rutile

/* Conductivity References:                                         *
 *  1 value from Harrell [U. WA, Ph.D. Thesis, 2002]                *
 *  2 from Xu et al. [PEPI, 2004]                                   *
 *  3 fit to model for augite rather than diopside by Hofmeister    *
 *    & Pertermann [Eur. J. Mineral., 2008] because their values    *
 *    for diopside seem unrealistically high, especially at low     *
 *    temperatures                                                  *
 *  4 determined from 2nd order fit to data in Geisting &           *
 *    Hofmeister [Phys.  Rev. B, 2002]                              *
 *  5 fit to data in Hofmeister [Am. Mineral., 2001]                *
 *  6 fit to data in Petrunin et al. [Izvestiya, 2004]              *
 *  7 fit to thermal diffusivity in Hofmeister [Am. Mineral., 2007] *
 *    using heat capacity and density (thermal expansion) from      *
 *    Holland and Powell [J. Metamorph. Petrol., 1998]              *
 *  8 mean of values from Harrell [U. WA, Ph.D. Thesis, 2002],      *
 *    Horai [JGR, 1971], and Diment [USGS, OFR, 1989] for augite    * 
 *    and diopside                                                  *
 *  9 mean of values from Clauser and Hunges [1995] and references  *
 *    therein                                                       *
 * 10 mean of amphiboles from Diment & Pratt [USGS OFR, 1988] and   *
 *    Horai [JGR, 1971]                                             */

/* temperature exponent */
static const double n[NMODE] = {
     1.48,  //          a-quartz Branlund & Hofmeister [2007]
    -1.00,  //          b-quartz Branlund & Hofmeister [2007]
//     1.32,  //          a-quartz Hoefer & Schilling [2002]
//    -3.35,  //          b-quartz Hoefer & Schilling [2002]
    -0.27,  // f(or)    orthoclase
    -0.21,  // f(plag)  plagioclase
     1.54,  // f(an)    muscovite
     1.54,  // f(phlog) biotite
     0.5,   // f(ep)    epidote
     1.54,  // f(clin)  chlorite
     0.5,   //          amphibole (hbl, tr, fact, parg)
     0.54,  // f(di)    clinopyroxene
     0.5,   // f(wo)    wollastonite
     0.3,   // f(en)    orthopyroxene
     0.49,  // f(fo)    olivine
     1.24,  // f(sp)    spinel
     0.37,  // f(py)    garnet
     0.78,  // f(cc)    carbonate
     1,     //          andalucite (sill, ky, and)
     1.5,   // f(rt)    rutile
    -0.30}; // f(san)   sanidine

/* Exponent references                                              *
 * 1 determined in a fit to Harrell [Ph.D. Thesis, Univ. Wash.,     *
 *   2002]                                                          *
 * 2 from Xu et al. [PEPI, 2004]                                    *
 * 3 fit to data in Hofmeister & Pertermann [Eur. J. Mineral.,      *
 *   2008]                                                          *
 * 4 fit to data in Osako et al. [PEPI, 2004]                       *
 * 5 fit to data in Hofmeister [Am. Mineral., 2001]                 *
 * 6 fit to data in Petrunin et al. [Izvestiya, 2004]               *
 * 7 see 7 above                                                    *
 * 8 etimated.                                                      */

/* Fitting constants */
static const double kl[NMODE][4] = {
    { 8.79,     0,      0,      0},       //         a-quartz Branlund & Hofmeister [2007]
    { 0.99,     0,      0,      0},       //         b-quartz Branlund & Hofmeister [2007]
//    { 8.85,     0,      0,      0},       //         a-quartz Hoefer & Schilling [2002]
//    { 0.0456,   0,      0,      0},       //         b-quartz Hoefer & Schilling [2002]
    { 1.79,     0,      0,      0},       // f(or)   orthoclase
    { 2.2,     -2.18,   1.9,    0},       // f(plag) plagioclase
    { 2.80,     0,      0,      0},       // f(an)   muscovite  Clauser & Hunges [1995]
    { 2.27,     0,      0,      0},       // f(musc) biotite  Clauser & Hunges [1995]
    { 3.08,     0,      0,      0},       // f(ep)   epidote
    { 4.35,     0,      0,      0},       // f(clin) chlorite
    { 2.6485,   4.0654, 3.7834, 1.7066},  //         amphibole (hbl, tr, fact, parg)
    { 4.25,     0,      0,      0},       // f(di)   clinopyroxene
    { 4.03,     0,      0,      0},       // f(wo)   wollastonite
    { 3.37,     0,      0,      0},       // f(en)   orthopyroxene
    { 3.09,    -1.17,   3.35,   0},       // f(fo)   olivine
    {11.94,     5.47,   8.71,   0},       // f(sp)   spinel
    { 4.97,    -7.42,   8.45,   0},       // f(py)   garnet
    { 3.46,     0,      0,      0},       // f(cc)   carbonate
    {10.32,    12.87,   6.85,   0},       //         andalucite (sill, ky, and)
    { 4.95,     0,      0,      0},       // f(rt)   rutile
    { 1.86,     0,      0,      0}};      // f(san)  sanidine

    /* amphibole parameters - A*m */
    //A = {{0,2,3,1,2,7},
    //     {0,2,5,0,0,8},
    //     {0,2,0,5,0,8},
    //     {1,2,4,0,3,6}};
    //m = {-0.1760, 0.2008, 0.1070, 0.0506, -0.4312, 0.3911};

// Radiative conductivity using O(T^3)
// krad = kr[0] + kr[1]*T + kr[2]*T^2 + kr[3]*T^3
//
//static const double kr[NMODE][4] = {
//    {5.1832, -0.0355,    7.5468e-5,  -4.5357e-8 },  //         a-quartz
//    {0,       0,         4.6418e-6,   7.5449e-10},  //         b-quartz
//    {0.2521,  0,         2.0765e-6,  -6.8230e-10},  // f(or)   orthoclase
//    {0,       0,         0,           0         },  // f(plag) plagioclase
//    {0.2226, -0.0015,    3.0404e-6,  -1.6450e-9 },  // f(an)   muscovite
//    {0.2226, -0.0015,    3.0404e-6,  -1.6450e-9 },  // f(bt)   phlogopite
//    {0,       0,         0,           0         },  // f(ep)   epidote
//    {0.2226, -0.0015,    3.0404e-6,  -1.6450e-9 },  // f(clin) chlorite
//    {0,       0,         0,           0         },  //         amphibole
//    {0.0729, -6.4079e-4, 1.3520e-6,  -3.3319e-10},  // f(di)   clinopyroxene
//    {0.0729, -6.4079e-4, 1.3520e-6,  -3.3319e-10},  // f(wo)   wollastonite
//    {0.0729, -6.4079e-4, 1.3520e-6,  -3.3319e-10},  // f(en)   orthopyroxene     *estimated from Harrell [2002] (data from Katsura [GJI, 1995])
//    {0.0729, -6.4079e-4, 1.3520e-6,  -3.3319e-10},  // f(fo)   olivine           Harrell [2002] (1.12*Katsura [GJI, 1995])
//    //{0.0175, -1.037e-4,  2.245e-7,   -3.4047e-11},  // f(fo)   olivine           erroneous? Hofmeister [1999]
//    //{-0.7670, 2.70e-03, -2.6384e-06,  9.7755e-10},  // f(en)   orthopyroxene     Schatz & Simmons [1972]
//    //{-0.778,  0.0016,    3.0973e-6,  -1.8627e-9 },  // f(fo)   olivine           Shankland et al. [1979] and references therein
//    {0,       0,         0,           0         },  // f(sp)   spinel
//    {0,       0,         0,           0         },  // f(py)   garnet
//    {1.5831, -0.0108,    2.2516e-5,  -1.3987e-8 },  // f(cc)   carbonate
//    {0,       0,         0,           0         },  //         andalucite
//    {0,       0,         0,           0         }}; // f(rt)   rutile


// Radiative conductivity using O(T^3)
// krad = kr[0](1 + erf(kr[1]*(T - kr[2])))
//
static const double kr[NMODE][4] = {
    {0.8107, 6.553e-3,  558},  //         a-quartz  Branlund & Hofmeister [2007]
    {0.8107, 6.553e-3,  558},  //         b-quartz  Branlund & Hofmeister [2007]
//    {0.8102, 6.916e-3,  592},  //         a-quartz Hoefer & Schilling [2002]
//    {0.8102, 6.916e-3,  592},  //         b-quartz Hoefer & Schilling [2002]
    {0.1949, 3.782e-3,  864},  // f(or)   orthoclase
    {0.0000, 0.000e-3,    0},  // f(plag) plagioclase
    {0.1166, 2.749e-3,  520},  // f(an)   muscovite
    {0.1166, 2.749e-3,  520},  // f(bt)   phlogopite
    {0.0000, 0.000e-3,    0},  // f(ep)   epidote
    {0.1166, 2.749e-3,  520},  // f(clin) chlorite
//    {0.3797, 1.164e-3, 1141},  //         amphibole
//    {0.3797, 1.164e-3, 1141},  // f(di)   clinopyroxene
//    {0.3797, 1.164e-3, 1141},  // f(wo)   wollastonite
//    {0.3797, 1.164e-3, 1141},  // f(en)   orthopyroxene
//    {0.5612, 2.746e-3, 1150},  // f(fo)   olivine
    {0.1725, 3.900e-3,  762},  //         amphibole
    {0.1725, 3.900e-3,  762},  // f(di)   clinopyroxene
    {0.1725, 3.900e-3,  762},  // f(wo)   wollastonite
    {0.1725, 3.900e-3,  762},  // f(en)   orthopyroxene
    {0.1725, 3.900e-3,  762},  // f(fo)   olivine
    {0.0000, 0.000e-3,    0},  // f(sp)   spinel
    {0.0000, 0.000e-3,    0},  // f(py)   garnet
    {0.1288, 1.067e-3,  550},  // f(cc)   carbonate
    {0.0000, 0.000e-3,    0},  //         andalucite
    {0.0000, 0.000e-3,    0},  // f(rt)   rutile
    {1.9890, 2.738e-3, 1089}}; // f(san)  sanidine

/* *not enough data/range of compositions analyized to make an          *
 *  accurate assessment of mixtures.  Hence averages of measured        *
 *  values are used for these phases. The values here are k298 values   *
 *  from Harrell [U. WA, Ph.D. Thesis, 2002] for opx and cpx.  The      *
 *  values are not end-member compositions and are both significantly   *
 *  lower than average values of En(86-100): 4.23 W/m/K [Horai, JGR,    *
 *  1971] and Di: 5.16 W/m/K Hofmeister & Pertermann [Eur. J. Mineral., *
 *  2008] actually for their augite.                                    *
 * ^amphiboles are treated specially because they are a complex series  *
 *  solution with multiple numerous cation and anion substitutions      *
 * +mixture of Al-spinel (X = 1) with magnetite (X = 0)                 *
 * ~Gr and Spess not included.                                          *
 *  opx, cpx, ol, amph, and wo all set to radiaive conductivity         *
 *  from difference estimate of lherzolite (Gibert et al. 2003) model.  */

double ref_density(int i, double Tref) {

    double rhor;

    rhor = rho298[i]*exp(-ithermexp(i,Tref,298));

    return(rhor);
}


/* ADJ_MF - Fix mass fractions for transitions                          */
void adj_mf(double T, double mf[], double *cmf) {

    int i;

    for (i = 0; i < NMIN; i++) {
        cmf[i] = mf[i];
    }

    /* alpha-beta quartz transition */
    if (T >= 845) {
        cmf[1] = mf[0];
        cmf[0] = 0;
    }
    /* other transitions to be added later */
}


/* ADJ_ALPHA - Adjust thermal expansion.                                *
 *                                                                      *
 * Took the easy way out, numerical fit to curves in Korenaga [JGR,     *
 * 2007].  Curves are similar, but not exact, maybe could be better but *
 * the exact numerical solution is complex... figure it out later.      */
double adj_alpha(double t, double T) {

    double u,R;

    //R = 0.85;
    //return(R);

    /* unstable if t < 1, thus for all t < 1, t = 1 */
    if (t > 1) {
        u = T - 552.1*pow(t,-0.25);
    } else {
        u = T - 552.1;
    }
    R = 0.5707 + u*(4.777e-4 - u*(5.823e-7 - u*5.403e-10));

    /* effective expansivity is bounded between the instantaneous       *
     * relaxation R = 1 and infinite relaxation time R = 0.55.          */
    if (R > 1) {
        R = 1;
    } else if (R < 0.55) {
        R = 0.55;
    }

    return(R);
}


/* COMPUTEDENSITY - Computes density using effective thermal        *
 * expansion.                                                           */
void computedensity(int aflag, double t, double Tref, double P, double T,
  double mf[], double *rho, double *alpha, double *vf) {

    int i;
    double alphaPT[NMIN],rhoPT[NMIN];
    double w[NMODE],vmode[NMODE];
    double R;

    /* ratio effective expansivity/expansivity */
    if (aflag) {
        R = adj_alpha(t,T);
    } else {
        R = 1;
    }

    /* density and thermal expansion of minerals */
    for (i = 0; i < NMIN; i++) {
        if (mf[i] != 0) {
            densityPT(i,Tref,P,T,&alphaPT[i],&rhoPT[i],R);
            //printf("%lg %lg\n",alphaPT[i],rhoPT[i]);
        } else {
            rhoPT[i] = 0;
            alphaPT[i] = 0;
        }
    }

    /* compute volume fraction from mass fraction */
    mf2vf(mf,vf,rhoPT);
    //for (i = 0; i < NMIN; i++) {
    //    fprintf(stdout," %8.4lf  %8.4lf\n",mf[i]*100,vf[i]*100);
    //}
    //fprintf(stdout,"\n");
    //exit(0);

    /* compute mode weights                       */
    modefrac(vf,w);

    /* rock expansivity */
    modevalue(alphaPT,vf,vmode);
    //for (i = 0; i < NMODE; i++) {
    //    printf("%lg: %lg   ",w[i],vmode[i]);
    //}
    //printf("\n");
    *alpha = wtmean(w,vmode);

    /* rock density */
    modevalue(rhoPT,vf,vmode);
    *rho = wtmean(w,vmode);
}


/* ******************************************************************** *
 * COMPUTEHEATCAP - Compute rock heat capacity.                         *
 *                                                                      *
 * Inputs:                                                              *
 *    P         pressure                                                *
 *    T         temperature                                             *
 *    vf[NMIN]  volume fraction of mineral assemblage                   *
 *   *Cp        heat capacity (a single value)                          *
 *                                                                      *
 * Last Modified: 18 Nov 2009 by D. Hasterok                            *
 * ******************************************************************** */
void computeheatcap(int aflag, double t, double Tref, double P, double T, double vf[],
  double *Cp) {

    int i;
    double CpPT[NMIN];
    double w[NMODE],vmode[NMODE];
    double R;

    modefrac(vf,w);

    /* ratio effective expansivity/expansivity */
    if (aflag) {
        R = adj_alpha(t,T);
    } else {
        R = 1;
    }

    /* compute mineral heat capacity */
    for (i = 0; i < NMIN; i++) {
        if (vf[i] != 0) {
            CpPT[i] = heatcap(i,Tref,P,T,R);
        } else {
            CpPT[i] = 0;
        }
    }

    /* compute rock heat capacity */
    modevalue(CpPT,vf,vmode);
    *Cp = wtmean(w,vmode);
}



/* ******************************************************************** *
 * COMPUTETHERMCOND - Compute rock thermal conductivity.                *
 *                                                                      *
 * Inputs:                                                              *
 *    P         pressure                                                *
 *    T         temperature                                             *
 *    vf[NMIN]  volume fraction of mineral assemblage                   *
 *   *lambda    thermal conductivity (a single value)                   *
 *    opt       conductivity option                                     *
 *                (0: lat + rad   1: lat only   2: rad only)            *
 *                                                                      *
 * Last Modified: 18 Nov 2009 by D. Hasterok                            *
 * ******************************************************************** */
void computethermcond(double P, double T, double vf[], double *lambda, int opt) {

    int i;
    double lambda_lat[NMODE],lambda_rad[NMODE];
    double lambda_eff[NMODE];
    double w[NMODE];

    modefrac(vf,w);

    /* compute mode conductivity */
    thermcond(P,T,vf,lambda_lat,lambda_rad);

    /* compute rock conductivity */

    /* Note: radiative conductivity                                     *
     *                                                                  *
     * A radiative component is only added to ol, and px components.    *
     *                                                                  *
     * The radiative contribution to thermal conductivity results from  *
     * blackbody emmissions and how optically transparent a sample is   *
     * to various wavelengths of light.  A few studies have been done   *
     * on olivine and few on other mantle phases (En Schatz & Simmons   *
     * [JGR, 1972]).  However the radiative conductivity for pyroxene   *
     * may be very similar to olivine [Hofmeister & Pertermann, Eur.    *
     * J. Mineral., 2008]                                               *
     *                                                                  *
     * Garnet and spinel may not contribute to radiative conductivity,  *
     * so the true radiative conductivity should probably be reduced    *
     * by a fraction proportional to the mode fraction [Tom Shankland,  *
     * pers. comm.].                                                    */

    for (i = 0; i < NMODE; i++) {
        lambda_eff[i] = lambda_lat[i] + lambda_rad[i];
        //printf("%lg ",lambda_rad[i]);
    }
    //printf("\n ");
    switch (opt) {
        case 0:
            *lambda = geomean(w,lambda_eff);
            break;
        case 1 :
            *lambda = geomean(w,lambda_lat);
            break;
        case 2 :
            /* use weighted mean for radiative conductivity since       *
             * very small values have unrealistic weight on the result  */
            *lambda = wtmean(w,lambda_rad);
            break;
        default :
            fprintf(stderr,"Warning (computethermcond): option not recognized, computing effective conductivity.\n");
            exit(0);
    }

    /* debugging... */
    //if (P == 0 && vf[17] + vf[18] + vf[19] + vf[20] != 0) {
    //    printf("%lg %lg %lg %lg %7.4lf  %lg %7.4lf\n",
    //        vf[17],vf[18],vf[19],vf[20],
    //        (lambda_lat[NMODE-1]+lambda_rad),w[7],*lambda);
    //}
} /* THERMCOND */



/* DENSITYPT - Computes P-T-dependent properties of several minerals.   *
 *                                                                      *
 *    int DENSITYPT(double P, double T,                                 *
 *    double *RHO, double *ALPHA, double *Cp) computes the density,     *
 *    RHO, thermal expansivity, ALPHA; and specific heat capacity, Cp   *
 *    of several mantle minerals with respect to pressure, P, and       *
 *    temperature, T.                                                   *
 *                                                                      *
 * See THERMEXP and COMPRESS for the approximations of the P, T         *
 * dependence and resulting integral solutions.                         *
 *                                                                      *
 * Last modified: 27-Feb 2009 by D. Hasterok                            */
void densityPT(int i, double Tref, double P, double T, double *alphaPT, double *rhoPT, double R) {

    double alphaT, ialphaT;
    double rhor, rhoP;

    /* reference density */
    if (fabs(Tref - 298) < 1e-4) {
        rhor = rho298[i];
    } else {
        rhor = ref_density(i,Tref);
    }
    // printf("298 K - %lg    %lg K - %lg\n",rho298[i],Tref,rhor);

    /* thermal expansivity */
    alphaT = R * thermexp(i,T);
    ialphaT = R * ithermexp(i,T,Tref);
    //printf("%lg %lg\n",alphaT,ialphaT);

    /* debugging */
    //for (i = 0; i < NMIN; i++) {
    //    fprintf(stdout,"mode %i:  %lg  %lg\n",i,alphaT[i],ialphaT[i]);
    //}
    //exit(0);

    /* Solve for P-dependent density */
    rhoP = densityP(i, P, rhor);
    //printf("%lg\n",rhoP);

    *alphaPT = alphaT * pow(rhoP/rhor, -deltaT[i]);
    //printf("%.3lg   %.2lg   %.3lg   %.3lg\n",P,T,alphaT,*alphaPT);

    /* Compute P-T-dependent density and expansivity */
    *rhoPT = rhoP*exp(-ialphaT);
    //printf("%lg\n",*rhoPT);
    // fprintf(stdout,"%lg %lg\n",alphaPT[i],rhoPT[i]);
} /* DENSITYPT */



/* DENSITYP - Computes the pressure dependent density.                  *
 *                                                                      *
 *    double DENSITYP(double P, double rho0, double KT, double Kp)      *
 *    computes the pressure dependent density using Newton's method.    *
 *                                                                      *
 * Last Modified: 23 Sept 2009 by D. Hasterok                           */
double densityP(int i, double P, double rhor) {

    double tol = 1e-5;          /* tolerance */
    double count = 1;           /* counter */
    double maxit = 100;         /* maximum number of iterations */
    double rho_o,rho_n;

    if (P == 0) {
        return(rhor);
    }

    rho_o = rhor*(1 + P/K[i][0]);    // initial guess (1/K[i][1] is compressibility)
    while (1) {
        rho_n = rho_o - (P - pressure(i,rhor,rho_o))/(-dpressure(i,rhor,rho_o));
        /* debugging... */
        //p = pressure(rho0,rho_old,K[i][1],K[i][2]);
        //dp = dpressure(rho0,rho_old,K[i][1],K[i][2]);
    
        if (fabs(rho_n - rho_o) < tol) {
            rho_o = rho_n;
            break;
        }
    
        // printf("%lg %lg, %lg %lg\n",P,rho_o,rhor,rho_n);
        if (count > maxit) {
            fprintf(stderr,"WARNING (density): Did not converge.\n");
            exit(0);
        }
    
        rho_o = rho_n;
        count++;
    }
    /* debugging... */
    //rhoP-rho0
    return(rho_n);
} /* DENSITYP */


/* Logarithmic Equation of State (EoS)                                  *
 *                                                                      *
 * PRESSURE - Logarithmic (EoS).                                        *
 *                                                                      *
 *    This logarithmic EoS formulation is from Poirier & Tarantola      *
 *    [PEPI, 1998] and should be superior to the Birch-Murnaghan        *
 *    equation of state at high pressures.                              *
 *                                                                      *
 * Last Modified: 22 Sept 2009 by D. Hasterok                           */
double pressure(int i, double rho0, double rhoP) {

    double lnrho;
    double p;

    lnrho = log(rhoP/rho0);

    p = K[i][0]*(rhoP/rho0)*lnrho*(1 + 0.5*(K[i][1] - 2)*lnrho);

    return(p);
} /* FUNC */



/* DPRESSURE - Derivative of Logarithmic EoS.                           *
 *                                                                      *
 * Last Modified: 27-Feb 2009 by D. Hasterok                            */
double dpressure(int i, double rho0, double rhoP) {

    double lnrho;
    double dp;

    //double zeta;
    //zeta = 0.5*(K[i][1] - 2);
    lnrho = log(rhoP/rho0);

    //dp = K[i][0]/rho0*((1 + lnrho)*(1 + zeta*lnrho) + zeta*lnrho);
    dp = K[i][0]/rho0*(1 + K[i][1]*lnrho*(1 + 0.5*lnrho));

    return(dp);
} /* DFUNC */



/* Thermal Expansion - The next three functions compute the temperature *
 * dependent thermal expansion, and its derivative and integral.        *
 *                                                                      *
 * THERMEXP - Temperature dependent thermal expansivity.                *
 *                                                                      *
 *    double *THERMEXP(double T) computes the thermal conductivity for  *
 *    the several minerals given the temperature, T in [K].             *
 *                                                                      *
 * ITHERMEXP - Integral of thermal expansivity w.r.t. temperature.      *
 *                                                                      *
 *    double *ITHERMEXP(double T) computes the integral of thermal      *
 *    expansivity as a function of temperature for the several          *
 *    minerals given the temperature, T in [K].                         *
 *                                                                      *
 * DTHERMEXP - Temperature derivative of thermal expansivity.           *
 *                                                                      *
 *    double *THERMEXP(double T) computes the derivative of thermal     *
 *    expansivity w.r.t. temperature for several minerals given the     *
 *    temperature, T in [K].                                            *
 *                                                                      *
 * Last modified: 22 Sept 2009 by D. Hasterok                           */
/* thermal expansivity */
double thermexp(int i, double T) {

    double alpha;
    
    alpha = a[i][0] + T*a[i][1] + a[i][2]/(T*T);

    return(alpha);
} /* THERMEXP */

/* integral of expansivity */
double ithermexp(int i, double T, double Tref) {
    
    //const double Tref = 293;
    double ialpha;
    
    ialpha = (T - Tref)*a[i][0] 
        + 0.5*(T*T - Tref*Tref)*a[i][1] 
        - (1/T - 1/Tref)*a[i][2];

    return(ialpha);
} /* ITHERMEXP */



/* derivative of expansivity */
double dthermexp(int i, double T) {

    double dalpha;
    
    dalpha = a[i][1] - 2/(T*T*T)*a[i][2];

    return(dalpha);
} /* DTHERMEXP */



/* HEATCAP - Computes heat capacity for mantle minerals.                *
 *                                                                      *     
 *    Cp = HEATCAP(P,T,C,ALPHA,DALPHAT,RHO) computes the heat capacity  *
 *    of minerals associated with the density, RHO; expansivity,        *
 *    ALPHA; temperature derivative of expansivity, DALPHAT; using the  *
 *    constants in C at the pressure, P, and temperature T.  Note that  *
 *    T must be in Kelvin and P in GPa.                                 *
 *                                                                      *
 * Last modified: 27-Feb 2009 by D. Hasterok                            */
double heatcap(int i, double Tref, double P, double T, double R) {

    double rhoPT;
    double alphaPT,dalphaT;
    double dCp_dP;
    double Cp;

    densityPT(i,Tref,P,T,&alphaPT,&rhoPT,R);
    dalphaT = dthermexp(i,T);
    /* compute the pressure influence on heat capacity see Osako et *
     * al. [PEPI, 2004], must multiply by 1e9 to convert            *
     * J/(kg-K-Pa) to J/(kg-K-GPa)                                  */
    dCp_dP = -T*(alphaPT*alphaPT + dalphaT)/rhoPT*1e9;
    
    /* heat capacity */
    Cp = c[i][0] + c[i][1]/sqrt(T) + (c[i][2] + c[i][3]/T)/(T*T) + dCp_dP*P;

    return(Cp);
} /* HEATCAP */



/* PHASEFRAC - Fraction of each mineral phases nomalized by mode.       */
void phasefrac(double vf[], double *X) {

    int i,j;
    double s = 0;

    for (i = 0; i < NMODE; i++) {
        s = 0;
        for (j = ind[i][0]; j <= ind[i][1]; j++) {
            s += vf[j];  /* sum of weighted modes */
        }

        for (j = ind[i][0]; j <= ind[i][1]; j++) {
            if (s == 0) { X[j] = 0; continue; };
            X[j] = vf[j]/s;  /* sum of weighted modes */
        }
    }
    /* debugging */
    //for (i = 0; i < NMIN; i++) {
    //    fprintf(stdout,"%lg  ",X[i]);
    //}
    //fprintf(stdout,"\n");
} /* PHASEFRAC */


/* BVALUE - Pressure term for computing thermal conducivity.            */
void bvalue(double vf[], double *b) {

    int i,j;
    double X[NMIN];

    phasefrac(vf,X);

    for (i = 0; i < NMODE; i++) {
        //printf("%lg \n",X[i]);
        b[i] = 0;
        for (j = ind[i][0]; j <= ind[i][1]; j++) {
            b[i] += X[j]*K[j][1]/K[j][0];  /* sum of weighted modes */
        }
    }
    /* debugging */
    //for (i = 0; i < NMODE; i++) {
    //    fprintf(stdout,"%6.4lf ",b[i]);
    //}
    //fprintf(stdout,"\n");
} /* BVALUE */



/* CONDCONST */
void condconst(double vf[], double T, double *lambda0) {

    double X;

    /* compute mole fraction of each mode end-member. This is   *
     * different from the actual phase fraction because of the  *
     * sparce data on end-member thermal conductivities.  As    *
     * data are added this section will likely need updating.   */

    /* alpha-quartz     aqz                                     */
    if (vf[0] != 0) { lambda0[0] = kl[0][0]; } else { lambda0[0] = 0; }
    /* beta-quartz      bqz                                     */
    if (vf[1] != 0) { lambda0[1] = kl[1][0]; } else { lambda0[1] = 0; }

    /* orthoclase       or                                      */
    if (vf[2] != 0) {
        lambda0[2] = kl[2][0];
    } else { lambda0[2] = 0; }

    /* plagioclase      ab/(an + ab)                            */
    if (vf[4] + vf[5] != 0) {
        X = vf[4]/(vf[4] + vf[5]);
        lambda0[3] = kl[3][0] + X*(kl[3][1] + X*kl[3][2]);
    } else { lambda0[3] = 0; }

    /* muscovite        musc                                    */
    if (vf[6] + vf[7] != 0) { lambda0[4] = kl[4][0]; } else { lambda0[4] = 0; }

    /* biotite          bt                                      */
    if (vf[8] + vf[9] != 0) { lambda0[5] = kl[5][0]; } else { lambda0[5] = 0; }

    /* epidote          ep                                      */
    if (vf[10] != 0) { lambda0[6] = kl[6][0]; } else { lambda0[6] = 0; }

    /* chlorite         clin/(clin + daph)                      */
    if (vf[11] + vf[12] != 0) { lambda0[7] = kl[7][0]; } else { lambda0[7] = 0; }

    /* amphibolite      hbl/(hbl + trem + fact + parg)          */
    if (vf[13] + vf[14] + vf[15] + vf[16] != 0) {
        lambda0[8] = (vf[13]*kl[8][0] + vf[14]*kl[8][1] 
            + vf[15]*kl[8][3] + vf[16]*kl[8][4])
            / (vf[13] + vf[14] + vf[15] + vf[16]);
    } else { lambda0[8] = 0; }

    /* clinopyroxene    di/(di + hd + jd)                       */
    if (vf[17] + vf[18] + vf[19] != 0) { lambda0[9] = kl[9][0]; } else { lambda0[9] = 0; }

    /* wollastonite     wo                                      */
    if (vf[20] != 0) { lambda0[10] = kl[10][0]; } else { lambda0[10] = 0; }

    /* orthopyroxene    en/(en + fs)                            */
    if (vf[21] + vf[22] != 0) { lambda0[11] = kl[11][0]; } else { lambda0[11] = 0; }

    /* olivine          fo/(fo + fa)                            */
    if (vf[23] + vf[24] != 0) {
        X = vf[23]/(vf[23] + vf[24]);
        lambda0[12] = kl[12][0] + X*(kl[12][1] + X*kl[12][2]);
    } else { lambda0[12] = 0; }

    /* spinel           sp/(sp + mt)                            */
    if (vf[25] + vf[26] + vf[27] != 0) { 
        lambda0[13] = (vf[25]*kl[13][0] + vf[26]*kl[13][1] + vf[27]*kl[13][2])
            / (vf[25] + vf[26] + vf[27]);
    } else { lambda0[13] = 0; }

    /* garnet           (py + gr)/(py + al + gr + spess)        *
     *                                                          *
     * Note: 13-Oct 2009                                        *
     *                                                          *
     * While there are few measurements to make an accurate     *
     * assessment, gr is assumed ~= py and spess ~= al because  *
     * of the relative similarity in substituted cation size.   */
    if (vf[28] + vf[29] + vf[30] + vf[31] != 0) { 
        X = (vf[28] + vf[30])/(vf[28] + vf[29] + vf[30] + vf[31]);
        lambda0[14] =  kl[14][0] + X*(kl[14][1] + X*kl[14][2]);
    } else { lambda0[14] = 0; }

    /* carbonate        cc/(cc + arag + do)                     */
    if (vf[32] + vf[33] + vf[34] != 0) { lambda0[15] = kl[15][0]; } else { lambda0[15] = 0; }

    /* aluminosilicate                                          */
    if (vf[35] + vf[36] + vf[37] != 0) {
        lambda0[16] = (vf[35]*kl[16][0] + vf[36]*kl[16][1] + vf[37]*kl[16][2])
            / (vf[35] + vf[36] + vf[37]);
    } else { lambda0[16] = 0; }

    /* rutile           rt                                      */
    if (vf[38] != 0) { lambda0[17] = kl[17][0]; } else { lambda0[17] = 0; }

    /* sanidine */
    if (vf[3] != 0) {
        lambda0[18] = kl[18][0];
    } else { lambda0[18] = 0; }



    /* debugging */
    //for (i = 0; i < NMODE; i++) {
    //    fprintf(stdout,"%lg  ",lambda0[i]);
    //}
    //fprintf(stdout,"\n");
}




/* THERMCOND - Thermal conductivity of mantle phases.                   *
 *                                                                      *
 * Inputs:                                                              *
 *     type     -1: derivative, 0: conductivity, 1: integral            *
 *     P        pressure                                                *
 *     T        temperature                                             *
 *     vf[]     volume fractions                                        *
 *                                                                      *
 * Outputs:                                                             *
 *     val      returned value, type                                    *
 *                                                                      *
 * [KLAT,KRAD] = THERMCOND(P,T) computes the thermal conductivity of    *
 * mantle phases for ol, opx, cpx, sp, and gt given the pressure, P,    *
 * and temperature, T.  Temperature must be in Kelvin and pressure in   *
 * GPa.  The effective thermal conductivity is the sum of the lattace   *
 * conductivity, KLAT, and radiative conductivity, KRAD.                *
 *                                                                      *
 * KRAD is only computed for olivine and should be similar for          *
 * pyroxenes and much smaller for garnet.                               *
 *                                                                      *
 * Last modified: 19-Feb 2009 by D. Hasterok                            */
void thermcond(double P, double T, double vf[], double *lambda_lat, double *lambda_rad) {

    int i;
    double b[NMODE],lambda0[NMODE];

    bvalue(vf,b);
    condconst(vf,T,lambda0);
    //for (i = 0; i < NMODE; i++) {
    //    fprintf(stdout,"k298: %6.4lf   b: %6.4lf\n",lambda0[i],b[i]);
    //}

    /* b-values: Fo90     ~0.032 from Xu et al. [PEPI, 2001]            *
     *           Al75Py25 ~0.045 from Osako et al. [PEPI, 2004]         *
     *                                                                  *
     * Hofmeister & Pertermann [Eur. J. Mineral., 2008] suggests that   *
     * the type of experiment performed by Osako et al. may             *
     * artificially inflate the pressure dependence.                    *
     *                                                                  *
     * Lattace thermal conductivity                                     *
     * ----------------------------                                     *
     * Note that Hofmeister [Science, 1999] mathmatically represents    *
     * the thermal conductivity with an additional term,                *
     * exp(1/3 + 4*gammaT int(alpha dT)).  However, she also assumes    *
     * that the value of a is 0.33 for all silicates.  Since I am       *
     * using a value for a that fits the temperature dependent          *
     * conductivity data, the extra term is unnecessary (? I think).    *
     *                                                                  *
     * I am using the pressure term from Hofmeister [Science, 1999]     *
     * and the temperature dependence model from Xu et al. [PEPI,       *
     * 2004].  Xu suggests that the exponent is likely 0.5 for olivine  *
     * and may be for other minerals as well.  Although my fits show    *
     * that they can range between 0.56 and 0.3 for mantle minerals     *
     * and is actually -0.21 for plagioclase.                           *
     *                                                                  *
     * Since solid solutions of mineral do not mix linearly, it is      *
     * more useful to compute the mode immediately rather than compute  *
     * each mineral independent of each other.                          *
     *                                                                  *
     * The exact formula that I am using is:                            *
     *                                                                  *
     *                                (298)^n       K'                  *
     *    klat = (k0 + k1*X + k2*X^2) (---)  (1 + P*--)                 *
     *                                ( T )         KT                  *
     *                                                                  *
     * where klat is the lattice contribution to thermal conductivity,  *
     * k0-2 are constants used to determine the thermal conductivity    *
     * of a mixture at 298 K, k298; X is the mole fraction of the mode  *
     * end-member; P is pressure, and KT and K' are the isothermal      *
     * bulk modulus and its first pressure derivative.  If the mixing   *
     * constants are unknown, then a linear combination is assumed.     */

    /* If pressure dependence is to be ignored, uncomment line below */
    //b = zeros(size(k_298));
    for (i = 0; i < NMODE; i++) {
        if (lambda0[i] != 0) {
            lambda_lat[i] = lambda0[i]*pow(298/T,n[i])*(1 + b[i]*P);
            lambda_rad[i] = kr[i][0]*(1 + erf(kr[i][1]*(T - kr[i][2])));
            //lambda_rad[i] = kr[i][0] + T*(kr[i][1] + T*(kr[i][2] - T*kr[i][3]));
        } else {
            lambda_lat[i] = 0;
            lambda_rad[i] = 0;
        }
        //printf("%i  %lg  %lg\n",i,lambda_lat[i],lambda_rad[i]);
    }

    /* debugging... */
    //for (i = 0; i < NMODE; i++) {
    //    fprintf(stdout,"%6.4lf   ",lambda_lat[i]);
    //}
    //fprintf(stdout,"%6.4lf\n",*lambda_rad);
}



/* ITHERMCOND - Integral of thermal conductivity.                       */
//void ithermcond(double P, double T, double vf[], double *G_lat, double *G_rad) {
//
//    int i;
//    int Tref = 273;                  /* [K] reference temperature */
//    double b[NMODE],lambda0[NMODE];
//
//    double a = 0, g = 0;
//    double w[NMODE];
//    double Glat,Grad;
//
//    modefrac(vf,w);
//
//    bvalue(vf,b);
//    condconst(vf,T,lambda0);
//
//    // Conductivity integral
//    for (i = 0; i < NMODE; i++) {
//        if (lambda0[i] != 0) {
//            a = a - n[i]*w[i];
//            g = g*pow(lambda0[i]*pow(298,n[i])*(1 + b[i]*P),w[i]);
//    
//            /* Integral of klat */
//            if (a != -1) {
//                G_lat[i] = g*(pow(T,a + 1) - pow(Tref,a + 1))/(a + 1);
//            } else {
//                G_lat[i] = g*(log(T) - log(Tref));
//            }
//        
//            /* Integral of krad */
//            G_rad[i] = T*(kr[i][0] + T*(0.5*kr[i][1] 
//                + T*(kr[i][2]/3 - 0.25*T*kr[i][3])))
//                - Tref*(kr[i][0] + Tref*(0.5*kr[i][1] 
//                + Tref*(kr[i][2]/3 - 0.25*Tref*kr[i][3])));
//        } else {
//        }
//    }
//} /* ITHERMCOND */



/* Dervative of thermal conductivity.                                   */
void dthermcond(double P, double T, double vf[], double *dlambda) {

    int i;
    double b[NMODE],lambda0[NMODE];
    double S,W,w[NMODE];
    double u;

    double lambda_eff[NMODE],lambda_lat[NMODE],lambda_rad[NMODE];
    double dlambda_lat[NMODE],dlambda_rad[NMODE];
    double lambda;

    /* mode weights and sum of weights */
    modefrac(vf,w);
    for (i = 0; i < NMODE; i++) {
        W = W + w[i];
    //    printf("%2i  %lg\n",i,w[i]);
    }

    //for (i = 0; i < NMIN; i++) {
    //    printf("%2i  %lg\n",i,vf[i]);
    //}


    /* compute derivatives first */
    bvalue(vf,b);
    condconst(vf,T,lambda0);
    for (i = 0; i < NMODE; i++) {
        if (lambda0[i] != 0) {
            /* derivative of lattice conductivity */
            dlambda_lat[i] = n[i]*lambda0[i]*pow(298,n[i])*pow(T,-(n[i]+1))*(1 + b[i]*P);
        
            /* derivative of radiative conductivity */
            u = kr[i][1]*(T - kr[i][2]);
            dlambda_rad[i] = kr[i][0]*kr[i][1]*exp(-u*u)*ISQRTPI;
            //dlambda_rad[i] = kr[i][1] + T*(2*kr[i][2] - T*3*kr[i][3]);
        } else {
            dlambda_lat[i] = 0;
            dlambda_rad[i] = 0;
        }
    }

    /* compute thermal conductivity */
    thermcond(P,T, vf, lambda_lat,lambda_rad);

    S = 0;
    for (i = 0; i < NMODE; i++) {
        //printf("%2i  %lg  %lg  %lg %lg\n",i,lambda_lat[i],lambda_rad[i],dlambda_lat[i],dlambda_rad[i]);
        /* effective conductivity of mode */
        lambda_eff[i] = lambda_lat[i] + lambda_rad[i];

        if (lambda0[i] == 0) { continue; }

        /* sum of derivative from each mode */
        S = S + w[i]
            * (dlambda_lat[i] + dlambda_rad[i])
            / lambda_eff[i];
    }

    /* effective conductivity */
    lambda = geomean(w,lambda_eff);

    /* derivative of effective conductivity */
    *dlambda = lambda*S/W;

} /* DTHERMCOND */


double adiabat(double Tp, double alphaPT, double Cp) {

    double Gamma;

    Gamma = gravacc * alphaPT * Tp / Cp;

    return(Gamma);
}



/* SP2GT - Spinel to Garnet transition.                                 *
 *                                                                      *
 *    double SP2GT(double T) computes the pressure, P, of the spinel    *
 *    to garnet transition given the temperature, T in Kelvins.         *
 *                                                                      *
 *    This boundary is estimated from an inversion of the data in       *
 *    Robinson and Wood [EPSL, 1998], Walter et al. [GCA, 2002], and    *
 *    Klemme and O'Neill [Contrib. Mineral. Pet., 2000a].               *
 *                                                                      *
 * Last modified: 22-Sept 2009 by D. Hasterok                           */
double sp2gt(double T) {

    double P;

    //P = 1.3783 + exp(3.6825e-3*T - 6.3009); /* unweighted inversion result */
    P = 1.4209 + exp(3.9073e-3*T - 6.8041);   /* weighted inversion result   */

    return(P);
} /* SP2GT */



/* MODEVALUE - Computes mode average.                                   *
 *                                                                      *
 *    double *MODEVALUE(double VMIN[NMIN],double MR[NMIN]) computes     *
 *    the average mode value, VMODE[NMODE], for a given property for    *
 *    minerals associated with the properties, VMIN, for the phases     *
 *    with volume fractions, VF.                                        *
 *                                                                      *
 *    If a particular mineral is not needed the column must still be    *
 *    included, but can be any value.                                   *
 *                                                                      *
 *    Note chromite, rutile and quartz and micas should be added for    *
 *    a more complete modeling of igneous rocks.                        *
 *                                                                      *
 * Last modified: 22-Sept 2009 by D. Hasterok                           */
void modevalue(double vmin[], double vf[], double *vmode) {

    int i,j;
    double s;

    for (i = 0; i < NMODE; i++) {
        vmode[i] = 0;
        s = 0;
        for (j = ind[i][0]; j <= ind[i][1]; j++) {
            if(vf[j]==0) continue;
            s += vf[j];                 /* sum of weights */
            vmode[i] += vmin[j]*vf[j];  /* sum of weighted modes */
        }

        /* if all mole fractions are 0, then set mode value to 0 */
        if (s == 0) {
            vmode[i] = 0;
            continue;
        }

        /* divide by sum of weights */
        vmode[i] = vmode[i]/s;
    }

    /* debugging */
    //fprintf(stdout,"mode value:\n");
    //for (i = 0; i < NMODE; i++) {
    //    fprintf(stdout,"mode %i: %lg\n",i,vmode[i]);
    //}
} /* MODEVALUE */



/* MODEFRAC - Computes the mode fractions.                              */
void modefrac(double vf[], double *w) {

    int i,j;
    double s = 0;

    for (i = 0; i < NMIN; i++) {
        s += vf[i];
    }

    for (i = 0; i < NMODE; i++) {
        w[i] = 0;
        for (j = ind[i][0]; j <= ind[i][1]; j++) {
            w[i] += vf[j];  /* sum of weighted modes */
        }
        w[i] = w[i]/s;
    }

    /* debugging */
    //fprintf(stdout,"mode fraction:\n");
    //for (i = 0; i < NMODE; i++) {
    //    fprintf(stdout,"mode %i: %lg\n",i,w[i]);
    //}
} /* MODEFRAC */



/* CHECKWTS - Checks sum of mode and mineral fractions.                 */
int checkwts(double vf[]) {

    int i;
    static double tol = 1e-4;
    double s = 0;

    /* sum mineral fractions */
    /* make sure mineral fractions sum to 1 */
    for (i = 0; i < NMIN; i++) {
        s += vf[i];
    }

    if (fabs(s - 1) < tol) {
        return(1);
    }

    fprintf(stderr,"Warning (checkwts): Mineral volume fractions sum to %lg.\n",s);
    return(-1);
}


    
/* WTMEAN - Computes weighted mean.                                     *
 *                                                                      *
 *    double WTMEAN(double W, double VMODE) computes the weighted mean  *
 *    of the mode properties, VMODE, with weights, W.                   *
 *                                                                      *
 * Last Modified: 22 Sept 2009 by D. Hasterok                           */
double wtmean(double w[], double vmode[]) {

    double m = 0;
    int i;

    for (i = 0; i < NMODE; i++) {
        m += w[i]*vmode[i];
    }

    return(m);
} /* WTMEAN */



/* GEOMEAN - Geometric mean                                             *
 *                                                                      *
 *    double GEOMEAN(double W, double VMODE) computes the geometric     *
 *    mean of the mode properties, VMODE, with weights, W.              *
 *                                                                      *
 * Last Modified: 22 Sept 2009 by D. Hasterok                           */
double geomean(double w[], double vmode[]) {

    int i;
    double S = 0;
    double W = 0;

    for (i = 0; i < NMODE; i++) {
        W += w[i];
        if (w[i] > 0) {
            S += log(vmode[i])*w[i];
        }
    }

    return(exp(S/W));
} /* GEOMEAN */



/* VF2MF - Convert from volume percent to weight percent.               */
void vf2mf(double vf[], double mf[], double rho[]) {

    int i;
    double s = 0;

    // Compute weight of each component
    for (i = 0; i < NMIN; i++) {
        if (vf[i] == 0) {
            mf[i] = 0;
            continue;
        }
        mf[i] = vf[i] * rho[i];
        s += mf[i];
    }

    // Renormalize
    for (i = 0; i < NMIN; i++) {
        mf[i] = mf[i]/s;
    }

    //for (i = 0; i < NMIN; i++) {
    //    fprintf(stdout," %8.4lf\n",mf[i]*100);
    //}
    //fprintf(stdout,"\n");
} /* VF2MF */


/* MF2VF - Convert from weight percent to volume percent.               */
void mf2vf(double mf[], double vf[], double rho[]) {

    int i;
    double s = 0;

    // Compute mass of each component
    for (i = 0; i < NMIN; i++) {
        if (mf[i] == 0) {
            vf[i] = 0;
            continue;
        }
        vf[i] = mf[i]/rho[i];
        s += vf[i];
    }

    // Renormalize
    for (i = 0; i < NMIN; i++) {
        vf[i] = vf[i]/s;
    }

    //for (i = 0; i < NMIN; i++) {
    //    fprintf(stdout," %8.4lf\n",vf[i]*100);
    //}
    //fprintf(stdout,"\n");
} /* MF2VF */
