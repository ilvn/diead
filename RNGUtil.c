/*
*    Project:    DIEAD
*    File:       RNGUtil.c
*    Author:     Ilya Levin
*
*    A small suite to test 32-bit pseudorandom generators.
*
*    Copyright 2003,2016 Literatecode
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
*    The code in this file is modified from the code attributed to the Center
*    for Information Security & Cryptography (CISC), Department of Computer
*    Science, HKU.
*
*    The code is mostly modified to be properly compilable on a 64-bit platform
*    with a recent clang or gcc compiler.
*
*/
/*
 *   Filename: RNGUtil.c
 */
#include "RNGUtil.h"



/* Compute gamma(z) for the cases where 2z is an integer */
static double G(double z)
{
    int tmp = (int) (2. * z);

    if ( (z < 1) || (tmp != 2 * lround(z)) )
        printf("Error in calling G(z)!!!");

    switch(tmp) {
        case 1:
            return sqrt(PI);
        case 2:
            return 1.;
        default:
            return (z - 1.) * G(z - 1.);
    }
} /* G */


/* Compare 2 nos, returns -1, 0, 1 when <, ==, >, respectively */
int ucmpr( const void *i1, const void *i2)
{
    const uniform *u1 = (const uniform *)i1, *u2 = (const uniform *)i2;

    if( *u1 < *u2 ) return -1;
    if( *u1 == *u2 ) return 0;

    return 1;

} /* ucmpr */

/*read in a uniform random number from a file*/
/* Input DIM no. of unsigned integer from the file.

  */
uniform uni(char *filename)
{
    static FILE *infile;

    static char isopen = 'n';
    static counter count = DIM, m = 0;
    static uniform uniran[DIM];

    /* Close the file when filename = "close" */
    if( strcmp(filename, "close") == 0 ) {
        fclose(infile);
        isopen = 'n';
        count = DIM;

        return 0;
    }

    /* if (++count) == 0..DIM-1, return one no. in the buffer */
    if( (++count) < m ) {
        return uniran[count];
    }

    /* If it is the very first time, */
    if( isopen == 'n' ) {
        if( (infile = fopen(filename, "rb")) == NULL ) {
            printf("can't open file %s!!!\n", filename);
            exit(1);
        }

        isopen = 'y';
    }

    /* The buffer is empty, input DIM nos */
    m = fread( uniran, sizeof(uniform), DIM , infile );
    if (m == 0) {
        fclose(infile);
        if( (infile = fopen(filename, "rb")) == NULL ) {
            printf("can't open file %s!!!\n", filename);
            exit(1);
        }
        printf("\n### %s is reopened!\n", filename);
    }

    count = 0;

    return uniran[count];

} /* uni */


double erfc( double);

/* Error function: erf(x) = (2/sqrt(pi)) * int(0, x, e^(-t^2) dt )
   More accurate than erfc(x) when x <= 0.46875                    */
double erf(double x)
{

    const double p[] = {2.4266795523053175e+2,
                        2.1979261618294152e+1,
                        6.9963834996191355,
                        -3.5609843701815385E-2
                       };
    const double q[] = {2.1505887586986120e+2,
                        9.1164905404514901e+1,
                        1.5082797630407787e+1,
                        1.0000000000000000
                       };

    double nsum = 0., dsum = 0., xx = x * x, tmp = 1.;
    int i;

    if(x == 0.) return 0.;
    if(x >= .46875) return 1. - erfc(x);

    for(i = 0; i < 4; ++i) {
        nsum += p[i] * tmp;
        dsum += q[i] * tmp;
        tmp *= xx;
    }

    return x * nsum / dsum;

} /* erf */


/* a complement of erf(x). erfc(x) = 1. - erf(x)
   More accurate that erf(x) when x >= 0.46875    */
double erfc(double x)
{
    const double p2[] = {3.004592610201616005e+2,
                         4.519189537118729422e+2,
                         3.393208167343436870e+2,
                         1.529892850469404039e+2,
                         4.316222722205673530e+1,
                         7.211758250883093659,
                         5.641955174789739711e-1,
                         -1.368648573827167067e-7
                        };

    const double q2[] = {3.004592609569832933e+2,
                         7.909509253278980272e+2,
                         9.313540948506096211e+2,
                         6.389802644656311665e+2,
                         2.775854447439876434e+2,
                         7.700015293522947295e+1,
                         1.278272731962942351e+1,
                         1.000000000000000000
                        };

    const double p3[] = { -2.99610707703542174e-3,
                          -4.94730910623250734e-2,
                          -2.26956593539686930e-1,
                          -2.78661308609647788e-1,
                          -2.23192459734184686e-2
                        };
    const double q3[] = {1.06209230528467918e-2,
                         1.91308926107829841e-1,
                         1.05167510706793207,
                         1.98733201817135256,
                         1.00000000000000000
                        };

    const double *p = p2, *q = q2;
    size_t n = 8;
    double nsum = 0., dsum = 0., xx = x, tmp = 1.;
    unsigned int i;

    if(x < .46875) return 1. - erf(x);

    if(x >= 4) {
        p = p3;
        q = q3;
        n = 5;
        xx = x * x * x * x;
    }

    for(i = 0; i < n; ++i) {
        nsum += p[i] * tmp;
        dsum += q[i] * tmp;
        tmp *= xx;
    }

    tmp = nsum / dsum;
    if(x >= 4.) {
        tmp = tmp / (x * x) + 1 / sqrt(PI);
        tmp /= x;
    }

    return exp(-x * x) * tmp;

} /* erfc */

/* c. d. f. of standard normal */
double Phi(double x)
{
    if(x > 0.)
        return ( 1. + erf(x / sqrt(2.)) ) / 2.;
    else
        return erfc(-x / sqrt(2.)) / 2.;

} /* Phi */


/* p.d.f of Chi-square */
static double chisq(int df, double x)
{
    return ( pow(x / 2., (df - 2.) / 2.) * exp(-x / 2.) / (2.*G(df / 2.)) );

} /* chisq */


/* c.d.f of Chi-square */
double Chisq(int df, double x)
{
    switch(df) {
        case 1:
            return 2.*Phi(sqrt(x)) - 1.;
        case 2:
            return 1. - exp(-x / 2.);
        default:
            return ( Chisq(df - 2, x) - 2.*chisq(df, x) );
    }
} /* Chisq */



double Chisq2( double d, double z)
{
    double v, t;

    if (z <= 0.) return(0.);

    t = (exp(log(z / d) / 3.) - 1. + 2 / (9 * d)) / sqrt(2 / (9 * d));

    if (t > 5.) return (1.);

    if (t + 5. < 0.) return(0.);

    v = (t > 0.) ? 1. / (1. + .33267 * t) : 1. / (1. - .33267 * t);

    v = exp(-.5 * t * t) * (.1740120799 + (.3739278012 * v - .04793993633) * v) * v;

    if (t > 0) v = 1. - v;

    return(v);

}   /* Chisq */



/* p.d.f of Poisson distribution */
double Poisson(double lambda, int k)
{
    if(k == 0) return exp(-lambda);

    return exp(-lambda) * pow(lambda, k) / G(k + 1);

} /* Poisson */

/***************************************************************
  The following implements a modified Kolmogorov-Smirnov goodness-of-fit test.
  The test-statistic is (FN(X)-X)**2/(X*(1-X)) (Anderson-Darling)
  where X is a uniform under null hypothesis. FN(X) is the empirical
  distribution of X.*/


/* prototypes: */
double ad(double);           /* asymptotic AD distribution */
double adfix(double, int);   /* adjusts for particular n */
int comp(const void*, const void*); /* Compare function for qsort */

/*----------------------------------------*/
int comp(const void *i, const void *j)
{
    /* return *(float *)i - *(float *)j; }*/
    return *(const int *)i - *(const int *)j;
}
/*----------------------------------------*/
double  adfix(double p, int n)
{
    double c;
    c = 1. / (57. + .12 * n - 235. / n);
    if(p < c)return .7 * (p - sqrt(c * p)) / n;
    if(p < .8) return (p - c) * (.8 - p) * (.384 - .272 * p) / n;
    return (p - .8) * (1. - p) * (4. - 5.84 * p) / n;
}

/*-----------------------------------------*/
double ad(double z)
{
    if(z < .06) return 0.;
    if(z < 2.) return 2.*exp(-1.2337 / z) * (1. + z / 8. - .04985 * z * z / (1.325 + z)) / sqrt(z);
    if(z < 4.) return 1. - .6621361 * exp(-1.091638 * z) - .95059 * exp(-2.005138 * z);
    return 1. - .49385879 * exp(-1.050321 * z) - .5946211 * exp(-1.527198 * z);
}
/*------------------------------------------*/
float kstest(float x[], int n)
{
    double z, t, p;
    int i;

    if (n == 1) return x[0];

    qsort(x, (uint32_t)n, sizeof(float), comp);
    z = (double) - n * (n + 0.);
    for(i = 0; i < n; i++)   {
        t = (double)x[i] * (1. - (double)x[n - 1 - i]);
        if(t < 1.e-30) t = 1.e-30;
        z = z - (i + i + 1) * log(t);
    }
    p = ad(z / n);
    return (float) (p + adfix(p, n));
}
