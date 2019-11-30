/*
*    Project:    DIEAD
*    File:       RNGTest.c
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
 *   Filename: RNGTest.c
 */
#include <math.h>
#include "main.h"


/* -------------------------------------------------------------------------- */
/* -------Library of RNG tests -----------------------------------------------*/
/* ---------------------------------------------------------------------------*/

/* Gorilla test                                                               */
double gorilla( uint32_t (*rng)(void))
{
    const uint32_t a = 1, m = 26;
    const double n1 = 26., stdd = 4170;

    uint32_t k2 = (uint32_t)lround(pow(2., a * m - 5));    /* no. of entries in t[]     */
    uint32_t k3 = (uint32_t)lround(pow(2., a * m) - 1);    /* mask for getting one word */

    uint32_t bit[32];
    uint32_t *t;      /* table to keep a state of each urn: 0=empty, 1=filled */
    uint32_t i, ii, j, k, n, u, nclsn, nholes, x, y;

    double mean, p;
    float xa[33];

    printf("\n### Gorilla test starts...");

    n = (uint32_t)lround(pow(2., n1)); /* no. of words */  /* n1*/

    mean = pow(2., a * m) * exp(-pow(2., n1 - a * m));

    /* initialize bit[i] to 2^i */
    for (bit[0] = 1, i = 1; i < 32; ++i) {
        bit[i] = bit[i - 1] * 2;
    }

    /* allocate a table for keeping the status of cells: 0=empty, 1=occupied  */
    t = (uint32_t *) malloc( k2 * sizeof(uint32_t) );
    if (t == NULL) {
        printf("ERR: Out of memory\n");
        return 0.;
    }

    for (k = 32; k >= 1; --k) {
        j = k;

        /* reset the table of status */
        for (i = 0; i < k2; ++i) {
            t[i] = 0;
        }

        j = j - 1; /* j becomes the no. of bits needed to shift to the right */

        u = 0;
        for (ii = 1; ii <= a * (m - 1); ++ii)
            u = (u << 1) + ((rng() >> j) & 1);


        /* generate n words */
        nclsn = 0;
        for (i = 1; i <= n; ++i) {
            for (ii = 1; ii <= a; ++ii)
                u = (u << 1) + ((rng() >> j) & 1);
            u = u & k3;

            x = u >> 5;
            y = u & 31;
            if ( (t[x] & bit[y]) != 0) {
                ++nclsn;
            } else {
                t[x] = t[x] | bit[y];
            }
        }

        nholes = k2 * 32 - (n - nclsn);

        p = Phi( (nholes - mean) / stdd );

        xa[k - 1] = (float) p;

        if (k % 8 == 0) printf("\np:<b%2d>...<b%2d>|", k, k - 7);
        if (k % 8 != 1) printf(" %6.4f ", p);
        else printf("%6.4f", p);
    }

    p = (double)kstest(xa, 32);

    printf("\n### Gorilla test is complete. Overall p = %6.4f\n",  p);

    free(t);

    return p;

} /* gorilla() */



/***************************************************************************************/
/******** Collision Test ***************************************************************/

/* pcoll( m, n, c) computes the probability that <= c collisions occurs
   when n balls are thrown into n urns. It uses Knuth's algorithm when n or m
   is <= 2^16 (See D. Knuth, The art of computer programming, 3rd ed. page 71,
   1997. ). It uses an approximation formula otherwise (See Tsang et al.,
   Tuning the collision test for stringency). The largest absolute error in the
   values returned that are less than 0.05 or over 0.95 (regions possibly leading
   to rejections of hypotheses) is 0.000446. This upper limit of error occurs
   when m = n = 2^17 and c = 48404.
*/
double pcoll(uint32_t m, uint32_t n, uint32_t c)
{
    const int MN = 65536;
    const double epsilon = 1e-20;

    double *A, mm, q, r, mean, var, std, cdf;

    uint32_t i, j, j0, j1;

    mm = m;

    if ((m <= MN) || (n <= MN) ) {  /* Knuth's method */
        A = (double *) malloc( (n + 1) * sizeof(double));
        if (A == NULL) {
            printf("ERR: Out of memory at pcoll()\n\n");
            return -1.;
        }

        /* S1 */
        for (j = 0; j <= n; ++j) {
            A[j] = 0.;
        }
        A[1] = 1.;

        j0 = 1;
        j1 = 1;

        /* S2 */
        for (i = 1; i < n; ++i) {
            j1 = j1 + 1;

            for (j = j1; j >= j0; --j) {
                A[j] = (j / mm) * A[j] + ((1. + (1. / mm)) - (j / mm)) * A[j - 1];
            }

            if (A[j0] < epsilon) A[j0++] = 0.;
            if (A[j1] < epsilon) A[j1--] = 0.;
        }

        /* Compute the cdf */
        if (n - c > j1) {
            free(A);
            return 0.;
        }
        if (n - c < j0) {
            free(A);
            return 1.;
        }

        cdf = A[j1];
        while (n - c < j1)
            cdf = cdf + A[--j1];

        free(A);
        return cdf;
    }

    else { /* W. W. Tsang's method */
        q = exp(n * log(1. - 1. / mm));
        r = exp(n * log(1. - 2. / mm));
        mean = mm * q - mm + n;
        var = mm * (q + mm * r - r - mm * q * q);
        std = sqrt( var);

        return Phi( (c - mean) / std );
    }

}       /* pcoll() */



/*
   collision() conducts the collision test on the 0th to 31th bit sequences of
   the words generated by the 32-bit generator rng. j=0 indicates the sequence
   formed from the least significant bits of the words and j=31 indicates the
   most significant bit sequence. The test throws 2^20 balls into 2^20 urns.
   The destination of each ball is determined by log_2 m bits extracted from
   the bit sequence. The bits are obtained from the jth bit of the 32-bit word
   returned by calling rng(). (See Tsang's paper "Tuning the collision test
   for stringency" for details if needed.)
*/
double collision( uint32_t (*rng)(void) )
{
    const uint32_t logm = 24;

    uint32_t bit[32];
    uint32_t *t;        /* the table used for keeping track of occupancy */
    uint32_t m2;        /* the number of 32-bit words needed in t[] */
    uint32_t i, ii, k, j, j2, m1, m, n, nclsn, u, x, y;

    float xa[33];

    double p;


    printf("\n### Collision test (No. of holes = 2^%u; No. of balls = 2^%u) starts...",  logm, logm);

    m = n = 1;
    for (i = 1; i <= logm; ++i)
        m = n = 2 * n;


    /* initialize bit[i] to 2^i */
    for (bit[0] = 1, i = 1; i < 32; ++i) {
        bit[i] = bit[i - 1] * 2;
    }

    m1 = (uint32_t) (log(m) / log(2.) + 0.5);

    /* allocate a table for keeping the status of cells: 0 is empty and 1 is occupied */
    m2 = m / 32;
    t = (uint32_t *) malloc(  m2 * sizeof(uint32_t));
    if (t == NULL) {
        printf("ERR: Out of memory\n");
        return 0.;
    }

    for (k = 32; k >= 1; --k) {
        j = k - 1;

        /* the mask for extracting the jth least significant bit in a word */
        j2 = (uint32_t)lround(pow(2., j));


        /* Set the table entries to zero. */
        for (i = 0; i < m2; ++i) {
            t[i] = 0;
        }

        /* throw n^2 balls */
        nclsn = 0;
        for (i = 1; i <= n; ++i) {
            u = ((rng() & j2) >> j);

            for (ii = 1; ii < m1; ++ii) {
                u = ( (u * 2) + ((rng() & j2) >> j) );
            }

            x = u >> 5;
            y = u & 31;
            if ( (t[x] & bit[y]) != 0) {
                ++nclsn;
            } else {
                t[x] = t[x] | bit[y];
            }
        }

        p = pcoll(m, n, nclsn);

        xa[k - 1] = (float) p;

        if (k % 8 == 0) printf("\np:<b%2d>...<b%2d>|", k, k - 7);
        if (k % 8 != 1) printf(" %6.4f ", p);
        else printf("%6.4f", p);

    }  //for (k= ..

    p = (double)kstest( xa, 32);

    printf("\n### Collision test is complete. Overall p = %6.4f\n",  p);

    free(t);

    return p;

}       /* collision() */


/*********************************************************************/
/**** Birthday spacing test ******************************************/

double bday( uint32_t (*rng)(void) )
{
    size_t i, j, k, m = 4096;
    uint32_t t[4096], obs[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double w, x, ex[11], v = 0., p;


    printf("\n### Birthday spacings test (4096 birthdays, 2^32 days/year) starts...\n");

    /* loop for sample of size 5000 */
    for (k = 1; k <= 5000; k++) {
        for (i = 0; i < m; i++)
            t[i] = rng();    /*   get 4096 IUNI's */

        qsort(t, m, sizeof(uint32_t), ucmpr); /*  sort them */

        for (i = m - 1; i > 0; i--)
            t[i] = t[i] - t[i - 1]; /* get spacings */

        qsort(t, m, sizeof(uint32_t), ucmpr);  /* sort spacings */

        j = 0;

        for (i = 1; i < m; i++)
            j = j + (t[i] == t[i - 1]) ;

        if ( j > 10 ) j = 10;

        obs[j]++;      /* count duplicate spacings */

    }   /* end  sample of 5000 loop */


    printf("           Table of Expected vs. Observed counts:\n");

    ex[0] = 5000. * exp(-4.);

    for (i = 1; i < 11; i++)
        ex[i] = 4. * ex[i - 1] / i;

    for (i = 0; i < 10; i++)
        printf( "%6.0f", i + 0.);

    printf("   >=10\n ");

    for (i = 0; i < 10; i++)
        printf( "%6.1f", ex[i]);

    printf( "%6.1f\n", ex[10]);

    for (i = 0; i < 10; i++)
        printf( "%6.0f", obs[i] + 0.0);

    printf( "%6.0f\n", obs[10] + 0.0);

    for( i = 0; i < 11; i++)
        v += (obs[i] - ex[i]) * (obs[i] - ex[i]) / ex[i];

    x = w = .5 * v;

    v = 24. + x * (24. + x * (12. + x * (4. + x) ) ); /* evaluate chisq probability */

    p = 1. - exp(-w) * v / 24.;


    printf("### Birthday spacings test is complete. Overall p = %6.4f\n",  p);

    return p;

}   /* bday() */



/*  GCD Test */
double gcd( uint32_t (*rng)(void) )
{
    const unsigned long n = 10000000;

    unsigned long u, v, w, i, k;

    unsigned long t[47], gcd[101];

    double s, ch32, ch99, e, r = 1.0e-10, p;

    unsigned long kp[36] = {
        0, 0, 0, 5506, 29532, 144541, 590691, 2065333, 6277273, 16797273, 39964829,
        85160313, 163520964, 284315128, 449367802, 647663189, 853365977, 1030017639,
        1140689999, 1160424679, 1085307890, 933657875, 738971259, 538010076, 360199303,
        221583256, 125137422, 64787412, 30697488, 13285746, 5238570, 1876493, 608568,
        177920, 46632, 13674
    };

    printf("\n### GCD test (%lu pairs of 32-bit integers) starts...\n", n);

    for (i = 1; i < 47; i++)
        t[i] = 0;

    for (i = 0; i < 101; i++)
        gcd[i] = 0;

    for (i = 1; i <= n; i++) {
        k = 0;
        do {
            u = rng();
            v = rng();
        } while (u == 0 || v == 0);

        do {
            w = u % v;
            u = v;
            v = w;
            k++;
        } while (v > 0);

        if (u > 100) u = 100;

        gcd[u]++;

        if (k < 3) k = 3;

        if (k > 35) k = 35;

        t[k]++;
    }

    ch32 = 0.0;

    for (i = 3; i <= 35; i++) {
        e = (r * kp[i]) * n;
        s = (t[i] - e) * (t[i] - e) / e;
        ch32 = ch32 + s;
    }

    printf("Euclid's algorithm:\n");

    printf("p-value, steps to gcd:   %f\n", p = Chisq2( 16., .5 * ch32));

    e = n * .61097691e-2;

    ch99 = (gcd[100] - e) * (gcd[100] - e) / e;

    for (i = 1; i < 100; i++) {
        e = n * 0.6079271 / ( i * i);
        s = (gcd[i] - e) * (gcd[i] - e) / e;
        ch99 = ch99 + s;
    }

    printf("p-value, dist. of gcd's: %f\n", Chisq2( 49.5, .5 * ch99));

    printf("### GCD test is complete. Overall p = %6.4f\n", p);

    return p;

} /* gcd() */


/* This function performs a frequency test on the random nos generated by rng().
   It counts the nos. of 0 and 1 and then conducts a Chi-square goodness-of-fit test.
   n is the no. of uint32_tegers (assume 32 bits) being tested. The function
   returns a uniform random no. in [0,1]. The given random nos. are failed if the
   return value is equal to or larger than 0.999 */
double frequency( uint32_t (*rng)(void) )
{

    const unsigned long m = 10000000;

    const uint32_t Nbit = 32; /* no. of bits in each random integer */

    uint32_t n, n0, n1, i, j, u;

    double chisq, p, nd2;


    printf("\n### Frequency test (on %lu 32-bit numbers) starts...\n", m);

    n0 = 0;

    for (i = 1; i <= m; ++i) {
        u = rng();

        for (j = 1; j <= Nbit; ++j) {
            if ((u & 1) == 0) ++n0;

            u = u / 2;
        }
    }

    n = m * Nbit; /* total no. of bits */
    n1 = n - n0;
    nd2 = n / 2.;

    chisq = (nd2 - n0) * (nd2 - n0) / nd2 + (nd2 - n1) * (nd2 - n1) / nd2;

    p = Chisq( 1, chisq);

    printf("### Frequency test (%d bits) is complete with p = %f\n", n, p);

    return p;

} /* freq */


/* This function performs the Maurer's universal test on the random nos stored in filename.
   The function returns a uniform random no. in [0,1]. n is the no. of 32 bit integers being
   tested. The given random nos. are failed if the return value is equal to or larger
   than 0.9999.
   CRYPYTO '90, U.M. Maurer, "A universal statistical test for random bit generators". */
double maurer(  uint32_t (*rng)(void)  )
{
    const uint32_t n = 10000000;
    uint32_t *t;
    uint32_t nl, i, l, q, u = 0, v;
    double sum, x, mu = 15.167379, sigma = sqrt(3.421), p;

    printf("\n### Maurer's test (on %u 32-bit numbers) starts...\n", n);

    l = 16;                             /* no. of bits inspected each round */
    nl = (uint32_t)lround(pow(2., l));  /* nl = 2^l */
    q = 8 * nl;                         /* no. of rounds in initializaton */

    t = (uint32_t *) malloc( nl * sizeof(uint32_t));

    for (i = 0; i < nl; ++i) {
        t[i] = 0;
    }

    for (i = 1; i <= q; ++i) {
        if ((i & 1) == 1) {
            u = rng();     /* get a fresh 32-bit no. */
            v = u >> 16;   /* use the high order 16 bits this round */
        } else {
            v = u & 65535; /* use the low order 16 bits left */
        }
        t[v] = i;          /* record when the bit pattern occurs */
    }

    sum = 0.;
    for (i = q + 1; i <= q + 2 * n; ++i) {
        if ((i & 1) == 1) {
            u = rng();
            v = u >> 16;
        } else {
            v = u & 65535;
        }

        sum = sum + log( i - t[v]);
        t[v] = i;
    }

    free( t);

    sum = sum / log(2.);

    x = sum / (2.*n);

    p = Phi( (x - mu) / sigma);

    printf("### Maurer's universal test (32x%d bits) is complete with p = %f\n", n, p);

    return p;

} /* maurer() */


/****** End of RNG Tests ****************************************************************/
