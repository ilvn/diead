/*
*    Project:    DIEAD
*    File:       RNGUtil.h
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
#ifndef RNG_UTIL_H_
#define RNG_UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>


typedef uint32_t   uniform;
typedef uint64_t   counter;
typedef double     real;

#define DIM        4096
#define UNIMAX     4294967296.   /*pow(2,32)*/
#undef PI
#define PI         3.141592653589793

#define MAX(a, b)  ( (a)<(b) ? (b):(a) )
#define MIN(a, b)  ( (a)<(b) ? (a):(b) )
#define ABS(x)     ( (x)>0 ? (x):(-(x)) )
#define SIGN(x)    ( (x)>0 ? 1:(-1) )


uniform uni(char *filename);

double Phi(double z);
double Chisq(int df, double chsq);
double Chisq2( double d, double z);
double Poisson(double lambda, int k);
float kstest(float[], int); /* Anderson-Darling version, K-S test.*/
int    comp(const void *i, const void *j);
int    ucmpr( const void *i1, const void *i2);

#endif /* RNG_UTIL_H_ */
