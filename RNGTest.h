/*
*    Project:    DIEAD
*    File:       RNGTest.h
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
 *  Filename: RNGTest.h
 */
#ifndef RNG_TEST_H_
#define RNG_TEST_H_



/* Gorilla test                                                                         */
double gorilla( uint32_t (*rng)(void) );
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
double pcoll(uint32_t m, uint32_t n, uint32_t c);

/* collision() conducts the collision test on the 0th to 31th bit sequences of the
   words generated by the 32-bit generator rng. j=0 indicates the sequence formed from
   the least significant bits of the words and j=31 indicates the most significant bit
   sequence. The test throws 2^20 balls into 2^20 urns. The destination
   of each ball is determined by log_2 m bits extracted from the bit sequence. The bits are
   obtained from the jth bit of the 32-bit word returned by calling rng(). (See Tsang's paper
   "Tuning the collision test for stringency" for details if needed.)
*/
double collision( uint32_t (*rng)(void) );

/*********************************************************************/
/**** Birthday spacing test ******************************************/

double bday( uint32_t (*rng)(void) );

/*  GCD Test */
double gcd( uint32_t (*rng)(void) );

/* This function performs a frequency test on the random nos generated by rng().
   It counts the nos. of 0 and 1 and then conducts a Chi-square goodness-of-fit test.
   n is the no. of unsigned integers (assume 32 bits) being tested. The function
   returns a uniform random no. in [0,1]. The given random nos. are failed if the
   return value is equal to or larger than 0.999 */
double frequency( uint32_t (*rng)(void) );

/* This function performs the Maurer's universal test on the random nos stored in filename.
   The function returns a uniform random no. in [0,1]. n is the no. of 32 bit integers being
   tested. The given random nos. are failed if the return value is equal to or larger
   than 0.9999.
   CRYPYTO '90, U.M. Maurer, "A universal statistical test for random bit generators". */
double maurer(  uint32_t (*rng)(void) );



#endif /* RNG_TEST_H_ */
