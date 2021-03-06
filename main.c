/*
*    Project:    DIEAD
*    File:       main.c
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
#include "main.h"

#define TEST_RNG    sapparot2

/* ============================================================================
*  Implement your generator as a parameterless function which returns a 32-bit
*  unsigned integer for a pseudorandom value.
*
*  Modify TEST_RNG to point at your function to test.
*/

/* -----------------------------------------------------------------------------
*  sapparot2()
*  This sample function implements the Sapparot-2 pseudorandom number generator
*  as described at http://literatecode.com/sapparot2
*/
static uint32_t sapparot2(void)
{
    static uint32_t a = 0, b = 0, c = 0; /* {a,b,c} is a state */
    register uint32_t m;

    c += a;
    c = R(c, b >> 27);
    b = (b + ((a << 1) + 1)) ^ R(b, 5);
    a += 0x9e3779b9;
    a = R(a, 7);
    m = a;
    a = b;
    b = m;

    return (c ^ b ^ a);
} /* sapparot2 */


/* -------------------------------------------------------------------------- */
int main(void)
{
    struct {
        char *name;
        double p;
    } t[] = {
        {.name = "Maurer's",         .p = maurer(TEST_RNG)},
        {.name = "Frequency",        .p = frequency(TEST_RNG)},
        {.name = "GCD",              .p = gcd(TEST_RNG)},
        {.name = "Birthday spacing", .p = bday(TEST_RNG)},
        {.name = "Gorilla",          .p = gorilla(TEST_RNG)},
        {.name = "Collision",        .p = collision(TEST_RNG)}
    };
    const size_t tests = sizeof(t) / sizeof(t[0]);

    printf("\nTests Summary:\n");
    for (size_t i = 0; i < tests; i++) {
        printf("* %s test: p = %6.4f\n", t[i].name, t[i].p);
    }

    return 0;
} /* main */
