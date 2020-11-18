//******************************************************************************
//
// File:    Gamma.cu
// Author:  Alan Kaminsky
// Version: 02-Feb-2012
//
// This source file is copyright (C) 2012 by Parallel Crypto LLC. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// alan.kaminsky@parallelcrypto.com.
//
// This source file is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 3 of the License, or (at your option) any
// later version.
//
// This source file is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

#ifndef __GAMMA_CU_INCLUDED__
#define __GAMMA_CU_INCLUDED__

#include "Util.cu"

//------------------------------------------------------------------------------
// This file contains CUDA functions for the gamma function and related
// functions. This file is intended to be #included into a program source file.

static int GAMMA_ITMAX = 200;
static double GAMMA_EPS = 2.22e-16;
static double GAMMA_FPMIN = (2.23e-308 / GAMMA_EPS);

/**
 * Returns the incomplete gamma function P(a,x), evaluated by its series
 * representation. Assumes a > 0 and x >= 0.
 */
static double gser
(double a,
    double x)
{
    double ap, del, sum;
    int i;

    ap = a;
    del = 1.0 / a;
    sum = del;
    for (i = 1; i <= GAMMA_ITMAX; ++i)
    {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (fabs(del) < fabs(sum) * GAMMA_EPS)
        {
            return sum * exp(-x + a * log(x) - lgamma(a));
        }
    }
    return 1.0; // Too many iterations
}

/**
 * Returns the complementary incomplete gamma function Q(a,x), evaluated by its
 * continued fraction representation. Assumes a > 0 and x >= 0.
 */
static double gcf
(double a,
    double x)
{
    double b, c, d, h, an, del;
    int i;

    b = x + 1.0 - a;
    c = 1.0 / GAMMA_FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1; i <= GAMMA_ITMAX; ++i)
    {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < GAMMA_FPMIN) d = GAMMA_FPMIN;
        c = b + an / c;
        if (fabs(c) < GAMMA_FPMIN) c = GAMMA_FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < GAMMA_EPS)
        {
            return exp(-x + a * log(x) - lgamma(a)) * h;
        }
    }
    return 0.0; // Too many iterations
}

/**
 * Returns the incomplete gamma function P(a,x).
 */
static double gammp
(double a,
    double x)
{
    if (a <= 0.0) die("gammp(): a illegal");
    if (x < 0.0) die("gammp(): x illegal");
    return x == 0.0 ? 0.0 : x < a + 1.0 ? gser(a, x) : 1.0 - gcf(a, x);
}

/**
 * Returns the complementary incomplete gamma function Q(a,x) = 1 - P(a,x).
 */
static double gammq
(double a,
    double x)
{
    if (a <= 0.0) die("gammq(): a illegal");
    if (x < 0.0) die("gammq(): x illegal");
    return x == 0.0 ? 1.0 : x < a + 1.0 ? 1.0 - gser(a, x) : gcf(a, x);
}

#endif