//******************************************************************************
//
// File:    Util.cu
// Author:  Alan Kaminsky
// Version: 19-Jan-2012
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

#ifndef __UTIL_CU_INCLUDED__
#define __UTIL_CU_INCLUDED__

#include <stdlib.h>
#include <stdio.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include<math.h>
#include <sys/types.h>
#include <time.h>
#include <cuda_runtime_api.h>

/**
 * Program name (argv[0]).
 */
static char* progname;

/**
 * Print an error message and exit.
 *
 * @param  msg  Error message.
 */
static void die
(char* msg)
{
    fprintf(stderr, "%s: %s\n", progname, msg);
    exit(1);
}

/**
 * If necessary, print a CUDA related error message and exit.
 *
 * @param  err  CUDA error.
 * @param  msg  Error message.
 */
static void checkCuda
(cudaError_t err,
    char* msg)
{
    if (err != cudaSuccess)
    {
        fprintf(stderr, "%s: %s: %s (%d)\n",
            progname, msg, cudaGetErrorString(err), err);
        exit(1);
    }
}

/**
 * Set the CUDA device. The CUDA_DEVICE environment variable specifies the CUDA
 * device to use. If this variable is not set, CUDA device 0 is used. The CUDA
 * device must support compute capability 2.0 or higher.
 */
static void setCudaDevice()
{
    char* CUDA_DEVICE = getenv("CUDA_DEVICE");
    int dev;
    struct cudaDeviceProp prop;

    if (CUDA_DEVICE == NULL) CUDA_DEVICE = "0";
    if (sscanf(CUDA_DEVICE, "%d", &dev) != 1)
    {
        fprintf(stderr,
            "%s: Environment variable CUDA_DEVICE=\"%s\" invalid\n",
            progname, CUDA_DEVICE);
        exit(1);
    }
    if (cudaGetDeviceProperties(&prop, dev) != cudaSuccess)
    {
        fprintf(stderr,
            "%s: Could not get properties for CUDA device %d\n",
            progname, dev);
        exit(1);
    }
    if (prop.major < 2 || prop.major == 9999)
    {
        fprintf(stderr,
            "%s: CUDA device %d: %s, compute capability %d.%d, 2.0 required\n",
            progname, dev, prop.name, prop.major, prop.minor);
        exit(1);
    }
    printf("CUDA device %d: %s, compute capability %d.%d\n",
        dev, prop.name, prop.major, prop.minor);
    checkCuda(cudaSetDevice(dev), "Could not set CUDA device");
}
#ifndef _MSC_VER
/**
 * Returns the system clock in milliseconds.
 */
static u_int64_t currentTimeMillis()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000LLU + tv.tv_usec / 1000LLU;
}
#endif
#endif