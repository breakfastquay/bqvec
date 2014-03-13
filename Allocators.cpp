/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/* Copyright Chris Cannam - All Rights Reserved */

#include "Allocators.h"

#ifdef HAVE_IPP
#include <ipps.h>
#endif

#include <iostream>
using std::cerr;
using std::endl;

namespace breakfastquay {

#ifdef HAVE_IPP

template <>
float *allocate(size_t count)
{
    float *ptr = ippsMalloc_32f(count);
    if (!ptr) throw (std::bad_alloc());
    return ptr;
}

template <>
double *allocate(size_t count)
{
    double *ptr = ippsMalloc_64f(count);
    if (!ptr) throw (std::bad_alloc());
    return ptr;
}

template <>
void deallocate(float *ptr)
{
    if (ptr) ippsFree((void *)ptr);
}

template <>
void deallocate(double *ptr)
{
    if (ptr) ippsFree((void *)ptr);
}

#endif

}

