/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/* Copyright Chris Cannam - All Rights Reserved */

#ifndef BQ_VECTOR_OPS_H
#define BQ_VECTOR_OPS_H

#ifdef HAVE_IPP
#ifndef _MSC_VER
#include <inttypes.h>
#endif
#include <ipps.h>
#include <ippac.h>
#endif

#ifdef HAVE_VDSP
#include <vecLib/vDSP.h>
#include <vecLib/vForce.h>
#endif

#include <cstring>

#include "Restrict.h"

namespace breakfastquay {

/**
 * General principle:
 * 
 * Write basic vector-manipulation loops in such a way as to promote
 * the likelihood that a good current C++ compiler can auto-vectorize
 * them (e.g. gcc-4.x with -ftree-vectorize). Provide calls out to
 * supported vector libraries (e.g. IPP, Accelerate) where useful.
 * No intrinsics or assembly.
 *
 * Argument ordering:
 *
 * If a function operates on one or more vector sequences in memory,
 * they appear as pointers at the start of the argument list. If one
 * vector is the "destination", it appears first, with "source" second
 * if present (i.e. argument ordering follows memcpy).
 * 
 * The final argument is always a count of the number of elements in
 * each vector. 
 *
 * Some functions operate on a set of vectors at once: their names all
 * contain the text _channels, and the number of channels (i.e. number
 * of vectors) is the argument before last.
 *
 * Any additional arguments, e.g. scale factors, appear between the
 * vector pointers at the start and the counts at the end.
 *
 * The caller takes responsibility for ensuring that vector pointer
 * arguments are not aliased, i.e. that vectors provided as separate
 * arguments to the same function are distinct and do not overlap,
 * except where documented.
 */

/**
 * v_zero
 *
 * Zero the elements in the given vector, of length \arg count.
 */
template<typename T>
inline void v_zero(T *const R__ vec,
                   const int count)
{
    const T value = T(0);
    for (int i = 0; i < count; ++i) {
        vec[i] = value;
    }
}

#if defined HAVE_IPP
template<> 
inline void v_zero(float *const R__ vec, 
                   const int count)
{
    ippsZero_32f(vec, count);
}
template<> 
inline void v_zero(double *const R__ vec,
                   const int count)
{
    ippsZero_64f(vec, count);
}
#elif defined HAVE_VDSP
template<> 
inline void v_zero(float *const R__ vec, 
                   const int count)
{
    vDSP_vclr(vec, 1, count);
}
template<> 
inline void v_zero(double *const R__ vec,
                   const int count)
{
    vDSP_vclrD(vec, 1, count);
}
#endif

/**
 * v_zero_channels
 *
 * Zero the elements in the given set of \arg channels vectors, each
 * of length \arg count.
 */
template<typename T>
inline void v_zero_channels(T *const R__ *const R__ vec,
                            const int channels,
                            const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_zero(vec[c], count);
    }
}

/**
 * v_set
 *
 * Set all of the elements in the given vector, of length \arg count,
 * to \arg value.
 */
template<typename T>
inline void v_set(T *const R__ vec,
                  const T value,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        vec[i] = value;
    }
}

/**
 * v_copy
 *
 * Copy the contents of the vector \arg src to the vector \arg dst,
 * both of length \arg count.
 *
 * Caller guarantees that \arg src and \arg dst are non-overlapping.
 */
template<typename T>
inline void v_copy(T *const R__ dst,
                   const T *const R__ src,
                   const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = src[i];
    }
}

#if defined HAVE_IPP
template<>
inline void v_copy(float *const R__ dst,
                   const float *const R__ src,
                   const int count)
{
    ippsCopy_32f(src, dst, count);
}
template<>
inline void v_copy(double *const R__ dst,
                   const double *const R__ src,
                   const int count)
{
    ippsCopy_64f(src, dst, count);
}
#endif

/**
 * v_copy_channels
 *
 * Copy the contents of the individual vectors in the set \arg src to
 * the corresponding vectors in the set \arg dst. All vectors have
 * length \arg count and there are \arg channels vectors in each set.
 *
 * Caller guarantees that all of the \arg src and \arg dst vectors are
 * non-overlapping with each other.
 */
template<typename T>
inline void v_copy_channels(T *const R__ *const R__ dst,
                            const T *const R__ *const R__ src,
                            const int channels,
                            const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_copy(dst[c], src[c], count);
    }
}

/**
 * v_move
 *
 * Copy the contents of vector \arg src to the vector \arg dst, both
 * of length \arg count. The two vectors may overlap. (If you know
 * that they cannot overlap, use v_copy instead.)
 */
template<typename T>
inline void v_move(T *const dst,       // not R__ (aliased)
                   const T *const src, // not R__ (aliased)
                   const int count)
{
    memmove(dst, src, count * sizeof(T));
}

#if defined HAVE_IPP
template<>
inline void v_move(float *const dst,
                   const float *const src,
                   const int count)
{
    ippsMove_32f(src, dst, count);
}
template<>
inline void v_move(double *const dst,
                   const double *const src,
                   const int count)
{
    ippsMove_64f(src, dst, count);
}
#endif

/**
 * v_convert
 *
 * Copy the contents of vector \arg src to the vector \arg dst, both
 * of length \arg count. The two vectors may have different value
 * types, e.g. double and float. If they have the same type, this is
 * equivalent to v_copy.
 *
 * Caller guarantees that \arg src and \arg dst are non-overlapping.
 */
template<typename T, typename U>
inline void v_convert(U *const R__ dst,
                      const T *const R__ src,
                      const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = U(src[i]);
    }
}

template<>
inline void v_convert(float *const R__ dst,
                      const float *const R__ src,
                      const int count)
{
    v_copy(dst, src, count);
}
template<>
inline void v_convert(double *const R__ dst,
                      const double *const R__ src,
                      const int count)
{
    v_copy(dst, src, count);
}

#if defined HAVE_IPP
template<>
inline void v_convert(double *const R__ dst,
                      const float *const R__ src,
                      const int count)
{
    ippsConvert_32f64f(src, dst, count);
}
template<>
inline void v_convert(float *const R__ dst,
                      const double *const R__ src,
                      const int count)
{
    ippsConvert_64f32f(src, dst, count);
}
#elif defined HAVE_VDSP
template<>
inline void v_convert(double *const R__ dst,
                      const float *const R__ src,
                      const int count)
{
    vDSP_vspdp((float *)src, 1, dst, 1, count);
}
template<>
inline void v_convert(float *const R__ dst,
                      const double *const R__ src,
                      const int count)
{
    vDSP_vdpsp((double *)src, 1, dst, 1, count);
}
#endif

/**
 * v_convert_channels
 *
 * Copy the contents of the individual vectors in the set \arg src to
 * the corresponding vectors in the set \arg dst. All vectors have
 * length \arg count, and there are \arg channels vectors in each
 * set. The source and destination vectors may have different value
 * types, e.g. double and float. If they have the same type, this is
 * equivalent to v_copy_channels.
 *
 * Caller guarantees that all of the \arg src and \arg dst vectors are
 * non-overlapping with each other.
 */
template<typename T, typename U>
inline void v_convert_channels(U *const R__ *const R__ dst,
                               const T *const R__ *const R__ src,
                               const int channels,
                               const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_convert(dst[c], src[c], count);
    }
}

/**
 * v_add
 *
 * Add the elements in the vector \arg src to the corresponding
 * elements in the vector \arg dst, both of length arg \count, leaving
 * the result in \arg dst.
 *
 * Caller guarantees that \arg src and \arg dst are non-overlapping.
 */
template<typename T>
inline void v_add(T *const R__ dst,
                  const T *const R__ src,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] += src[i];
    }
}

/**
 * v_add
 *
 * Add the constant \arg value to every element of the vector \arg
 * dst, of length arg \count, leaving the result in \arg dst.
 */
template<typename T>
inline void v_add(T *const R__ dst,
                  const T value,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] += value;
    }
}

#if defined HAVE_IPP
template<>
inline void v_add(float *const R__ dst,
                  const float *const R__ src,
                  const int count)
{
    ippsAdd_32f_I(src, dst, count);
}    
inline void v_add(double *const R__ dst,
                  const double *const R__ src,
                  const int count)
{
    ippsAdd_64f_I(src, dst, count);
}    
#endif

/**
 * v_add_channels
 *
 * Add the elements in the individual vectors in the set \arg src to
 * the corresponding elements of the corresponding vectors in \arg
 * dst, leaving the results in \arg dst. All vectors have length \arg
 * count, and there are \arg channels vectors in each set.
 *
 * Caller guarantees that all of the \arg src and \arg dst vectors are
 * non-overlapping with each other.
 */
template<typename T>
inline void v_add_channels(T *const R__ *const R__ dst,
                           const T *const R__ *const R__ src,
                           const int channels, const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_add(dst[c], src[c], count);
    }
}

/**
 * v_add_with_gain
 *
 * Add the elements in the vector \arg src, each multiplied by the
 * constant factor \arg gain, to the corresponding elements in the
 * vector \arg dst, both of length arg \count, leaving the result in
 * \arg dst.
 *
 * Caller guarantees that \arg src and \arg dst are non-overlapping.
 */
template<typename T, typename G>
inline void v_add_with_gain(T *const R__ dst,
                            const T *const R__ src,
                            const G gain,
                            const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] += src[i] * gain;
    }
}

/**
 * v_add_channels_with_gain
 *
 * Add the elements in the individual vectors in the set \arg src,
 * each multiplied by the constant factor \arg gain, to the
 * corresponding elements of the corresponding vectors in \arg dst,
 * leaving the results in \arg dst. All vectors have length \arg
 * count, and there are \arg channels vectors in each set.
 *
 * Caller guarantees that all of the \arg src and \arg dst vectors are
 * non-overlapping with each other.
 */
template<typename T, typename G>
inline void v_add_channels_with_gain(T *const R__ *const R__ dst,
                                     const T *const R__ *const R__ src,
                                     const G gain,
                                     const int channels,
                                     const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_add_with_gain(dst[c], src[c], gain, count);
    }
}

template<typename T>
inline void v_subtract(T *const R__ dst,
                       const T *const R__ src,
                       const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] -= src[i];
    }
}

#if defined HAVE_IPP
template<>
inline void v_subtract(float *const R__ dst,
                       const float *const R__ src,
                       const int count)
{
    ippsSub_32f_I(src, dst, count);
}    
inline void v_subtract(double *const R__ dst,
                       const double *const R__ src,
                       const int count)
{
    ippsSub_64f_I(src, dst, count);
}    
#endif

template<typename T, typename G>
inline void v_scale(T *const R__ dst,
                    const G gain,
                    const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] *= gain;
    }
}

#if defined HAVE_IPP 
template<>
inline void v_scale(float *const R__ dst,
                    const float gain,
                    const int count)
{
    ippsMulC_32f_I(gain, dst, count);
}
template<>
inline void v_scale(double *const R__ dst,
                    const double gain,
                    const int count)
{
    ippsMulC_64f_I(gain, dst, count);
}
#endif

template<typename T>
inline void v_multiply(T *const R__ dst,
                       const T *const R__ src,
                       const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] *= src[i];
    }
}

#if defined HAVE_IPP 
template<>
inline void v_multiply(float *const R__ dst,
                       const float *const R__ src,
                       const int count)
{
    ippsMul_32f_I(src, dst, count);
}
template<>
inline void v_multiply(double *const R__ dst,
                       const double *const R__ src,
                       const int count)
{
    ippsMul_64f_I(src, dst, count);
}
#endif

template<typename T>
inline void v_multiply(T *const R__ dst,
                       const T *const R__ src1,
                       const T *const R__ src2,
                       const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = src1[i] * src2[i];
    }
}

template<typename T>
inline void v_divide(T *const R__ dst,
                     const T *const R__ src,
                     const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] /= src[i];
    }
}

#if defined HAVE_IPP 
template<>
inline void v_divide(float *const R__ dst,
                     const float *const R__ src,
                     const int count)
{
    ippsDiv_32f_I(src, dst, count);
}
template<>
inline void v_divide(double *const R__ dst,
                     const double *const R__ src,
                     const int count)
{
    ippsDiv_64f_I(src, dst, count);
}
#endif

#if defined HAVE_IPP 
template<>
inline void v_multiply(float *const R__ dst,
                       const float *const R__ src1,
                       const float *const R__ src2,
                       const int count)
{
    ippsMul_32f(src1, src2, dst, count);
}    
template<>
inline void v_multiply(double *const R__ dst,
                       const double *const R__ src1,
                       const double *const R__ src2,
                       const int count)
{
    ippsMul_64f(src1, src2, dst, count);
}
#endif

template<typename T>
inline void v_multiply_and_add(T *const R__ dst,
                               const T *const R__ src1,
                               const T *const R__ src2,
                               const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] += src1[i] * src2[i];
    }
}

#if defined HAVE_IPP
template<>
inline void v_multiply_and_add(float *const R__ dst,
                               const float *const R__ src1,
                               const float *const R__ src2,
                               const int count)
{
    ippsAddProduct_32f(src1, src2, dst, count);
}
template<>
inline void v_multiply_and_add(double *const R__ dst,
                               const double *const R__ src1,
                               const double *const R__ src2,
                               const int count)
{
    ippsAddProduct_64f(src1, src2, dst, count);
}
#endif

template<typename T>
inline T v_sum(const T *const R__ src,
               const int count)
{
    T result = T();
    for (int i = 0; i < count; ++i) {
        result += src[i];
    }
    return result;
}

template<typename T>
inline void v_log(T *const R__ dst,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = log(dst[i]);
    }
}

#if defined HAVE_IPP
template<>
inline void v_log(float *const R__ dst,
                  const int count)
{
    ippsLn_32f_I(dst, count);
}
template<>
inline void v_log(double *const R__ dst,
                  const int count)
{
    ippsLn_64f_I(dst, count);
}
#elif defined HAVE_VDSP
// no in-place vForce functions for these -- can we use the
// out-of-place functions with equal input and output vectors? can we
// use an out-of-place one with temporary buffer and still be faster
// than doing it any other way?
template<>
inline void v_log(float *const R__ dst,
                  const int count)
{
    float tmp[count];
    vvlogf(tmp, dst, &count);
    v_copy(dst, tmp, count);
}
template<>
inline void v_log(double *const R__ dst,
                  const int count)
{
    double tmp[count];
    vvlog(tmp, dst, &count);
    v_copy(dst, tmp, count);
}
#endif

template<typename T>
inline void v_exp(T *const R__ dst,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = exp(dst[i]);
    }
}

#if defined HAVE_IPP
template<>
inline void v_exp(float *const R__ dst,
                  const int count)
{
    ippsExp_32f_I(dst, count);
}
template<>
inline void v_exp(double *const R__ dst,
                  const int count)
{
    ippsExp_64f_I(dst, count);
}
#elif defined HAVE_VDSP
// no in-place vForce functions for these -- can we use the
// out-of-place functions with equal input and output vectors? can we
// use an out-of-place one with temporary buffer and still be faster
// than doing it any other way?
template<>
inline void v_exp(float *const R__ dst,
                  const int count)
{
    float tmp[count];
    vvexpf(tmp, dst, &count);
    v_copy(dst, tmp, count);
}
template<>
inline void v_exp(double *const R__ dst,
                  const int count)
{
    double tmp[count];
    vvexp(tmp, dst, &count);
    v_copy(dst, tmp, count);
}
#endif

template<typename T>
inline void v_sqrt(T *const R__ dst,
                   const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = sqrt(dst[i]);
    }
}

#if defined HAVE_IPP
template<>
inline void v_sqrt(float *const R__ dst,
                   const int count)
{
    ippsSqrt_32f_I(dst, count);
}
template<>
inline void v_sqrt(double *const R__ dst,
                   const int count)
{
    ippsSqrt_64f_I(dst, count);
}
#elif defined HAVE_VDSP
// no in-place vForce functions for these -- can we use the
// out-of-place functions with equal input and output vectors? can we
// use an out-of-place one with temporary buffer and still be faster
// than doing it any other way?
template<>
inline void v_sqrt(float *const R__ dst,
                   const int count)
{
    float tmp[count];
    vvsqrtf(tmp, dst, &count);
    v_copy(dst, tmp, count);
}
template<>
inline void v_sqrt(double *const R__ dst,
                   const int count)
{
    double tmp[count];
    vvsqrt(tmp, dst, &count);
    v_copy(dst, tmp, count);
}
#endif

template<typename T>
inline void v_square(T *const R__ dst,
                   const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = dst[i] * dst[i];
    }
}

#if defined HAVE_IPP
template<>
inline void v_square(float *const R__ dst,
                   const int count)
{
    ippsSqr_32f_I(dst, count);
}
template<>
inline void v_square(double *const R__ dst,
                   const int count)
{
    ippsSqr_64f_I(dst, count);
}
#endif

template<typename T>
inline void v_abs(T *const R__ dst,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = fabs(dst[i]);
    }
}

#if defined HAVE_IPP
template<>
inline void v_abs(float *const R__ dst,
                  const int count)
{
    ippsAbs_32f_I(dst, count);
}
template<>
inline void v_abs(double *const R__ dst,
                  const int count)
{
    ippsAbs_64f_I(dst, count);
}
#elif defined HAVE_VDSP
template<>
inline void v_abs(float *const R__ dst,
                  const int count)
{
    float tmp[count];
#if (MACOSX_DEPLOYMENT_TARGET <= 1070 && MAC_OS_X_VERSION_MIN_REQUIRED <= 1070)
    vvfabf(tmp, dst, &count);
#else
    vvfabsf(tmp, dst, &count);
#endif
    v_copy(dst, tmp, count);
}
#endif

template<typename T>
inline void v_interleave(T *const R__ dst,
                         const T *const R__ *const R__ src,
                         const int channels, 
                         const int count)
{
    int idx = 0;
    switch (channels) {
    case 2:
        // common case, may be vectorized by compiler if hardcoded
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < 2; ++j) {
                dst[idx++] = src[j][i];
            }
        }
        return;
    case 1:
        v_copy(dst, src[0], count);
        return;
    default:
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < channels; ++j) {
                dst[idx++] = src[j][i];
            }
        }
    }
}

#if defined HAVE_IPP 
template<>
inline void v_interleave(float *const R__ dst,
                         const float *const R__ *const R__ src,
                         const int channels, 
                         const int count)
{
    ippsInterleave_32f((const Ipp32f **)src, channels, count, dst);
}
// IPP does not (currently?) provide double-precision interleave
#endif

template<typename T>
inline void v_deinterleave(T *const R__ *const R__ dst,
                           const T *const R__ src,
                           const int channels, 
                           const int count)
{
    int idx = 0;
    switch (channels) {
    case 2:
        // common case, may be vectorized by compiler if hardcoded
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < 2; ++j) {
                dst[j][i] = src[idx++];
            }
        }
        return;
    case 1:
        v_copy(dst[0], src, count);
        return;
    default:
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < channels; ++j) {
                dst[j][i] = src[idx++];
            }
        }
    }
}

#if defined HAVE_IPP
template<>
inline void v_deinterleave(float *const R__ *const R__ dst,
                           const float *const R__ src,
                           const int channels, 
                           const int count)
{
    ippsDeinterleave_32f((const Ipp32f *)src, channels, count, (Ipp32f **)dst);
}
// IPP does not (currently?) provide double-precision deinterleave
#endif

template<typename T>
inline void v_fftshift(T *const R__ vec,
                       const int count)
{
    const int hs = count/2;
    for (int i = 0; i < hs; ++i) {
        T t = vec[i];
        vec[i] = vec[i + hs];
        vec[i + hs] = t;
    }
}

template<typename T>
inline T v_mean(const T *const R__ vec, const int count)
{
    T t = T(0);
    for (int i = 0; i < count; ++i) {
        t += vec[i];
    }
    t /= T(count);
    return t;
}

template<typename T>
inline T v_mean_channels(const T *const R__ *const R__ vec,
                         const int channels,
                         const int count)
{
    T t = T(0);
    for (int c = 0; c < channels; ++c) {
        t += v_mean(vec[c], count);
    }
    t /= T(channels);
    return t;
}

}

#endif
