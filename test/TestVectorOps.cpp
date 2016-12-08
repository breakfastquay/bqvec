/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "bqvec/VectorOps.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace breakfastquay;

BOOST_AUTO_TEST_SUITE(TestVectorOps)

#define COMPARE_ARRAY(a, b)						\
    for (int cmp_i = 0; cmp_i < (int)(sizeof(a)/sizeof(a[0])); ++cmp_i) { \
        BOOST_CHECK_SMALL(a[cmp_i] - b[cmp_i], 1e-14);			\
    }

#define COMPARE_N(a, b, n)						\
    for (int cmp_i = 0; cmp_i < n; ++cmp_i) { \
        BOOST_CHECK_SMALL(a[cmp_i] - b[cmp_i], 1e-14);			\
    }

BOOST_AUTO_TEST_CASE(add)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { -1.0, 3.0, -4.5 };
    double expected[] = { 0.0, 5.0, -1.5 };
    v_add(a, b, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(subtract)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { -1.0, 3.0, -4.5 };
    double expected[] = { 2.0, -1.0, 7.5 };
    v_subtract(a, b, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(increment)
{
    double a[] = { -1.0, 3.0, -4.5 };
    double incr = -0.5;
    double expected[] = { -1.5, 2.5, -5.0 };
    v_increment(a, incr, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(scale)
{
    double a[] = { -1.0, 3.0, -4.5 };
    double scale = -0.5;
    double expected[] = { 0.5, -1.5, 2.25 };
    v_scale(a, scale, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(multiply)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { -1.0, 3.0, -4.5 };
    double expected[] = { -1.0, 6.0, -13.5 };
    v_multiply(a, b, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(divide)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { -1.0, 3.0, -4.5 };
    double expected[] = { -1.0, 2.0/3.0, 3.0/-4.5 };
    v_divide(a, b, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(sum)
{
    double a[] = { 1.0, 2.0, -3.5 };
    double s = v_sum(a, 3);
    BOOST_CHECK_EQUAL(s, -0.5);
}

BOOST_AUTO_TEST_CASE(multiply_and_sum)
{
    double a[] = { 2.0, 0.0, -1.5 };
    double b[] = { 3.0, 4.0, 5.0 };
    double s = v_multiply_and_sum(a, b, 3);
    BOOST_CHECK_EQUAL(s, -1.5);
}

BOOST_AUTO_TEST_CASE(log)
{
    double a[] = { 1.0, 1.0 / M_E, M_E };
    double expected[] = { 0.0, -1.0, 1.0 };
    v_log(a, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(exp)
{
    double a[] = { 0.0, -1.0, 2.0 };
    double expected[] = { 1.0, 1.0 / M_E, M_E * M_E };
    v_exp(a, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(sqrt)
{
    double a[] = { 0.0, 1.0, 4.0 };
    double expected[] = { 0.0, 1.0, 2.0 };
    v_sqrt(a, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(square)
{
    double a[] = { 0.0, 1.5, -2.0 };
    double expected[] = { 0.0, 2.25, 4.0 };
    v_square(a, 3);
    COMPARE_N(a, expected, 3);
}

BOOST_AUTO_TEST_CASE(abs)
{
    double a[] = { -1.9, 0.0, 0.01, -0.0 };
    double expected[] = { 1.9, 0.0, 0.01, 0.0 };
    v_abs(a, 4);
    COMPARE_N(a, expected, 4);
}

BOOST_AUTO_TEST_CASE(mean)
{
    double a[] = { -1.0, 1.6, 3.0 };
    double s = v_mean(a, 3);
    BOOST_CHECK_EQUAL(s, 1.2);
}

BOOST_AUTO_TEST_CASE(interleave_1)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double *ch[] = { a };
    double o[3];
    double expected[] = { 1.0, 2.0, 3.0 };
    v_interleave(o, ch, 1, 3);
    COMPARE_N(o, expected, 3);
}

BOOST_AUTO_TEST_CASE(interleave_2)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { 4.0, 5.0, 6.0 };
    double *ch[] = { a, b };
    double o[6];
    double expected[] = { 1.0, 4.0, 2.0, 5.0, 3.0, 6.0 };
    v_interleave(o, ch, 2, 3);
    COMPARE_N(o, expected, 6);
}

BOOST_AUTO_TEST_CASE(interleave_3)
{
    double a[] = { 1.0, 2.0 };
    double b[] = { 3.0, 4.0 };
    double c[] = { 5.0, 6.0 };
    double *ch[] = { a, b, c };
    double o[6];
    double expected[] = { 1.0, 3.0, 5.0, 2.0, 4.0, 6.0 };
    v_interleave(o, ch, 3, 2);
    COMPARE_N(o, expected, 6);
}

BOOST_AUTO_TEST_CASE(deinterleave_1)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double o[3];
    double *oo[] = { o };
    double *expected[] = { a };
    v_deinterleave(oo, a, 1, 3);
    COMPARE_N(oo[0], expected[0], 3);
}

BOOST_AUTO_TEST_CASE(deinterleave_2)
{
    double a[] = { 1.0, 4.0, 2.0, 5.0, 3.0, 6.0 };
    double o1[3], o2[3];
    double *oo[] = { o1, o2 };
    double e1[] = { 1.0, 2.0, 3.0 }, e2[] = { 4.0, 5.0, 6.0 };
    double *expected[] = { e1, e2 };
    v_deinterleave(oo, a, 2, 3);
    COMPARE_N(oo[0], expected[0], 3);
    COMPARE_N(oo[1], expected[1], 3);
}

BOOST_AUTO_TEST_CASE(deinterleave_3)
{
    double a[] = { 1.0, 3.0, 5.0, 2.0, 4.0, 6.0 };
    double o1[2], o2[2], o3[2];
    double *oo[] = { o1, o2, o3 };
    double e1[] = { 1.0, 2.0 }, e2[] = { 3.0, 4.0 }, e3[] = { 5.0, 6.0 };
    double *expected[] = { e1, e2, e3 };
    v_deinterleave(oo, a, 3, 2);
    COMPARE_N(oo[0], expected[0], 2);
    COMPARE_N(oo[1], expected[1], 2);
    COMPARE_N(oo[2], expected[2], 2);
}

BOOST_AUTO_TEST_CASE(mix_1)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double *ch[] = { a };
    double o[3];
    double expected[] = { 1.0, 2.0, 3.0 };
    v_mix(o, ch, 1, 3);
    COMPARE_N(o, expected, 3);
}

BOOST_AUTO_TEST_CASE(mix_2)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { 4.0, 5.0, 6.0 };
    double *ch[] = { a, b };
    double o[6];
    double expected[] = { 2.5, 3.5, 4.5 };
    v_mix(o, ch, 2, 3);
    COMPARE_N(o, expected, 3);
}

BOOST_AUTO_TEST_CASE(mix_3)
{
    double a[] = { 1.0, 2.0 };
    double b[] = { 3.0, 4.0 };
    double c[] = { 5.0, 6.0 };
    double *ch[] = { a, b, c };
    double o[6];
    double expected[] = { 3.0, 4.0 };
    v_mix(o, ch, 3, 2);
    COMPARE_N(o, expected, 2);
}

BOOST_AUTO_TEST_CASE(reconfigure_1_2)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double *aa[] = { a };
    double o1[3], o2[3];
    double *oo[] = { o1, o2 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 1.0, 2.0, 3.0 };
    double *expected[] = { e1, e2 };
    v_reconfigure_channels(oo, 2, aa, 1, 3);
    COMPARE_N(oo[0], expected[0], 3);
    COMPARE_N(oo[1], expected[1], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_2_1)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[] = { 4.0, 5.0, 6.0 };
    double *aa[] = { a1, a2 };
    double o1[3];
    double *oo[] = { o1 };
    double e1[] = { 2.5, 3.5, 4.5 };
    double *expected[] = { e1 };
    v_reconfigure_channels(oo, 1, aa, 2, 3);
    COMPARE_N(oo[0], expected[0], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_3_1)
{
    double a1[] = { 1.0, 2.0 };
    double a2[] = { 3.0, 4.0 };
    double a3[] = { 5.0, 6.0 };
    double *aa[] = { a1, a2, a3 };
    double o1[2];
    double *oo[] = { o1 };
    double e1[] = { 3.0, 4.0 };
    double *expected[] = { e1 };
    v_reconfigure_channels(oo, 1, aa, 3, 2);
    COMPARE_N(oo[0], expected[0], 2);
}

BOOST_AUTO_TEST_CASE(reconfigure_1_3)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double *aa[] = { a };
    double o1[3], o2[3], o3[3];
    double *oo[] = { o1, o2, o3 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 1.0, 2.0, 3.0 };
    double e3[] = { 1.0, 2.0, 3.0 };
    double *expected[] = { e1, e2, e3 };
    v_reconfigure_channels(oo, 3, aa, 1, 3);
    COMPARE_N(oo[0], expected[0], 3);
    COMPARE_N(oo[1], expected[1], 3);
    COMPARE_N(oo[2], expected[2], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_2_3)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[] = { 4.0, 5.0, 6.0 };
    double *aa[] = { a1, a2 };
    double o1[3], o2[3], o3[3];
    double *oo[] = { o1, o2, o3 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 4.0, 5.0, 6.0 };
    double e3[] = { 0.0, 0.0, 0.0 };
    double *expected[] = { e1, e2, e3 };
    v_reconfigure_channels(oo, 3, aa, 2, 3);
    COMPARE_N(oo[0], expected[0], 3);
    COMPARE_N(oo[1], expected[1], 3);
    COMPARE_N(oo[2], expected[2], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_3_2)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[] = { 4.0, 5.0, 6.0 };
    double a3[] = { 7.0, 8.0, 9.0 };
    double *aa[] = { a1, a2, a3 };
    double o1[3], o2[3];
    double *oo[] = { o1, o2 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 4.0, 5.0, 6.0 };
    double *expected[] = { e1, e2 };
    v_reconfigure_channels(oo, 2, aa, 3, 3);
    COMPARE_N(oo[0], expected[0], 3);
    COMPARE_N(oo[1], expected[1], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_3_3)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[] = { 4.0, 5.0, 6.0 };
    double a3[] = { 7.0, 8.0, 9.0 };
    double *aa[] = { a1, a2, a3 };
    double o1[3], o2[3], o3[3];
    double *oo[] = { o1, o2, o3 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 4.0, 5.0, 6.0 };
    double e3[] = { 7.0, 8.0, 9.0 };
    double *expected[] = { e1, e2, e3 };
    v_reconfigure_channels(oo, 3, aa, 3, 3);
    COMPARE_N(oo[0], expected[0], 3);
    COMPARE_N(oo[1], expected[1], 3);
    COMPARE_N(oo[2], expected[2], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_1_2_inplace)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[3];
    double *aa[] = { a1, a2 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 1.0, 2.0, 3.0 };
    double *expected[] = { e1, e2 };
    v_reconfigure_channels_inplace(aa, 2, 1, 3);
    COMPARE_N(aa[0], expected[0], 3);
    COMPARE_N(aa[1], expected[1], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_2_1_inplace)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[] = { 4.0, 5.0, 6.0 };
    double *aa[] = { a1, a2 };
    double e1[] = { 2.5, 3.5, 4.5 };
    double *expected[] = { e1 };
    v_reconfigure_channels_inplace(aa, 1, 2, 3);
    COMPARE_N(aa[0], expected[0], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_3_1_inplace)
{
    double a1[] = { 1.0, 2.0 };
    double a2[] = { 3.0, 4.0 };
    double a3[] = { 5.0, 6.0 };
    double *aa[] = { a1, a2, a3 };
    double e1[] = { 3.0, 4.0 };
    double *expected[] = { e1 };
    v_reconfigure_channels_inplace(aa, 1, 3, 2);
    COMPARE_N(aa[0], expected[0], 2);
}

BOOST_AUTO_TEST_CASE(reconfigure_1_3_inplace)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[3], a3[3];
    double *aa[] = { a1, a2, a3 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 1.0, 2.0, 3.0 };
    double e3[] = { 1.0, 2.0, 3.0 };
    double *expected[] = { e1, e2, e3 };
    v_reconfigure_channels_inplace(aa, 3, 1, 3);
    COMPARE_N(aa[0], expected[0], 3);
    COMPARE_N(aa[1], expected[1], 3);
    COMPARE_N(aa[2], expected[2], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_2_3_inplace)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[] = { 4.0, 5.0, 6.0 };
    double a3[3];
    double *aa[] = { a1, a2, a3 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 4.0, 5.0, 6.0 };
    double e3[] = { 0.0, 0.0, 0.0 };
    double *expected[] = { e1, e2, e3 };
    v_reconfigure_channels_inplace(aa, 3, 2, 3);
    COMPARE_N(aa[0], expected[0], 3);
    COMPARE_N(aa[1], expected[1], 3);
    COMPARE_N(aa[2], expected[2], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_3_2_inplace)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[] = { 4.0, 5.0, 6.0 };
    double a3[] = { 7.0, 8.0, 9.0 };
    double *aa[] = { a1, a2, a3 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 4.0, 5.0, 6.0 };
    double *expected[] = { e1, e2 };
    v_reconfigure_channels_inplace(aa, 2, 3, 3);
    COMPARE_N(aa[0], expected[0], 3);
    COMPARE_N(aa[1], expected[1], 3);
}

BOOST_AUTO_TEST_CASE(reconfigure_3_3_inplace)
{
    double a1[] = { 1.0, 2.0, 3.0 };
    double a2[] = { 4.0, 5.0, 6.0 };
    double a3[] = { 7.0, 8.0, 9.0 };
    double *aa[] = { a1, a2, a3 };
    double e1[] = { 1.0, 2.0, 3.0 };
    double e2[] = { 4.0, 5.0, 6.0 };
    double e3[] = { 7.0, 8.0, 9.0 };
    double *expected[] = { e1, e2, e3 };
    v_reconfigure_channels_inplace(aa, 3, 3, 3);
    COMPARE_N(aa[0], expected[0], 3);
    COMPARE_N(aa[1], expected[1], 3);
    COMPARE_N(aa[2], expected[2], 3);
}

BOOST_AUTO_TEST_CASE(fftshift)
{
    double a[] = { 0.1, 2.0, -0.3, 4.0 };
    double e[] = { -0.3, 4.0, 0.1, 2.0 };
    v_fftshift(a, 4);
    COMPARE_N(a, e, 4);
}

BOOST_AUTO_TEST_SUITE_END()

