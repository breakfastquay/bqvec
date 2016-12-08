/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "bqvec/VectorOpsComplex.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace breakfastquay;

BOOST_AUTO_TEST_SUITE(TestVectorOpsComplex)

#define COMPARE_N(a, b, n)						\
    for (int cmp_i = 0; cmp_i < n; ++cmp_i) { \
        BOOST_CHECK_SMALL(a[cmp_i].re - b[cmp_i].re, 1e-14);		\
        BOOST_CHECK_SMALL(a[cmp_i].im - b[cmp_i].im, 1e-14);		\
    }

BOOST_AUTO_TEST_CASE(multiply)
{
    bq_complex_t a[] = { { 1.0, 2.0 }, { 3.0, -4.0 } };
    bq_complex_t b[] = { { -1.0, 3.0 }, { -4.5, 0.0 } };
    bq_complex_t o[2];
    bq_complex_t expected[] = { { -7.0, 1.0 }, { -13.5, 18.0 } };
    v_multiply(o, a, b, 2);
    COMPARE_N(o, expected, 2);
}

BOOST_AUTO_TEST_SUITE_END()

