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

BOOST_AUTO_TEST_CASE(add)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { -1.0, 3.0, -4.5 };
    double expected[] = { 0.0, 5.0, -1.5 };
    v_add(a, b, 3);
    COMPARE_ARRAY(a, expected);
}

BOOST_AUTO_TEST_CASE(subtract)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { -1.0, 3.0, -4.5 };
    double expected[] = { 2.0, -1.0, 7.5 };
    v_subtract(a, b, 3);
    COMPARE_ARRAY(a, expected);
}

BOOST_AUTO_TEST_CASE(increment)
{
    double a[] = { -1.0, 3.0, -4.5 };
    double incr = -0.5;
    double expected[] = { -1.5, 2.5, -5.0 };
    v_increment(a, incr, 3);
    COMPARE_ARRAY(a, expected);
}

BOOST_AUTO_TEST_CASE(scale)
{
    double a[] = { -1.0, 3.0, -4.5 };
    double scale = -0.5;
    double expected[] = { 0.5, -1.5, 2.25 };
    v_scale(a, scale, 3);
    COMPARE_ARRAY(a, expected);
}

BOOST_AUTO_TEST_CASE(multiply)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { -1.0, 3.0, -4.5 };
    double expected[] = { -1.0, 6.0, -13.5 };
    v_multiply(a, b, 3);
    COMPARE_ARRAY(a, expected);
}

BOOST_AUTO_TEST_CASE(divide)
{
    double a[] = { 1.0, 2.0, 3.0 };
    double b[] = { -1.0, 3.0, -4.5 };
    double expected[] = { -1.0, 2.0/3.0, 3.0/-4.5 };
    v_divide(a, b, 3);
    COMPARE_ARRAY(a, expected);
}

BOOST_AUTO_TEST_SUITE_END()

