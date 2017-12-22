/*******************************************************************************
 * tests/hypergeometric_distribution_test.cpp
 *
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <sampling/hypergeometric_distribution.hpp>

#include <gtest/gtest.h>

using namespace sampling;

TEST(Hypergeometric, Simple) {
    hypergeometric hyp(42);

    // Test a known sequence
    ASSERT_EQ(9, hyp(10, 6, 14));
    ASSERT_EQ(8, hyp(10, 6, 14));
    ASSERT_EQ(9, hyp(10, 6, 14));
    ASSERT_EQ(9, hyp(10, 6, 14));
    ASSERT_EQ(10, hyp(10, 6, 14));
    ASSERT_EQ(9, hyp(10, 6, 14));
    ASSERT_EQ(8, hyp(10, 6, 14));
    ASSERT_EQ(10, hyp(10, 6, 14));
    ASSERT_EQ(8, hyp(10, 6, 14));
    ASSERT_EQ(9, hyp(10, 6, 14));

    // Test with nbad = 0
    ASSERT_EQ(7, hyp(10, 0, 7));
    ASSERT_EQ(7, hyp(10, 0, 7));
    ASSERT_EQ(7, hyp(10, 0, 7));
    ASSERT_EQ(7, hyp(10, 0, 7));

    ASSERT_EQ(42, hyp(100, 0, 42));
    ASSERT_EQ(42, hyp(100, 0, 42));
    ASSERT_EQ(42, hyp(100, 0, 42));
    ASSERT_EQ(42, hyp(100, 0, 42));

    // Test with ngood = 0
    ASSERT_EQ(0, hyp(0, 10, 5));
    ASSERT_EQ(0, hyp(0, 10, 5));
    ASSERT_EQ(0, hyp(0, 10, 5));
    ASSERT_EQ(0, hyp(0, 10, 5));
}

/******************************************************************************/
