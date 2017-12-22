/*******************************************************************************
 * tests/statistics_test.cpp
 *
 * Copyright (C) 2015 Timo Bingmann <tb@panthema.net>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <sampling/benchmark.hpp>

#include <gtest/gtest.h>

using namespace sampling;

TEST(Statistics, Simple) {
    statistics stats;

    // Test that empty statistics are handled correctly
    ASSERT_EQ(0.0, stats.avg());
    ASSERT_EQ(0.0, stats.stddev());

    stats.push(1.0);
    // Test that stddev() handles a single value
    ASSERT_EQ(0.0, stats.stddev());

    for (size_t i = 2; i <= 1000; ++i) {
        stats.push(1.0 / static_cast<double>(i));
    }

    ASSERT_EQ(1000, stats.count);
    ASSERT_DOUBLE_EQ(0.0074854708605503447, stats.avg());
    ASSERT_DOUBLE_EQ(0.039868430925506362, stats.stddev());
    ASSERT_DOUBLE_EQ(0.039848491723996423, stats.stddev(0));
}

/******************************************************************************/
