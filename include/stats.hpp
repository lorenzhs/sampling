/*******************************************************************************
 * include/stats.hpp
 *
 * Copyright (C) 2016-2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cmath>

namespace sampling {

// template-based loop unrolling
template <size_t N> struct faux_unroll {
    template <typename F> static void call(F &&f) {
        faux_unroll<N-1>::call(f);
        f(N-1);
    }
};

template <> struct faux_unroll<0> {
    template <typename F> static void call(F &&) {}
};


struct statistics {
    // Single-pass standard deviation calculation as described in Donald Knuth:
    // The Art of Computer Programming, Volume 2, Chapter 4.2.2, Equations 15&16
    double mean;
    double nvar; // approx n * variance; stddev = sqrt(nvar / (count-1))
    size_t count;

    statistics() : mean(0.0), nvar(0.0), count(0) {}

    void push(double t) {
        ++count;
        if (count == 1) {
            mean = t;
        } else {
            double oldmean = mean;
            mean += (t - oldmean) / count;
            nvar += (t - oldmean) * (t - mean);
        }
    }

    double avg() {
        return mean;
    }
    double stddev() {
        //assert(count > 1);
        return sqrt(nvar / (count - 1));
    }
};

struct global_stats {
    void push_sum(double t) { s_sum.push(t); }
    void push_gen(double t) { s_gen.push(t); }
    void push_prefsum(double t) { s_prefsum.push(t); }
    void push_fix(double t) { s_fix.push(t); }
    statistics s_sum, s_gen, s_prefsum, s_fix;
};

/*
// Print SIMD registers
void print128_num(__m128i var) {
    uint32_t *val = (uint32_t*) &var;
    printf("Numerical: %i %i %i %i \n", val[0], val[1], val[2], val[3]);
}

void print256_num(__m256i var) {
    uint32_t *val = (uint32_t*) &var;
    printf("Numerical: %i %i %i %i %i %i %i %i \n",
           val[0], val[1], val[2], val[3], val[4], val[5],
           val[6], val[7]);
}
*/

} // namespace sampling
