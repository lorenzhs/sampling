/*******************************************************************************
 * include/methodR.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef METHOD_R_HEADER
#define METHOD_R_HEADER

#include "definitions.hpp"
#include "hypergeometric_distribution.hpp"
#include "methodH.hpp"
#include "methodSH.hpp"

namespace sampling {

template <typename BaseSampler = HashSampling<>>
class SeqDivideSampling {
public:
    typedef BaseSampler base_type;

    SeqDivideSampling(BaseSampler &base_sampler, ULONG base_size, ULONG seed)
        : hyp(seed), base_sampler(std::move(base_sampler)), base_size(base_size) {
    }

    template <typename F>
    void sample(ULONG N, ULONG n, F &&callback, ULONG offset = 0) {
        if (n <= base_size) {
            base_sampler.sample(N, n, [&](ULONG sample) { callback(sample + offset); });
            return;
        }

        ULONG N_split = N / 2;
        // std::cout << "draw " << n << " from " << N << " with " << N_split << std::endl;
        ULONG x = hyp(n, N - n, N_split); //stocc.Hypergeometric(N_split, n, N);
        sample(N_split, x, callback, offset);
        sample(N - N_split, n - x, callback, offset + N_split);
    }

private:
    hypergeometric_distribution<ULONG> hyp;
    BaseSampler base_sampler;
    ULONG base_size;
};

}

#endif // METHOD_R_HEADER
