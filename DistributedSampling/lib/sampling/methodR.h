/******************************************************************************
 * methodR.h
 *
 * Source of the sampling routine
 ******************************************************************************
 * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/


#ifndef _METHOD_R_H_
#define _METHOD_R_H_

#include "definitions.h"
#include "stocc.h"
#include "var_gen.h"
#include "sampling/methodD.h"
#include "sampling/methodH.h"
#include "sampling/methodSH.h"

template <typename Stocc = StochasticLib1, typename BaseSampler = HashSampling<>>
class SeqDivideSampling {
    public:
        typedef BaseSampler base_type;

        SeqDivideSampling(SamplingConfig &config, BaseSampler &base_sampler, ULONG base_size, ULONG seed) 
            : config(config),
              stocc(seed),
              base_sampler(std::move(base_sampler)),
              base_size(base_size)
        { }

        template <typename F>
        void sample(ULONG N, 
                    ULONG n, 
                    F &&callback,
                    ULONG offset = 0) {
            if (n <= base_size) {
                base_sampler.sample(N, n, [&](ULONG sample) { callback(sample + offset); });
                return;
            } 
            
            ULONG N_split = N/2;
            ULONG x = 0;
            if (config.use_binom)
                x = stocc.Binomial(n, (double)N_split/N); 
            else 
                x = stocc.Hypergeometric(n, N_split, N); 
            sample(N_split, x, callback, offset);
            sample(N-N_split, n-x, callback, offset + N_split);
        }

    private:
        SamplingConfig &config;
        Stocc stocc;
        BaseSampler base_sampler;
        ULONG base_size;
};

#endif 

