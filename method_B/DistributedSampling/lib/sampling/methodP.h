/******************************************************************************
 * methodP.h
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


#ifndef _METHOD_P_H_
#define _METHOD_P_H_

#include <vector>
#include "definitions.h"
#include "stocc.h"
#include "randomc.h"
#include "tools/crc_hash.h"
#include "tools/spooky_hash.h"
#include "sampling/methodR.h"

template <typename Stocc = StochasticLib1, typename LocalSampler = SeqDivideSampling<>, typename H = Spooky>
class ParDivideSampling {
    public:
        ParDivideSampling(SamplingConfig &config, ULONG seed, PEID size) 
            : config(config),
              stocc(seed),
              hash_seed(seed)
        { 
            // Compute input distribution
            rem = config.N % size;
            div = config.N / size;
        }

        template <typename F>
        void sample(ULONG n, 
                    PEID j, 
                    PEID k, 
                    PEID i,
                    F &&callback,
                    ULONG offset = 0) {
            if (j - k == 0) {
                ULONG h = H::hash(hash_seed + i);
                typename LocalSampler::base_type base_sampler(h, config.k);
                // Allocate hash table for base case

                LocalSampler local_sampler(base_sampler, config.k, h);
                local_sampler.sample(N(i+1) - N(i), n, callback, offset);
                return;
            } 
            
            ULONG m = (j + k) / 2;
            ULONG h = H::hash(hash_seed + j + k);
            stocc.RandomInit(h);
            ULONG N_split = N(m) - N(j-1);
            // TODO: Merge hashing with random variates to reduce initialization overhead
            ULONG x = stocc.Hypergeometric(N_split, n, N(k) - N(j-1)); 
            if (i < m) sample(x, j, m, i, callback, offset);
            else sample(n-x, m + 1, k, i, callback, offset + N_split);
        }

    private:
        SamplingConfig &config;
        Stocc stocc;
        ULONG hash_seed;
        ULONG div;
        PEID rem;

        inline ULONG N(PEID i) {
            return i * div + std::min(i, rem);
        }

};

#endif 

