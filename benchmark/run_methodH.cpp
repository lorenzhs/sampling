/*******************************************************************************
 * benchmark/run_methodH.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "arg_parser.hpp"
#include "parse_parameters.hpp"

#include <sampling/benchmark.hpp>
#include <sampling/methodB.hpp>
#include <sampling/methodH.hpp>
#include <sampling/sampling_config.hpp>
#include <sampling/timer.hpp>
#include <sampling/rng/mkl.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace sampling;

int main(int argn, char **argv) {
    // Read command-line args
    SamplingConfig config;
    parse_parameters(argn, argv, config);

    // Main algorithm
    FILE *fp;
    std::cout << "sample (n=" << config.n <<
        ", N=" << config.N <<
        ", k=" << config.k <<
        ", s=" << config.seed <<
        ", i=" << config.iterations << ")" << std::endl;
    std::string filename = config.output_file;
    fp = fopen(filename.c_str(), "w+");

    // Resulting samples
    std::vector<ULONG> sample;
    sample.reserve(config.n);
    // ULONG samples_taken = 0;

    // Statistics
    timer t;
    statistics stats;

    std::cout << "warmup" << std::endl;
    for (ULONG iteration = 0; iteration < std::min((ULONG)100, config.iterations); ++iteration) {
        sample.clear();
        // samples_taken = 0;

        // Compute sample
        HashSampling<> hs(iteration, config.n);
        hs.sample(config.N,
                  config.n,
                  [&](ULONG elem) {
                      // fprintf(fp, "%lld\n", elem);
                      sample.push_back(elem);
                      // samples_taken++;
                  });
    }

    std::cout << "measurements" << std::endl;
    for (ULONG iteration = 0; iteration < config.iterations; ++iteration) {
        sample.clear();
        // samples_taken = 0;
        t.reset();

        // Compute sample
        HashSampling<> hs(config.seed + iteration, config.n);
        hs.sample(config.N,
                  config.n,
                  [&](ULONG elem) {
                      // fprintf(fp, "%lld\n", elem);
                      sample.push_back(elem);
                      // samples_taken++;
                  });

        stats.push(t.get());
    }

    std::cout << "RESULT runner=H"
              << " time=" << stats.avg()
              << " stddev=" << stats.stddev()
              << " iterations=" << config.iterations
              << " gen=" << rng::select_t::name << std::endl;
    fprintf(fp, "RESULT runner=H time=%f stddev=%f iterations=%llu\n",
            stats.avg(), stats.stddev(), config.iterations);
    fclose(fp);
}
