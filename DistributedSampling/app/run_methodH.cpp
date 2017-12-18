/******************************************************************************
 * run_methodH.cpp 
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

#include <vector>
#include <algorithm>

#include "timer.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "sampling_config.h"
#include "sampling/methodH.h"
#include "tools/benchmark.h"
#include "tools/arg_parser.h"

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
        t.restart();

        // Compute sample
        HashSampling<> hs(config.seed + iteration, config.n);
        hs.sample(config.N,
                  config.n,
                  [&](ULONG elem) {
                      // fprintf(fp, "%lld\n", elem);
                      sample.push_back(elem);
                      // samples_taken++;
                  });

        double time = t.elapsed();
        stats.push(time);    
    }

    std::cout << "RESULT runner=H" 
              << " time=" << stats.avg()
              << " stddev=" << stats.stddev()
              << " iterations=" << config.iterations << std::endl;
    fprintf(fp, "RESULT runner=H time=%f stddev=%f iterations=%llu\n", stats.avg(), stats.stddev(), config.iterations);
    fclose(fp);
}

