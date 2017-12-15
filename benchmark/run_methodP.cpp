/*******************************************************************************
 * benchmark/run_methodP.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "arg_parser.hpp"
#include "parse_parameters.hpp"

#include <include/benchmark.hpp>
#include <include/methodP.hpp>
#include <include/sampling_config.hpp>
#include <include/timer.hpp>

#include <mpi.h>

#include <iostream>
#include <vector>

using namespace sampling;

int main(int argn, char **argv) {
    // Init MPI
    MPI_Init(&argn, &argv);
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Read command-line args
    SamplingConfig config;
    parse_parameters(argn, argv, config);

    // Main algorithm
    FILE *fp = nullptr; // set to nullptr to shut up -Wmaybe-uninitialized
    if (rank == ROOT) {
        std::cout << "sample (n=" << config.n <<
            ", N=" << config.N <<
            ", k=" << config.k <<
            ", s=" << config.seed <<
            ", i=" << config.iterations <<
            ", p=" << size << ")" << std::endl;
        std::string filename = config.output_file;
        fp = fopen(filename.c_str(), "w+");
    }

    // Resulting samples
    std::vector<ULONG> sample;
    sample.reserve((config.n / size) * 1.2);
    // ULONG samples_taken = 0;

    // Timers
    timer t;
    statistics stats;

    std::cout << "warmup" << std::endl;
    for (ULONG iteration = 0; iteration < std::min((ULONG)100, config.iterations); ++iteration) {
        sample.clear();
        // samples_taken = 0;
        MPI_Barrier(MPI_COMM_WORLD);

        // Compute sample
        ParDivideSampling<> pds(config.N, config.k, config.seed + iteration, size);
        pds.sample(config.n,
                   1,
                   size,
                   rank,
                   [&](ULONG elem) {
                       // fprintf(fp, "%lld\n", elem);
                       sample.push_back(elem);
                       // samples_taken++;
                   });
    }

    std::cout << "ss " << sample.size() << std::endl;

    std::cout << "measurements" << std::endl;
    for (ULONG iteration = 0; iteration < config.iterations; ++iteration) {
        sample.clear();
        // samples_taken = 0;
        MPI_Barrier(MPI_COMM_WORLD);
        t.reset();

        // Compute sample
        ParDivideSampling<> pds(config.N, config.k, config.seed + iteration, size);
        pds.sample(config.n,
                   1,
                   size,
                   rank,
                   [&](ULONG elem) {
                       // fprintf(fp, "%lld\n", elem);
                       sample.push_back(elem);
                       // samples_taken++;
                   });

        double time = t.get();
        double max_time = 0.0;
        MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == ROOT) stats.push(max_time);
    }

    if (rank == ROOT) {
        std::cout << "RESULT runner=P"
                  << " n=" << config.n
                  << " N=" << config.N
                  << " p=" << size
                  << " time=" << stats.avg()
                  << " stddev=" << stats.stddev()
                  << " iterations=" << config.iterations
                  << " gen=" << rng::select_t::name << std::endl;
        fprintf(fp, "RESULT runner=P n=%llu N=%llu p=%d time=%f stddev=%f iterations=%llu\n",
                config.n, config.N, size, stats.avg(), stats.stddev(), config.iterations);
        fclose(fp);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
