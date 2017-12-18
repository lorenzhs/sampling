/******************************************************************************
 * run_methodP.cpp 
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

#include <argtable2.h>
#include <mpi.h>

#include <vector>

#include "timer.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "sampling_config.h"
#include "sampling/methodP.h"
#include "tools/benchmark.h"

int main(int argn, char **argv) {
    // Init MPI
    MPI_Init(&argn, &argv);    
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Read command-line args
    SamplingConfig config;
    int ret_code = parse_parameters(argn, argv, 
                                    config); 
    if (ret_code) { MPI_Finalize(); return 0; }

    // Main algorithm
    FILE *fp;
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
    sample.reserve(config.n);

    // Timers
    timer t;
    statistics stats;

    std::cout << "warmup" << std::endl;
    for (ULONG iteration = 0; iteration < std::min((ULONG)100, config.iterations); ++iteration) {
        sample.clear();
        MPI_Barrier(MPI_COMM_WORLD);

        // Compute sample
        ParDivideSampling<> pds(config, config.seed + iteration, size);
        pds.sample(config.n,
                   1,
                   size,
                   rank,
                   [&](ULONG elem) {
                       // fprintf(fp, "%lld\n", elem);
                       sample.push_back(elem);
                   });
    }

    std::cout << "measurements" << std::endl;
    for (ULONG iteration = 0; iteration < config.iterations; ++iteration) {
        sample.clear();
        MPI_Barrier(MPI_COMM_WORLD);
        t.restart();

        // Compute sample
        ParDivideSampling<> pds(config, config.seed + iteration, size);
        pds.sample(config.n,
                   1,
                   size,
                   rank,
                   [&](ULONG elem) {
                       // fprintf(fp, "%lld\n", elem);
                       sample.push_back(elem);
                   });

        double time = t.elapsed();
        double max_time = 0.0;
        MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == ROOT) stats.push(max_time);
    }


    if (rank == ROOT) {
        std::cout << "RESULT runner=P" 
                  << " time=" << stats.avg()
                  << " stddev=" << stats.stddev()
                  << " iterations=" << config.iterations << std::endl;
        fprintf(fp, "RESULT runner=P time=%f stddev=%f iterations=%llu\n", stats.avg(), stats.stddev(), config.iterations);
        fclose(fp);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
