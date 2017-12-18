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

#include <algorithm>
#include <argtable2.h>
#include <mpi.h>

#include "timer.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "sampling_config.h"
#include "sampling/methodD.h"
#include "sampling/methodH.h"
#include "sampling/methodR.h"
#include "sampling/methodP.h"

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
    if (rank == ROOT) std::cout << "sample (n=" << config.n << ", N=" << config.N << ", k=" << config.k << ", s=" << config.seed << ", p=" << size << ")\n\n";
    
    // Resulting samples
    std::vector<ULONG> sample;
    sample.reserve(config.n);

    std::vector<std::pair<std::string, double>> results(4);
    results[0] = std::make_pair("vitter", 0.0);
    results[1] = std::make_pair("hashing", 0.0);
    results[2] = std::make_pair("sequential divide-and-conquer", 0.0);
    results[3] = std::make_pair("parallel divide-and-conquer", 0.0);

    // Timers
    timer t;
    t.restart();
    double time_taken;

    if (rank == ROOT) std::cout << "perform measurements over " << config.iterations << " iterations\n";
    if (rank == ROOT) std::cout << "===============================================\n";
    for (ULONG i = 0; i < config.iterations; ++i) {
        // Sampling methods
        Vitter<> vi(config.seed);

        HashSampling<> hs(config.seed);
        hs.resizeTable(config.N, config.n);
        
        HashSampling<> shs(config.seed);
        shs.resizeTable(config.N, config.k);
        SeqDivideSampling<> sds(shs, config.k, config.seed);

        ParDivideSampling<> pds(config, config.seed, size);

        // Measurements
        if (rank == ROOT) std::cout << "vitter\n";
        if (rank == ROOT) std::cout << "warmup... ";
        vi.sample(config.N, config.n, [&](ULONG elem) { sample.push_back(elem); });
        if (rank == ROOT) std::cout << "done\n";
        sample.clear();

        if (rank == ROOT) std::cout << "sampling... ";
        t.restart();
        vi.sample(config.N, config.n, [&](ULONG elem) { sample.push_back(elem); });
        time_taken = t.elapsed();
        if (rank == ROOT) std::cout << "done\t\t" << "[" << time_taken << "]\n\n";
        results[0].second += time_taken;
        sample.clear();

        if (rank == ROOT) std::cout << "hashing\n";
        if (rank == ROOT) std::cout << "warmup... ";
        hs.sample(config.N, config.n, [&](ULONG elem) { sample.push_back(elem); });
        if (rank == ROOT) std::cout << "done\n";
        sample.clear();

        if (rank == ROOT) std::cout << "sampling... ";
        t.restart();
        hs.sample(config.N, config.n, [&](ULONG elem) { sample.push_back(elem); });
        time_taken = t.elapsed();
        if (rank == ROOT) std::cout << "done\t\t" << "[" << time_taken << "]\n\n";
        results[1].second += time_taken;
        sample.clear();

        if (rank == ROOT) std::cout << "sequential divide-and-conquer\n";
        if (rank == ROOT) std::cout << "warmup... ";
        sds.sample(config.N, config.n, [&](ULONG elem) { sample.push_back(elem); });
        if (rank == ROOT) std::cout << "done\n";
        sample.clear();
        
        if (rank == ROOT) std::cout << "sampling... ";
        t.restart();
        sds.sample(config.N, config.n, [&](ULONG elem) { sample.push_back(elem); });
        time_taken = t.elapsed();
        if (rank == ROOT) std::cout << "done\t\t" << "[" << time_taken << "]\n\n";
        results[2].second += time_taken;
        sample.clear();

        if (rank == ROOT) std::cout << "parallel divide-and-conquer\n";
        if (rank == ROOT) std::cout << "warmup... ";
        pds.sample(config.n, 1, size, rank, [&](ULONG elem) { sample.push_back(elem); });
        if (rank == ROOT) std::cout << "done\n";
        sample.clear();

        if (rank == ROOT) std::cout << "sampling... ";
        t.restart();
        pds.sample(config.n, 1, size, rank, [&](ULONG elem) { sample.push_back(elem); });
        time_taken = t.elapsed();
        if (rank == ROOT) std::cout << "done\t\t" << "[" << time_taken << "]\n\n";
        results[3].second += time_taken;
        sample.clear();

        // Increase seed
        config.seed++;
    }

    // Sort results by average speed
    for (ULONG i = 0; i < results.size(); i++) results[i].second /= config.iterations;
    std::sort(results.begin(), results.end(), [](const std::pair<std::string, double> &left, const std::pair<std::string, double> &right) {
            return left.second < right.second; 
        });

    // Output
    if (rank == ROOT) std::cout << "ranking on average speed over " << config.iterations << " iterations\n";
    if (rank == ROOT) std::cout << "===============================================\n";
    for (ULONG i = 0; i < results.size(); ++i) {
        if (rank == ROOT) std::cout << i+1 << ". " << results[i].first << "\n";
        if (rank == ROOT) std::cout << "\t\t\t\t[" << results[i].second << "]\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
