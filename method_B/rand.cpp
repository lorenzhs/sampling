#include <algorithm>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>

#include "include/arg_parser.h"
#include "include/util.h"
#include "include/sampler.h"

#ifndef USE64BIT
#include "include/MKL_gen.h"
#endif

#ifndef NOSTD
#include "include/std_gen.h"
#endif

#ifndef VERSION
#define VERSION "n/a"
#endif

template <typename T, typename F>
void run(F&& runner, const std::vector<std::unique_ptr<T[]>> &data,
         int num_threads, int iterations, std::string name,
         const std::string extra = "", const bool verbose = false,
         const bool quiet = false) {

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    global_stats stats;
    // stats for multi-threaded version including sync / load imbalance
    statistics mt_stats;
    timer t;

    for (int iteration = 0; iteration < iterations; ++iteration) {
        threads.clear();

        t.reset();
        if (num_threads > 1) {
            for (int thread = 0; thread < num_threads; ++thread) {
                threads.emplace_back(runner, data[thread].get(), thread,
                                     iteration, &stats);
            }

            for (auto &thread : threads) {
                thread.join();
            }
        } else {
            runner(data[0].get(), 0, iteration, &stats);
        }

        double time = t.get();
        mt_stats.push(time);

        // we can't print the other info because it's encapsulated in the runner
        if (verbose)
            std::cout << "RESULT runner=" << name << " time=" << time
                      << " num_threads=" << num_threads
                      << " iteration=" << iteration
                      << extra << std::endl;
    }

    if (quiet) return;
    std::cout << "RESULT runner=" << name
              << " time=" << stats.s_sum.avg()
              << " stddev=" << stats.s_sum.stddev();
    if (num_threads > 1)
        std::cout << " mt_time=" << mt_stats.avg()
                  << " mt_dev=" << mt_stats.stddev();
    std::cout << " t_gen=" << stats.s_gen.avg()
              << " t_gen_dev=" << stats.s_gen.stddev()
              << " t_fix=" << stats.s_fix.avg()
              << " t_fix_dev=" << stats.s_fix.stddev()
              << " t_pref=" << stats.s_prefsum.avg()
              << " t_pref_dev=" << stats.s_prefsum.stddev()
              << " numthreads=" << num_threads
              << " iterations=" << iterations
#ifdef FIX_STABLE
              << " fixer=stable"
#else
              << " fixer=fast"
#endif
              << extra << std::endl;
}

#ifndef USE64BIT
using T = int32_t;
#else
using T = int64_t;
#endif
static std::mutex cout_mutex;

int main(int argc, char** argv) {
    arg_parser args(argc, argv);
    size_t universe = args.get<size_t>("n", 1<<30);
    size_t k = args.get<size_t>("k", 1<<20); // sample size

    int num_threads = args.get<int>("t", 1);
    int iterations = args.get<int>("i", 1);
    const bool verbose = args.is_set("v") || args.is_set("vv");
    const bool very_verbose = args.is_set("vv");
    const bool quiet = args.is_set("q");

    double p; size_t ssize;
    std::tie(p, ssize) = sampler::calc_params(universe, k);

    std::cout << "Geometric sampler, " << k << " samples per thread "
              << "(p = " << p << ") from universe of size " << universe
              << ", using " << num_threads << " thread(s), "
              << iterations << " iteration(s), git: " << VERSION << std::endl;

    auto data = std::vector<std::unique_ptr<T[]>>(num_threads);
    // initialize in parallel
    run([ssize, &data]
        (auto /* dataptr */, int thread, int /* iteration */, auto /*stats*/) {
            data[thread] = std::make_unique<T[]>(ssize); // weak scaling
            // ensure that the memory is initialized
            std::fill(data[thread].get(), data[thread].get() + ssize, 0);
        }, data, num_threads, 1, "init", "", false, quiet);


    // Define the lambdas here because of an ICC bug that causes it to fail to
    // compile with nested lamdbas.  Ugly preprocessor hackery for feature
    // detection.
#ifndef USE64BIT
    auto mkl_gen = [](auto begin, auto end, double p, unsigned int seed)
        { MKL_gen::generate_block(begin, end-begin, p, seed); };
#else
    auto mkl_gen = [](auto, auto, double, unsigned int) {};
#endif
#ifndef NOSTD
    auto std_gen = [](auto begin, auto end, double p, unsigned int seed)
        { std_gen::generate_block(begin, end, p, seed); };
#else
    auto std_gen = [](auto, auto, double, unsigned int) {};
#endif
    auto prefsum = [](auto begin, auto end)
        { sampler::inplace_prefix_sum_disp<true>(begin, end); };

    // warmup
    size_t k_warmup = std::min<size_t>(1<<16, k);
    std::cout << "Running warmup (" << k_warmup << " samples)" << std::endl;
    run([k_warmup, universe, mkl_gen, std_gen, prefsum]
        (auto data, int /*thread*/, int /*iteration*/, auto stats) {
            double p_warmup; size_t ssize_warmup;
            std::tie(p_warmup, ssize_warmup) =
                sampler::calc_params(universe, k_warmup);
#ifndef USE64BIT
            // MKL_gen
            sampler::sample(data, ssize_warmup, k_warmup, p_warmup, universe,
                            mkl_gen, prefsum, stats);
#endif

#ifndef NOSTD
            // std_gen
            sampler::sample(data, ssize_warmup, k_warmup, p_warmup, universe,
                            std_gen, prefsum, stats);
#endif
        }, data, num_threads, 100, "warmup", "", verbose, quiet);

    std::stringstream extra_stream;
    extra_stream << " k=" << k << " b=" << ssize
                 << " p=" << p << " N=" << universe;
    auto extra = extra_stream.str();

    std::cout << "Running measurements..." << std::endl;

#ifndef USE64BIT
    // Measure MKL_gen
    run([mkl_gen, prefsum, universe, k, p, ssize, num_threads,
         verbose, very_verbose]
        (auto data, int thread_id, int iteration, auto stats){
            auto msg = sampler::sample(data, ssize, k, p, universe,
                                       mkl_gen, prefsum, stats, very_verbose);

            if (verbose) {
                cout_mutex.lock();
                // timer output is in milliseconds
                std::cout << msg << " method=mkl"
                          << " thread_id=" << thread_id
                          << " num_threads=" << num_threads
                          << " iteration=" << iteration << std::endl;
                cout_mutex.unlock();
            }
        }, data, num_threads, iterations, "mkl", extra);
#endif

#ifndef NOSTD
    // Measure std_gen
    run([std_gen, prefsum, universe, k, p, ssize, num_threads,
         verbose, very_verbose]
        (auto data, int thread_id, int iteration, auto stats){
            auto msg = sampler::sample(data, ssize, k, p, universe,
                                       std_gen, prefsum, stats, very_verbose);

            if (verbose) {
                cout_mutex.lock();
                // timer output is in milliseconds
                std::cout << msg << " method=std"
                          << " thread_id=" << thread_id
                          << " num_threads=" << num_threads
                          << " iteration=" << iteration << std::endl;
                cout_mutex.unlock();
            }
        }, data, num_threads, iterations, "std", extra);
#endif
}
