#include "timer.h"
#include "macros_assertions.h"
#include "sampling_config.h"
#include "sampling/methodD.h"
#include "tools/benchmark.h"
#include "tools/arg_parser.h"

#include "stocc.h"
#include "var_gen.h"

int main(int argn, char **argv) {
    // Read command-line args
    SamplingConfig config;

    arg_parser args(argn, argv);
    config.N = (ULONG) 1 << args.get<ULONG>("N", 30);
    config.n = (ULONG) 1 << args.get<ULONG>("n", 20); // sample size
    config.seed = args.get<ULONG>("seed", 1);
    config.iterations = args.get<ULONG>("i", 1);
    config.output_file = args.get<std::string>("output", "tmp");

    // Main algorithm
    FILE *fp;
    std::cout << "sample (n=" << config.n << 
                          ", N=" << config.N << 
                          ", s=" << config.seed <<
                          ", i=" << config.iterations << ")" << std::endl;
    std::string filename = config.output_file;
    fp = fopen(filename.c_str(), "w+");
    
    // Statistics
    timer t;
    statistics stats;

    VarGen bino(config.seed);
    StochasticLib1 stocc(config.seed);
    std::cout << "warmup" << std::endl;
    for (ULONG iteration = 0; iteration < std::min((ULONG)100, config.iterations); ++iteration) {
        std::cout << bino.Binomial(config.n, 0.5) << std::endl;
        std::cout << stocc.Binomial(config.n, 0.5) << std::endl;
    }

    std::cout << "measurements" << std::endl;
    for (ULONG iteration = 0; iteration < config.iterations; ++iteration) {
        t.restart();
        std::cout << bino.Binomial(config.n, 0.5) << std::endl;
        std::cout << stocc.Binomial(config.n, 0.5) << std::endl;
        double time = t.elapsed();
        stats.push(time);    
    }

    std::cout << "RESULT runner=test" 
              << " time=" << stats.avg()
              << " stddev=" << stats.stddev()
              << " iterations=" << config.iterations << std::endl;
    fprintf(fp, "RESULT runner=test time=%f stddev=%f iterations=%llu\n", stats.avg(), stats.stddev(), config.iterations);
    fclose(fp);
}

