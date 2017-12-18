#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include <string.h>
#include "sampling_config.h"
#include "tools/arg_parser.h"

void parse_parameters(int argn, char **argv, 
                      SamplingConfig & config) {

    arg_parser args(argn, argv);

    config.N = (ULONG) 1 << args.get<ULONG>("N", 30);
    config.n = (ULONG) 1 << args.get<ULONG>("n", 20); // sample size
    config.k = (ULONG) 1 << args.get<ULONG>("k", 10); // base case
    config.seed = args.get<ULONG>("seed", 1);
    config.iterations = args.get<ULONG>("i", 1);
    config.output_file = args.get<std::string>("output", "tmp");
    config.use_binom = args.is_set("binom");
}

#endif 
