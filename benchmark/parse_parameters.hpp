/*******************************************************************************
 * benchmark/parse_parameters.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef PARSE_PARAMETERS_HEADER
#define PARSE_PARAMETERS_HEADER

#include "arg_parser.hpp"
#include "configuration.hpp"

#include <include/sampling_config.hpp>

#include <cmath>
#include <cstring>
#include <regex.h>

namespace sampling {

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

} // namespace sampling

#endif
