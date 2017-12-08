/*******************************************************************************
 * benchmark/configuration.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef CONFIGURATION_HEADER
#define CONFIGURATION_HEADER

#include <include/sampling_config.hpp>

namespace sampling {

class configuration {
public:
    configuration() {} ;
    virtual ~configuration() {};

    void standard(SamplingConfig & config);
};

inline void configuration::standard(SamplingConfig & config) {
    config.seed          = 0;
    config.n             = 0;
    config.N             = 0;
    config.k             = 1;
    config.p             = 1.0;
    config.output_file   = "tmp";
    config.write_to_disk = false;
    config.iterations    = 3;
}

}

#endif // CONFIGURATION_HEADER
