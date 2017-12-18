/******************************************************************************
 * sampling_config.h      
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

#ifndef _SAMPLING_CONFIG_H_
#define _SAMPLING_CONFIG_H_

#include <string>
#include "definitions.h"

// Configuration for the generator.
struct SamplingConfig {
    SamplingConfig() {}

    // Seed for the PRNG
    int seed;
    // Number of samples
    ULONG n;
    // Size of population
    ULONG N;
    // Base case size
    ULONG k;
    // Sample probability (bernoulli)
    double p;
    // Output filename
    std::string output_file;
    // Write edges directly to disk
    bool write_to_disk;
    // Number of iterations
    ULONG iterations;

    // Shared Mem OMP
    void LogDump(FILE *out) const { }
};

#endif 

