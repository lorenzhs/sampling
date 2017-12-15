/*******************************************************************************
 * include/rng/mkl.hpp
 *
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef RNG_MKL_HEADER
#define RNG_MKL_HEADER

#include <tlx/logger.hpp>

#include <mkl.h>
#include <mkl_vsl.h>
#include "errcheck.inc"

#include <limits>

namespace sampling {
namespace rng {

/*!
 * MKL generator wrapper
 *
 * If you need more elements at a time than an int can hold, you need to
 * `#define MKL_INT size_t` before including this file
 */
class mkl {
public:
    //! Initialize new MKL random generator
    mkl(size_t seed) {
        vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
    }

    ~mkl() {
        vslDeleteStream(&stream);
    }

    //! Generate `size` uniform doubles from [0, 1)
    void generate_block(std::vector<double> &output, size_t size) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }

        MKL_INT count = static_cast<MKL_INT>(size);
        // VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
        int status = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream,
                                  count, output.data(), 0.0, 1.0);
        CheckVslError(status);
    }

    //! Generate `size` uniform integers from [min, max] (i.e., both inclusive).
    //! Unfortunately, MKL does not support generating 64-bit random integers.
    void generate_int_block(int min, int max, std::vector<int> &output,
                            size_t size) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }

        MKL_INT count = static_cast<MKL_INT>(size);
        // VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
        int status = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream,
                                  count, output.data(), min, max + 1);
        CheckVslError(status);
    }

    //! Generate `size` geometrically distributed integers with parameter p
    void generate_geometric_block(double p, std::vector<int> &output,
                                  size_t size) {
        check_size(size);

        MKL_INT count = static_cast<MKL_INT>(size);
        int status = viRngGeometric(VSL_RNG_METHOD_GEOMETRIC_ICDF, stream,
                                    count, output.data(), p);
        CheckVslError(status);
    }


private:
    //! Check that `size` fits into an MKL_INT
    bool check_size(size_t size) {
        if (size >= std::numeric_limits<MKL_INT>::max()) {
            LOG1
                << "Error: MKL generator block size exceeds value range of MKL_INT:"
                << size << " >= " << std::numeric_limits<MKL_INT>::max();
            return false;
        }
        return true;
    }

    //! MKL state object
    VSLStreamStatePtr stream;
};

} // namespace rng
} // namespace sampling

#endif // RNG_MKL_HEADER
