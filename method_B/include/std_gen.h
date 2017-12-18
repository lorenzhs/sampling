#pragma once

#include <cassert>
#include <limits>
#include <random>

struct std_gen {
    template <typename It>
    static void generate_block(It begin, It end, double p,
                               unsigned int seed = 0) {
        using value_type = typename std::iterator_traits<It>::value_type;
        assert(p > 0 && p < 1);

        if (seed == 0) {
            seed = std::random_device{}();
        }

        std::mt19937_64 gen(seed);
        std::geometric_distribution<value_type> dist(p);

        for (auto it = begin; it < end; ++it) {
            *it = dist(gen);
        }
    }

    template <typename It>
    static void generate_block(It dest, size_t size, double p,
                               unsigned int seed = 0) {
        generate_block(dest, dest+size, p, seed);
    }
};
