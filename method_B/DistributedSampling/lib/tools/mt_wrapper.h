/******************************************************************************
 * mt_wrapper.h
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


#ifndef _MT_WRAPPER_H_
#define _MT_WRAPPER_H_

#include "definitions.h"

extern "C" {
    #include "mt64.h"
    void init_genrand64(unsigned long long seed);
    void init_by_array64(unsigned long long init_key[], 
                 unsigned long long key_length);
    unsigned long long genrand64_int64(void);
    double genrand64_real2(void);
}

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(const char *ErrorText);// System-specific error reporting (userintf.cpp)

class MTWrapper {
    public:
        MTWrapper() {};

        MTWrapper(ULONG seed) {
            RandomInit(seed);
        };

        void RandomInit(ULONG seed) {
            init_genrand64(seed);
        }

        void RandomInitByArray(ULONG seeds[], ULONG NumSeeds) {
            init_by_array64(seeds, NumSeeds);
        }

        ULONG BRandom() {
            return genrand64_int64();
        }

        double Random() {
            return genrand64_real2();
        }

        ULONG IRandom(ULONG min, ULONG max) {
            if (max == min) return min;
            ULONG r = ULONG((double)(ULONG)(max - min + 1) * Random() + min); 
            if (r > max) r = max;
            return r;
        }
};

#endif
