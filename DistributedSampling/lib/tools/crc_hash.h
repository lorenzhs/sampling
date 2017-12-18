/******************************************************************************
 * crc_hash.h
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


#ifndef _CRC_HASH_H_
#define _CRC_HASH_H_

#include <x86intrin.h>
#include "definitions.h"

#if __GNUC__
  #if __x86_64__ || __ppc64__
    #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

#define UPPER_MASK 0xffff0000 /* most significant w-r bits */
#define LOWER_MASK 0x0000ffff /* least significant r bits */

#define SEEDA 28475421
#define SEEDB 52150599

class CRCHash {
    public:
        static ULONG hash(ULONG x) {
#ifdef ENV64BIT 
            ULONG hash = _mm_crc32_u64(SEEDA, x);
            hash = hash << 32; 
            hash += _mm_crc32_u64(SEEDB, x);
#else
            ULONG hash = 0;
            hash += _mm_crc32_u32(SEEDA, x & UPPER_MASK);
            hash += _mm_crc32_u32(SEEDA, x & LOWER_MASK);
            hash = hash << 32; 
            hash += _mm_crc32_u32(SEEDB, x & UPPER_MASK);
            hash += _mm_crc32_u32(SEEDB, x & LOWER_MASK);
#endif
            return hash;
        };

};

#endif 

