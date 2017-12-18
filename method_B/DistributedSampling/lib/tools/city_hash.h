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


#ifndef _CITY_HASH_H_
#define _CITY_HASH_H_

#include "city.h"
#include "definitions.h"

#if __GNUC__
  #if __x86_64__ || __ppc64__
    #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

#define SEEDA 28475421
#define SEEDB 52150599

class CityHash {
    public:
        static ULONG hash(ULONG x) {
#ifdef ENV64BIT 
            ULONG hash = CityHash64WithSeed(&x, 8, SEEDA);
#else
            ULONG hash = CityHash32(&x, 8);
#endif
            return hash;
        };

};

#endif 

