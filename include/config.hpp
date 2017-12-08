/*******************************************************************************
 * include/config.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef CONFIG_HEADER
#define CONFIG_HEADER

namespace sampling {

// MSVC doesn't define __SSE4_2__, so also check for __AVX__ // NOLINT
#if defined(__SSE4_2__) || defined(__AVX__)
#define SAMPLING_HAVE_SSE4_2
#endif

// 64 or 32 bit environment?
#if __GNUC__
  #if __x86_64__ || __ppc64__
    #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif


} // namespace sampling

#endif // CONFIG_HEADER
