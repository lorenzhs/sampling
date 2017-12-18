/**************************   mersenne.cpp   **********************************
* Author:        Agner Fog
* Date created:  2001
* Last modified: 2008-11-16
* Project:       randomc.h
* Platform:      Any C++
* Description:
* Random Number generator of type 'Mersenne Twister'
*
* This random number generator is described in the article by
* M. Matsumoto & T. Nishimura, in:
* ACM Transactions on Modeling and Computer Simulation,
* vol. 8, no. 1, 1998, pp. 3-30.
* Details on the initialization scheme can be found at
* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2001-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*******************************************************************************/

#include "randomc.h"

void CRandomMersenne::Init0(int64_t seed) {
   // Seed generator
   const uint64_t factor = 1812433253UL;
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(int64_t seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int64_t i = 0; i < 37; i++) BRandom();
}


void CRandomMersenne::RandomInitByArray(int64_t const seeds[], int64_t NumSeeds) {
   // Seed by more than 32 bits
   int64_t i, j, k;

   // Initialize
   Init0(19650218);

   if (NumSeeds <= 0) return;

   // Randomize mt[] using whole seeds[] array
   i = 1;  j = 0;
   k = (MERS_N > NumSeeds ? MERS_N : NumSeeds);
   for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + (uint64_t)seeds[j] + j;
      i++; j++;
      if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
      if (j >= NumSeeds) j=0;}
   for (k = MERS_N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
      if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
   mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

   // Randomize some more
   mti = 0;
   for (int64_t i = 0; i <= MERS_N; i++) BRandom();
}


uint64_t CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint64_t y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint64_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint64_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint64_t mag01[2] = {0, MERS_A};

      int64_t kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }
   y = mt[mti++];

   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;

   return y;
}


double CRandomMersenne::Random() {
   // Output random double number in the interval 0 <= x < 1
   // Multiply by 2^(-32)
   return (double)BRandom() * (1./(65536.*65536.));
}


int64_t CRandomMersenne::IRandom(int64_t min, int64_t max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int64_t r = int64_t((double)(uint64_t)(max - min + 1) * Random() + min); 
   if (r > max) r = max;
   return r;
}


int64_t CRandomMersenne::IRandomX(int64_t min, int64_t max) { // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
#ifdef  INT64_SUPPORTED
   // 64 bit integers available. Use multiply and shift method
   uint64_t interval;                    // Length of interval
   uint64_t longran;                     // Random bits * interval
   uint64_t iran;                        // Longran / 2^32
   uint64_t remainder;                   // Longran % 2^32

   interval = uint64_t(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder >= 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      RLimit = uint64_t(((uint64_t)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64_t)BRandom() * interval;
      iran = (uint64_t)(longran >> 32);
      remainder = (uint64_t)longran;
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return (int64_t)iran + min;

#else
   // 64 bit integers not available. Use modulo method
   uint64_t interval;                    // Length of interval
   uint64_t bran;                        // Random bits
   uint64_t iran;                        // bran / interval
   uint64_t remainder;                   // bran % interval

   interval = uint64_t(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when iran = 2^32 / interval
      // We can't make 2^32 so we use 2^32-1 and correct afterwards
      RLimit = (uint64_t)0xFFFFFFFF / interval;
      if ((uint64_t)0xFFFFFFFF % interval == interval - 1) RLimit++;
   }
   do { // Rejection loop
      bran = BRandom();
      iran = bran / interval;
      remainder = bran % interval;
   } while (iran >= RLimit);
   // Convert back to signed and return result
   return (int64_t)remainder + min;

#endif
}
