#ifndef _VAR_GEN_H_
#define _VAR_GEN_H_

#include <cstddef> 

#include <mpfr.h>
#include <gmp.h>

#include "mpreal.h"
#include "definitions.h"

class VarGen {
    public:
        VarGen(ULONG seed) {
            mpfr::mpreal::set_default_prec(128);

            bino_n_last = -1.0;
            bino_p_last = -1.0;

            bino_h = -1.0;
            bino_a = -1.0;
            bino_bound = -1.0;

            bino_mode = -1.0;
            bino_r1 = -1.0;
            bino_g = -1.0;

            hyp_n_last = -1.0;
            hyp_m_last = -1.0;
            hyp_N_last = -1.0;

            hyp_h = -1.0;
            hyp_a = -1.0;
            hyp_fm = -1.0;
            hyp_bound = -1.0;

            gmp_randinit_mt(rng);
            gmp_randseed_ui(rng, seed);
        }

        void RandomInit(ULONG seed) {
            gmp_randseed_ui(rng, seed);
        }

        ULONG Binomial(ULONG n, double p) {
            ULONG inv = 0;                    // invert
            ULONG x;                          // result

            if (p > 0.5) {                      // faster calculation by inversion
                p = 1. - p;  inv = 1;
            }

            mpfr::mpreal nf = n;
            mpfr::mpreal pf = p;
            // ratio of uniforms method
            mpfr::mpreal bin;
            BinomialRatioOfUniforms(nf, pf, bin);
            x = (ULONG)bin;
            if (inv) {
                x = n - x;      // undo inversion
            }
            return x;
        }

        ULONG Hypergeometric(ULONG n, ULONG m, ULONG N) {
            ULONG fak, addd;
            ULONG x;

            fak = 1;  addd = 0;
            if (m > N/2) {
               // invert m
               m = N - m;
               fak = -1;  addd = n;
            }    
            if (n > N/2) {
               // invert n
               n = N - n;
               addd += fak * m;  fak = - fak;
            }    
            if (n > m) {
               // swap n and m
               x = n;  n = m;  m = x;
            }    

            mpfr::mpreal nf = n;
            mpfr::mpreal mf = m;
            mpfr::mpreal Nf = N;

            mpfr::mpreal hyp;
            HypRatioOfUnifoms(nf, mf, Nf, hyp);
            x = (ULONG)hyp;

            return x * fak + addd;
        }

    private:
        mpfr::mpreal bino_n_last, bino_p_last;
        mpfr::mpreal bino_h, bino_a, bino_bound;
        mpfr::mpreal bino_mode, bino_r1, bino_g;

        mpfr::mpreal hyp_n_last, hyp_m_last, hyp_N_last;
        mpfr::mpreal hyp_h, hyp_a, hyp_fm, hyp_bound;

        gmp_randstate_t rng;

        const double SHAT1 = 2.943035529371538573;    // 8/e
        const double SHAT2 = 0.8989161620588987408;   // 3-sqrt(12/e)

        mpfr::mpreal LnFac(mpfr::mpreal n) {
           // log factorial function. gives natural logarithm of n!

           // define constants
           static const double                 // coefficients in Stirling approximation     
              C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
              C1 =  1./12., 
              C3 = -1./360.;
           // C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
           // C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
           // static variables
           static const ULONG FAK_LEN = 1024;
           static mpfr::mpreal fac_table[FAK_LEN];   // table of ln(n!):
           static bool initialized = false;         // remember if fac_table has been initialized

           if ((ULONG)n < FAK_LEN) {
              if (n <= 1) {
                 return 0;
              }
              if (!initialized) {              // first time. Must initialize table
                 // make table of ln(n!)
                  mpfr::mpreal sum = fac_table[0] = 0.;
                 for (ULONG i = 1; i < FAK_LEN; i++) {
                    sum += log(mpfr::mpreal(i));
                    fac_table[i] = sum;
                 }
                 initialized = 1;
              }
              return fac_table[(ULONG)n];
           }
           // not found in table. use Stirling approximation
           // float128  n1, r;
           mpfr::mpreal n1 = n;  
           mpfr::mpreal r  = 1. / n1;
           return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);
        }


        void fc_lnpk(mpfr::mpreal & k, mpfr::mpreal & L, mpfr::mpreal & m, mpfr::mpreal & n, mpfr::mpreal & result) {
            result = LnFac(k) + LnFac(m - k) + LnFac(n - k) + LnFac(L + k);
        }

        void BinomialRatioOfUniforms(mpfr::mpreal & n, mpfr::mpreal & p, mpfr::mpreal & result) {
            /* 
            Subfunction for Binomial distribution. Assumes p < 0.5.

            Uses the Ratio-of-Uniforms rejection method.

            The computation time hardly depends on the parameters, except that it matters
            a lot whether parameters are within the range where the LnFac function is 
            tabulated.

            Reference: E. Stadlober: "The ratio of uniforms approach for generating
            discrete random variates". Journal of Computational and Applied Mathematics,
            vol. 31, no. 1, 1990, pp. 181-189.
            */   
            mpfr::mpreal u;                           // uniform random
            mpfr::mpreal q1;                          // 1-p
            mpfr::mpreal np;                          // n*p
            mpfr::mpreal var;                         // variance
            mpfr::mpreal lf;                          // ln(f(x))
            mpfr::mpreal x;                           // real sample
            mpfr::mpreal k;                           // integer sample

            if (bino_n_last != n || bino_p_last != p) {    // Set_up
                bino_n_last = n;
                bino_p_last = p;
                q1 = 1.0 - p;
                np = n * p;
                bino_mode = floor(np + p);             // mode
                bino_a = np + 0.5;                         // hat center
                bino_r1 = log(p / q1);
                bino_g = LnFac(bino_mode) + LnFac(n-bino_mode);
                var = np * q1;                             // variance
                bino_h = sqrt(SHAT1 * (var+0.5)) + SHAT2;  // hat width
                bino_bound = floor(bino_a + 6.0 * bino_h);// safety-bound
                if (bino_bound > n) bino_bound = n;        // safety-bound
            }

            while (1) {                                   // rejection loop
                mpfr_urandomb(u.mpfr_ptr(), rng);
                if (u == 0) continue;                      // avoid division by 0
                mpfr_urandomb(x.mpfr_ptr(), rng);
                x = bino_a + bino_h * (x - 0.5) / u;
                if (x < 0.) continue;    // reject, avoid overflow
                k = floor(x);                            // truncate
                lf = (k-bino_mode)*bino_r1+bino_g-LnFac(k)-LnFac(n-k);// ln(f(k))
                if (u * (4.0 - u) - 3.0 <= lf) break;      // lower squeeze accept
                if (u * (u - lf) > 1.0) continue;          // upper squeeze reject
                if (2.0 * log(u) <= lf) break;             // final acceptance
            }

            // return k;
            result = k;
        }

        void HypRatioOfUnifoms (mpfr::mpreal & n, mpfr::mpreal & m, mpfr::mpreal & N, mpfr::mpreal & result) {
            /*
            Subfunction for Hypergeometric distribution using the ratio-of-uniforms
            rejection method.

            This code is valid for 0 < n <= m <= N/2.

            The computation time hardly depends on the parameters, except that it matters
            a lot whether parameters are within the range where the LnFac function is
            tabulated.

            Reference: E. Stadlober: "The ratio of uniforms approach for generating
            discrete random variates". Journal of Computational and Applied Mathematics,
            vol. 31, no. 1, 1990, pp. 181-189.
            */
             
            // convert inputs 
            mpfr::mpreal L;                          // N-m-n
            mpfr::mpreal mode;                       // mode
            mpfr::mpreal k;                          // integer sample
            mpfr::mpreal x;                           // real sample
            mpfr::mpreal rNN;                         // 1/(N*(N+2))
            mpfr::mpreal my;                          // mean
            mpfr::mpreal var;                         // variance
            mpfr::mpreal u;                           // uniform random
            mpfr::mpreal lf;                          // ln(f(x))
            mpfr::mpreal lnfac;                          

            L = N - m - n;
            if (hyp_N_last != N || hyp_m_last != m || hyp_n_last != n) {
                hyp_N_last = N;  hyp_m_last = m;  hyp_n_last = n;       // Set-up
                rNN = 1. / (N*(N+2));                                   // make two divisions in one
                my = n * m * rNN * (N+2);                               // mean = n*m/N
                mode = floor((n+1) * (m+1) * rNN * N);                  // mode = floor((n+1)*(m+1)/(N+2))
                var = n * m * (N-m) * (N-n) / (N*N*(N-1));              // variance
                hyp_h = sqrt(SHAT1 * (var+0.5)) + SHAT2;                // hat width
                hyp_a = my + 0.5;                                       // hat center
                fc_lnpk(mode, L, m, n, lnfac);
                hyp_fm = lnfac;
                hyp_bound = floor(hyp_a + 4.0 * hyp_h);                 // safety-bound
                if (hyp_bound > n) hyp_bound = n;
            }

            while(1) {
                mpfr_urandomb(u.mpfr_ptr(), rng);
                if (u == 0) continue;                       // avoid division by 0
                mpfr_urandomb(x.mpfr_ptr(), rng);
                x = hyp_a + hyp_h * (x-0.5) / u;            // generate hat distribution
                if (x < 0.) continue;                       // reject, avoid overflow
                k = floor(x);
                if (k > hyp_bound) continue;                // reject if outside range
                // lf = hyp_fm - fc_lnpk(k,L,m,n);          // ln(f(k))
                fc_lnpk(k, L, m, n, lnfac);
                lf = hyp_fm - lnfac;            // ln(f(k))
                if (u * (4.0 - u) - 3.0 <= lf) break;       // lower squeeze accept
                if (u * (u-lf) > 1.0) continue;             // upper squeeze reject
                if (2.0 * log(u) <= lf) break;              // final acceptance
            }

            // return k;
            result = k;
        }
};

#endif 
