/******************************************************************************
 * parse_parameters.h 
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


#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include <regex.h>
#include <string.h>
#include "configuration.h"

int parse_parameters(int argn, char **argv, 
                     SamplingConfig & config) {

    const char *progname = argv[0];

    // Setup argtable parameters.
    struct arg_lit *help                           = arg_lit0(NULL, "help", "Print help.");
    struct arg_int *user_seed                      = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
    struct arg_lit *ex_n                           = arg_lit0(NULL, "exact_n", "Use exact number of samples.");
    struct arg_lit *ex_N                           = arg_lit0(NULL, "exact_N", "Use exact size of population.");
    struct arg_int *n                              = arg_int0(NULL, "n", NULL, "Number of samples to the power of 2");
    struct arg_int *N                              = arg_int0(NULL, "N", NULL, "Size of population to the power of 2");
    struct arg_int *k                              = arg_int0(NULL, "k", NULL, "Base case size to the power of 2");
    struct arg_dbl *p                              = arg_dbl0(NULL, "p", NULL, "Sampling probability");
    struct arg_int *i                              = arg_int0(NULL, "i", NULL, "Iterations");
    struct arg_str *output                         = arg_str0(NULL, "output", NULL, "Output filename");
    struct arg_end *end                            = arg_end(100);

    // Define argtable.
    void* argtable[] = {
        help, user_seed, ex_n, ex_N, n, N, k, i, output, 
        end
    };

    // Parse arguments.
    int nerrors = arg_parse(argn, argv, argtable);

    // Catch case that help was requested.
    if (help->count > 0) {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        arg_print_glossary(stdout, argtable,"  %-40s %s\n");
        printf("This is the experimental graph generator program.\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

        return 1;
    }


    if (nerrors > 0) {
        arg_print_errors(stderr, end, progname);
        printf("Try '%s --help' for more information.\n",progname);
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1; 
    }

    configuration cfg;
    cfg.standard(config);

    if (user_seed->count > 0) {
        config.seed = user_seed->ival[0];
    }

    if (n->count > 0) {
        if (ex_n->count > 0)
            config.n = n->ival[0];
        else
            // config.n = pow(10,n->ival[0]);
            config.n = pow(2,n->ival[0]);
    }

    if (N->count > 0) {
        if (ex_N->count > 0)
            config.N = N->ival[0];
        else 
            config.N = pow(2,N->ival[0]);
    }
    
    if (k->count > 0) {
        config.k = pow(2, k->ival[0]);
    }

    if (i->count > 0) {
        config.iterations = i->ival[0];
    }

    if (p->count > 0) {
        config.p = p->dval[0];
    }

    if (output->count > 0) {
        config.output_file = output->sval[0];
        config.write_to_disk = true;
    }

    return 0;
}

#endif 

