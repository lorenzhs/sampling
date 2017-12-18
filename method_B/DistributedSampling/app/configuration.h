/******************************************************************************
 * configuration.h 
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


#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include "sampling_config.h"

class configuration {
        public:
                configuration() {} ;
                virtual ~configuration() {};

                void standard(SamplingConfig & config);
};

inline void configuration::standard(SamplingConfig & config) {
    config.seed                                   = 0;
    config.n                                      = 0;
    config.N                                      = 0;
    config.k                                      = 1;
    config.p                                      = 1.0;
    config.output_file                            = "tmp";
    config.write_to_disk                          = false;
    config.iterations                             = 3;
}

#endif 

