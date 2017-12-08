################################################################################
# benchmark/CMakeLists.txt
#
# CMake file for sampling benchmarks
#
# Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

build_prog(run_methodD)
build_prog(run_methodH)
build_prog(run_methodSH)
build_mpi_prog(run_methodP)
build_prog(run_methodR)
build_prog(run_methodSR)

################################################################################
