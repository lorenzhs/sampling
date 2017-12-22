################################################################################
# tests/CMakeLists.txt
#
# Copyright (C) 2015-2017 Timo Bingmann <tb@panthema.net>
# Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

include_directories(SYSTEM ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# macro for building test target programs with correct libraries
macro(sampling_build_target TARGETNAME)

  add_executable(${TARGETNAME} ${ARGN})
  target_link_libraries(${TARGETNAME} sampling gtest gtest_main)

endmacro(sampling_build_target)

# macro for building test programs, without main() in gtest_main
macro(sampling_build_plain PROGNAME)

  string(REPLACE "/" "_" TESTNAME "${PROGNAME}") # replace slashes

  sampling_build_target(${TESTNAME} ${PROGNAME}.cpp ${ARGN})

endmacro(sampling_build_plain)

# macro for building test programs, adding gtest runner in gtest_main
macro(sampling_build_only PROGNAME)

  # append gtest runner program.
  sampling_build_plain(${PROGNAME} ${ARGN})

endmacro(sampling_build_only)

# macro for registering test programs: maybe prepend valgrind
macro(sampling_test_only TESTNAME)

  set(TARGETNAME ${TESTNAME} ${ARGN})
  string(REPLACE "/" "_" TARGETNAME "${TARGETNAME}") # replace slashes
  string(REPLACE ";" "_" TARGETNAME "${TARGETNAME}") # stringify list

  if(USE_VALGRIND)
    # prepend valgrind call
    add_test(
      NAME ${TARGETNAME}
      COMMAND /usr/bin/valgrind ${SAMPLING_VALGRIND_OPTS}
      --xml=yes --xml-file=${TESTNAME}.xml
      ./${TESTNAME} ${ARGN})
  else()
    add_test(
      NAME ${TARGETNAME}
      COMMAND ${TESTNAME} ${ARGN})
  endif()

endmacro(sampling_test_only)

# macro for building and running test programs
macro(sampling_build_test PROGNAME)

  sampling_build_only(${PROGNAME})

  string(REPLACE "/" "_" TESTNAME "${PROGNAME}") # replace slashes
  sampling_test_only(${TESTNAME})

endmacro(sampling_build_test)

### list of tests in subdirectories

sampling_build_test(hash_test)
sampling_build_test(hypergeometric_distribution_test)
sampling_build_test(statistics_test)

if(MPI_FOUND)
#   sampling_build_only(mpi_test)
#   # run test with mpirun
#   add_test(mpi_test1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 mpi_test)
#   add_test(mpi_test2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 mpi_test)
#   add_test(mpi_test3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 mpi_test)
#   add_test(mpi_test7 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 7 mpi_test)
#   add_test(mpi_test8 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 mpi_test)
endif()

################################################################################
