cmake_minimum_required(VERSION 2.8.9)
project(sampling)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)

if(NOT CMAKE_BUILD_TYPE) 
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

find_package(MPI REQUIRED)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -Wall -march=native -O3 -g -std=c++11 -DHAVE_SSE2 -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -Wall -march=native -O3 -g -std=c++11 -mavx2 -DHAVE_SSE2 -ffast-math")

include_directories(
    app
    lib
    lib/sampling
    lib/tools
)

include_directories(${MPI_INCLUDE_PATH})
 
add_subdirectory(extern)

add_executable(run_methodD "app/run_methodD.cpp")
add_executable(run_methodH "app/run_methodH.cpp")
add_executable(run_methodR "app/run_methodR.cpp")
add_executable(run_methodSH "app/run_methodSH.cpp")
add_executable(run_methodSR "app/run_methodSR.cpp")
add_executable(run_methodP "app/run_methodP.cpp")

# add_executable(test_hyper "app/test_hyper.cpp")
# add_executable(test_binom "app/test_binom.cpp")

set(LIBS 
    stocc
    dsfmt
    spooky
    city
    mpfr
    gmp
    ${MPI_LIBRARIES}
)
 
target_link_libraries(run_methodD LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodH LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodR LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodSH LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodSR LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodP LINK_PUBLIC ${LIBS})

# target_link_libraries(test_hyper LINK_PUBLIC ${LIBS})
# target_link_libraries(test_binom LINK_PUBLIC ${LIBS})
