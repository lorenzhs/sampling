set(stocc_sources 
    "stoc1.cpp"
    "stoc2.cpp"
    "userintf.cpp"
)

set(stocc_headers 
    "stocc.h"
    "randomc.h"
)

add_library(stocc ${stocc_sources} ${stocc_headers})

target_include_directories(stocc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(stocc)
