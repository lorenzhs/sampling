set(dsfmt_sources 
    "dSFMT.c"
)

set(dsfmt_headers 
    "dSFMT.h"
)

add_library(dsfmt ${dsfmt_sources} ${dsfmt_headers})

target_include_directories(dsfmt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

