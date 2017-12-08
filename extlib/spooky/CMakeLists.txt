set(spooky_sources 
    "spooky.cpp"
)

set(spooky_headers 
    "spooky.h"
)

add_library(spooky ${spooky_sources} ${spooky_headers})

target_include_directories(spooky PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(spooky tlx)
