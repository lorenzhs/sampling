set(city_sources 
    "city.cc"
)

set(city_headers 
    "city.h"
)

add_library(city ${city_sources} ${city_headers})

target_include_directories(city PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

