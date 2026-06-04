# Utils library CMakeLists.txt

file(GLOB UTILS_SOURCES
    logger.cpp
)

file(GLOB UTILS_HEADERS
    logger.hpp
    array2d.hpp
    timer.hpp
)

add_library(heat_equation_utils ${UTILS_SOURCES} ${UTILS_HEADERS})

target_link_libraries(heat_equation_utils
    PUBLIC
        heat_equation_options
)

target_include_directories(heat_equation_utils
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Installation
install(TARGETS heat_equation_utils
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(FILES ${UTILS_HEADERS}
        DESTINATION include/utils)
