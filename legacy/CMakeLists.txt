# Legacy compatibility layer CMakeLists.txt

file(GLOB LEGACY_SOURCES
    heat_utils_compat.cpp
    interfaces_compat.cpp
)

file(GLOB LEGACY_HEADERS
    heat_utils_compat.hpp
    interfaces_compat.hpp
)

# Legacy library (optional, for backward compatibility)
add_library(heat_equation_legacy STATIC ${LEGACY_SOURCES} ${LEGACY_HEADERS})

target_include_directories(heat_equation_legacy
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(heat_equation_legacy
    PUBLIC
        heat_equation_utils
        heat_equation_mpi
)

# Legacy main executable (for testing backward compatibility)
if(BUILD_TESTS)
    add_executable(heat_equation_legacy_test legacy_main.cpp test_legacy_compat.cpp)
    target_link_libraries(heat_equation_legacy_test
        heat_equation_legacy
        ${MPI_CXX_LIBRARIES}
    )
    add_test(NAME LegacyBackwardCompatibility COMMAND heat_equation_legacy_test)
endif()

# Installation
install(TARGETS heat_equation_legacy
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(FILES ${LEGACY_HEADERS}
        DESTINATION include/legacy)
