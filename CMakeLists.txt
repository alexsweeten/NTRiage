cmake_minimum_required(VERSION 3.16)
project(NTRIAGE VERSION 0.3.0 LANGUAGES CXX)

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

add_executable(ntriage
    src/main.cpp
    src/core.cpp
    src/deep.cpp
    src/fasta.cpp
)
target_include_directories(ntriage PRIVATE include)
target_link_libraries(ntriage PRIVATE ZLIB::ZLIB Threads::Threads)
target_compile_features(ntriage PRIVATE cxx_std_17)
target_compile_definitions(ntriage PRIVATE NTRIAGE_VERSION="${PROJECT_VERSION}")

if(MSVC)
    target_compile_options(ntriage PRIVATE /W4 /permissive-)
else()
    target_compile_options(ntriage PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
    add_test(NAME ntriage_version COMMAND ntriage --version)
    set_tests_properties(ntriage_version PROPERTIES
        PASS_REGULAR_EXPRESSION "N-TRIAGE ${PROJECT_VERSION}"
    )
    add_test(NAME ntriage_self_test COMMAND ntriage --self-test)
    add_test(
        NAME ntriage_integration_test
        COMMAND ${CMAKE_COMMAND}
            -DNTRIAGE=$<TARGET_FILE:ntriage>
            -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration.cmake
    )
endif()
