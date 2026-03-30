cmake_minimum_required(VERSION 3.20)
project(dspar_exp1 CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release build by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0 -fsanitize=address,undefined")

# Dependencies
find_package(igraph REQUIRED)   # sudo apt install libigraph-dev
find_package(OpenMP REQUIRED)   # sudo apt install libomp-dev
find_package(CURL REQUIRED)     # sudo apt install libcurl4-openssl-dev
find_package(ZLIB REQUIRED)     # for .gz decompression

# toml11: header-only single-include at extern/toml11/single_include
include_directories(extern/toml11/single_include)

# Sources shared between exp1 and tests
set(CORE_SOURCES
    src/graph.cpp
    src/dspar.cpp
    src/sampler.cpp
    src/modularity.cpp
    src/leiden.cpp
    src/experiment.cpp
)

# Sources only needed by the main binary (require CURL/ZLIB)
set(MAIN_ONLY_SOURCES
    src/config.cpp
    src/downloader.cpp
)

# Main binary
add_executable(exp1
    src/main.cpp
    ${CORE_SOURCES}
    ${MAIN_ONLY_SOURCES})

target_link_libraries(exp1
    PRIVATE
    igraph::igraph
    OpenMP::OpenMP_CXX
    CURL::libcurl
    ZLIB::ZLIB)

target_include_directories(exp1 PRIVATE src)

# Tests (Catch2)
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.0)
FetchContent_MakeAvailable(Catch2)

add_executable(tests
    tests/test_graph.cpp
    tests/test_dspar.cpp
    tests/test_sampler.cpp
    tests/test_modularity.cpp
    tests/test_leiden.cpp
    tests/test_experiment.cpp
    ${CORE_SOURCES})

target_link_libraries(tests
    PRIVATE
    Catch2::Catch2WithMain
    igraph::igraph
    OpenMP::OpenMP_CXX)

target_include_directories(tests PRIVATE src)

enable_testing()
add_test(NAME unit_tests COMMAND tests)
