cmake_minimum_required(VERSION 3.14)
project(BambooFilter CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Release-only optimization flags. -O3 is the standard high-optimisation
# level; -march=native lets the compiler use every instruction supported by
# the host CPU (popcnt, BMI, AVX, …); -flto enables link-time optimisation
# so the optimiser can inline across translation units (e.g. hot helpers in
# bamboo_filter.cpp into the benchmark loop). Binary is no longer portable
# to older CPUs — fine for local benchmarking, would need relaxing for a
# distributable build.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -march=native -flto)
    add_link_options(-flto)
endif()

include_directories(${PROJECT_SOURCE_DIR}/src)

add_library(filter_lib
    src/murmurhash3.cpp
    src/hash.cpp
    src/bamboo_filter.cpp
)

add_executable(test_bamboo tests/test_bamboo.cpp)
target_link_libraries(test_bamboo filter_lib)

add_executable(benchmark_runner benchmark/benchmark_main.cpp)
target_link_libraries(benchmark_runner filter_lib)

# Optional comparison against the original Wang et al. implementation.
# Linux + x86_64 + AVX2 only — the reference uses _mm256_* intrinsics and
# GCC-only `-fpermissive` features. Enabled with -DBUILD_REFERENCE=ON.
option(BUILD_REFERENCE "Build comparison against wanghanchengchn/bamboofilters" OFF)
if(BUILD_REFERENCE)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
        message(FATAL_ERROR
            "BUILD_REFERENCE requires Linux. The reference impl uses AVX2 "
            "intrinsics and -fpermissive GCC-isms that don't port elsewhere.")
    endif()

    include(FetchContent)
    FetchContent_Declare(
        bamboo_reference
        GIT_REPOSITORY https://github.com/wanghanchengchn/bamboofilters.git
        GIT_TAG main
    )
    # Populate without calling add_subdirectory — their CMake builds a test
    # binary that pulls in OpenSSL/mman.h, which we don't need.
    FetchContent_GetProperties(bamboo_reference)
    if(NOT bamboo_reference_POPULATED)
        FetchContent_Populate(bamboo_reference)
    endif()

    add_library(bamboo_reference_bobhash STATIC
        ${bamboo_reference_SOURCE_DIR}/src/common/BOBHash.cpp
    )
    target_include_directories(bamboo_reference_bobhash PUBLIC
        ${bamboo_reference_SOURCE_DIR}/src
    )
    target_compile_options(bamboo_reference_bobhash PRIVATE
        -fpermissive -march=core-avx2 -fno-strict-aliasing
    )
    set_target_properties(bamboo_reference_bobhash PROPERTIES
        CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON
    )

    add_executable(reference_benchmark reference/reference_main.cpp)
    target_include_directories(reference_benchmark PRIVATE
        ${bamboo_reference_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/src
    )
    target_link_libraries(reference_benchmark PRIVATE
        bamboo_reference_bobhash
    )
    target_compile_options(reference_benchmark PRIVATE
        -fpermissive -march=core-avx2 -fno-strict-aliasing
    )
    set_target_properties(reference_benchmark PROPERTIES
        CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON
    )
endif()