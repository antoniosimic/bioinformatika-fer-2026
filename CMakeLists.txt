cmake_minimum_required(VERSION 3.14)
project(BambooFilter CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${PROJECT_SOURCE_DIR}/src)

add_library(filter_lib
    src/murmurhash3.cpp
    src/hash.cpp
    src/bloom_filter.cpp
    src/cuckoo_filter.cpp
)

add_executable(test_bloom tests/test_bloom.cpp)
target_link_libraries(test_bloom filter_lib)

add_executable(test_cuckoo tests/test_cuckoo.cpp)
target_link_libraries(test_cuckoo filter_lib)