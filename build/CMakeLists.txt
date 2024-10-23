cmake_minimum_required(VERSION 3.5)
project(wreath-math-mpz)

include_directories(
    ${PROJECT_SOURCE_DIR}/../include
    ${PROJECT_SOURCE_DIR}/../src
)

file(GLOB all_srcs
    "${PROJECT_SOURCE_DIR}/../include/*.h"
    "${PROJECT_SOURCE_DIR}/../src/*.c"
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(main ${all_srcs})
