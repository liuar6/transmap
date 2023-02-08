cmake_minimum_required(VERSION 3.10)
project(transmap C)

set(CMAKE_C_STANDARD 99)
add_subdirectory(bioidx)
add_executable(transmap transmap.c transmap_bed.c transmap_bam.c)
target_link_libraries(transmap hts bioidx)
