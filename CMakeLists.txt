cmake_minimum_required(VERSION 3.16)

project(cfd-benchmarks CXX C)

if (NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    message(FATAL_ERROR "Do not support non-Linux now!")
endif()

add_subdirectory(src)
