cmake_minimum_required(VERSION 3.16)

project(cfd-benchmarks CXX)

if (NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    message(FATAL_ERROR "Do not support non-Linux now!")
endif()




add_compile_options("-g")

add_subdirectory(src)

set_property(
	GLOBAL
    
	PROPERTY
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)