project(ROE CXX)

add_executable(_roe_sw_swp _roe_sw_swp.cpp)

add_executable(eno eno.cpp)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/kernels)