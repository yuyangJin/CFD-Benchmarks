project(ROE CXX)

add_executable(_roe_sw_swp _roe_sw_swp.cpp)

add_executable(eno eno.cpp)
add_executable(muscl muscl.cpp)
add_executable(roe roe.cpp)
add_executable(upwindtvd upwindtvd.cpp)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/kernels)