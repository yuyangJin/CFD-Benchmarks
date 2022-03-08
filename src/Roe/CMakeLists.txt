project(ROE CXX)

add_compile_options(-g -O3)

# add_executable(_roe_sw_swp _roe_sw_swp.cpp)

add_executable(eno eno.cpp)
add_executable(muscl muscl.cpp)
add_executable(roe roe.cpp)
add_executable(upwindtvd upwindtvd.cpp)
add_executable(symtvd symtvd.cpp)
add_executable(nnd nnd.cpp)
add_executable(weno weno.cpp)
# add_executable(compact compact.cpp)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/kernels)