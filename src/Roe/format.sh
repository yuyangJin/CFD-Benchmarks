#!/bin/bash
set -x
SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd) 
find ${SCRIPT_DIR}/include -iname '*.h' -o -iname '*.cpp' | xargs clang-format -i -style=file
find ${SCRIPT_DIR}/kernels -iname '*.h' -o -iname '*.cpp' | xargs clang-format -i -style=file
clang-format ${SCRIPT_DIR}/eno.cpp -i -style=file
clang-format ${SCRIPT_DIR}/muscl.cpp -i -style=file
clang-format ${SCRIPT_DIR}/roe.cpp -i -style=file
clang-format ${SCRIPT_DIR}/upwindtvd.cpp -i -style=file
clang-format ${SCRIPT_DIR}/symtvd.cpp -i -style=file
clang-format ${SCRIPT_DIR}/nnd.cpp -i -style=file
clang-format ${SCRIPT_DIR}/weno.cpp -i -style=file