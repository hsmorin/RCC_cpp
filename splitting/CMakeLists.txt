cmake_minimum_required(VERSION 4.3.3)
project(splitting)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_PREFIX_PATH 
    "../cpp_libraries/or-tools_x86_64_Ubuntu-24.04_cpp_v9.12.4544"
)

find_package(ortools REQUIRED)
add_executable(splitting src/splitting.cpp)
target_link_libraries(splitting ortools::ortools)

include(FetchContent)

FetchContent_Declare(
    treehh
    GIT_REPOSITORY https://github.com/kpeeters/tree.hh.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(treehh)

target_include_directories(splitting PRIVATE ${treehh_SOURCE_DIR}/src)
