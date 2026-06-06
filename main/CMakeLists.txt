cmake_minimum_required(VERSION 4.3.3)
project(main)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_PREFIX_PATH 
    "/home/linuxbrew/.linuxbrew"
    "/home/linuxbrew/.linuxbrew/opt/clp"
    "/home/linuxbrew/.linuxbrew/opt/cbc"
    "/home/linuxbrew/.linuxbrew/opt/glpk"
)

find_package(ortools REQUIRED)
add_executable(main src/findIndicatorMats.cpp)
target_link_libraries(main ortools::ortools)

include(FetchContent)

FetchContent_Declare(
    treehh
    GIT_REPOSITORY https://github.com/kpeeters/tree.hh.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(treehh)

target_include_directories(main PRIVATE ${treehh_SOURCE_DIR}/src)
