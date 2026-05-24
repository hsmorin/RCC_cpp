cmake_minimum_required(VERSION 4.3.3)
project(cpp)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ortools REQUIRED)
set(CMAKE_PREFIX_PATH 
    "/home/linuxbrew/.linuxbrew"
    "/home/linuxbrew/.linuxbrew/opt/clp"
    "/home/linuxbrew/.linuxbrew/opt/cbc"
    "/home/linuxbrew/.linuxbrew/opt/glpk"
)
add_executable(cpp src/findIndicatorMats.cpp)
target_link_libraries(cpp ortools::ortools)
