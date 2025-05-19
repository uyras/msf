cmake_minimum_required(VERSION 3.13)

project(MSF)

set (CMAKE_CXX_STANDARD 17)

set (ARGUMENTUM_BUILD_STATIC_LIBS ON)

find_package(OpenMP REQUIRED)

add_subdirectory(argumentum)
include_directories("argumentum/include")

add_executable(msf main.cpp MagneticSystem.cpp)
target_link_libraries(msf OpenMP::OpenMP_CXX Argumentum::argumentum)