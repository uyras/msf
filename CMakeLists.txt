cmake_minimum_required(VERSION 3.0)

project(MSF)

set (CMAKE_CXX_STANDARD 17)

set (ARGUMENTUM_BUILD_STATIC_LIBS ON)

find_package(OpenMP REQUIRED)

add_subdirectory(argumentum)
add_subdirectory(partsEngine)

include_directories("partsEngine" "argumentum/include")

add_executable(msf main.cpp)
target_link_libraries(msf partsEngine OpenMP::OpenMP_CXX Argumentum::argumentum)