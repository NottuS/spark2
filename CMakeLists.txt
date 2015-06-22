
cmake_minimum_required(VERSION 2.8)
 
project(TrabPi)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
add_executable(detector detector.cpp)
target_link_libraries(detector ${ITK_LIBRARIES})
