
cmake_minimum_required(VERSION 2.8)
 
project(ImageFileReader)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
add_executable(ImageFileReader easyitk.cpp)
target_link_libraries(ImageFileReader ${ITK_LIBRARIES})
