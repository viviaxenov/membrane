cmake_minimum_required(VERSION 2.8)
 
PROJECT(TestVTK)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(TestVTK MACOSX_BUNDLE test/testVTK.cpp)
add_library(lb source/membrane.cpp source/vec3.cpp)
 
if(VTK_LIBRARIES)
  target_link_libraries(TestVTK ${VTK_LIBRARIES} lb)
else()
  target_link_libraries(TestVTK vtkHybrid vtkWidgets lb)
endif()
