cmake_minimum_required(VERSION 2.8)
 
PROJECT(UTestVTK)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(UTestVTK MACOSX_BUNDLE UTestVTK.cxx)
add_library(lb source/membrane.cpp source/uniform.cpp)
 
if(VTK_LIBRARIES)
  target_link_libraries(UTestVTK ${VTK_LIBRARIES} lb)
else()
  target_link_libraries(UTestVTK vtkHybrid vtkWidgets lb)
endif()
