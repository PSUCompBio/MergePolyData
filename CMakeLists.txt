
cmake_minimum_required(VERSION 2.8)
 
PROJECT(MergePolyData)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(include)

add_executable(MergePolyData MACOSX_BUNDLE include/vtkCleanUnstructuredGrid.h vtkCleanUnstructuredGrid.cxx MergePolyData.cxx )
 
target_link_libraries(MergePolyData ${VTK_LIBRARIES})
set(CMAKE_MSVCIDE_RUN_PATH ${VTK_RUNTIME_LIBRARY_DIRS})
