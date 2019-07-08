
cmake_minimum_required(VERSION 2.8)

PROJECT(MergePolyData)

find_package(VTK REQUIRED)

#include_directories(include)

add_executable(MergePolyData MACOSX_BUNDLE include/vtkCleanUnstructuredGrid.h vtkCleanUnstructuredGrid.cxx MergePolyData.cxx )
add_executable(ImageCapture MACOSX_BUNDLE include/tinyply.h  tinyply.cxx ImageCapture.cxx)
add_executable(Quad2Hex MACOSX_BUNDLE Quad2Hex.cxx )

file(COPY examples DESTINATION .)

target_link_libraries(Quad2Hex ${VTK_LIBRARIES})
target_link_libraries(MergePolyData ${VTK_LIBRARIES})
target_link_libraries(ImageCapture ${VTK_LIBRARIES})
set(CMAKE_MSVCIDE_RUN_PATH ${VTK_RUNTIME_LIBRARY_DIRS})
