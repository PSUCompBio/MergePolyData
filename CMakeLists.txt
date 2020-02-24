
cmake_minimum_required(VERSION 2.8)

PROJECT(MergePolyData)

find_package(VTK REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})

#include_directories(include)
#Add external packages and headers
add_subdirectory(third-party)

#Add JsonCPP Include and Link Directories
ExternalProject_Get_Property(JsonCpp source_dir)
include_directories(${source_dir}/include)

ExternalProject_Get_Property(JsonCpp binary_dir)
link_directories(${binary_dir}/src/lib_json)


add_executable(MergePolyData MACOSX_BUNDLE include/vtkCleanUnstructuredGrid.h vtkCleanUnstructuredGrid.cxx MergePolyData.cxx )
add_executable(ImageCapture MACOSX_BUNDLE include/tinyply.h  tinyply.cxx ImageCapture.cxx)
add_executable(Quad2Hex MACOSX_BUNDLE Quad2Hex.cxx )
add_executable(MultipleViewPorts MACOSX_BUNDLE include/tinyply.h  tinyply.cxx MultipleViewPorts.cxx )
ADD_DEPENDENCIES(MultipleViewPorts JsonCpp)

file(COPY examples DESTINATION .)

target_link_libraries(Quad2Hex ${VTK_LIBRARIES})
target_link_libraries(MergePolyData ${VTK_LIBRARIES})
target_link_libraries(ImageCapture ${VTK_LIBRARIES})
target_link_libraries(MultipleViewPorts ${VTK_LIBRARIES} jsoncpp)
set(CMAKE_MSVCIDE_RUN_PATH ${VTK_RUNTIME_LIBRARY_DIRS})
