# MergePolyData
Merge multiple VTK files into a single VTK file with multiple parts

## To compile:
1) cd MergePolyData
2) mkdir build 
3) cd build
4) ccmake ../
5) You might need to supply cmake the path to your VTK build directory

## To use, at the command prompt:
MergePolyData.exe -in file1.vtk file2.vtk -out file12.vtk
