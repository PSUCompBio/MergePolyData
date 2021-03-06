# MergePolyData
Merge multiple VTK files into a single VTK file with multiple parts

# ImageCapture
Capture .ply mesh file with texture as png image

# Quad2Hex
Convert .ply mesh file with <b>Quad</b> elements to .vtk mesh file with <b>Hexahedron</b> elements

# MultipleViewPorts
Create post-processing image from simulation data

# To compile:
0) First you need VTK. Here are instructions on that: https://psucompbio.gitbook.io/femtech/external-programs/vtk
1) cd MergePolyData
2) mkdir build
3) cd build
4) ccmake ../
5) You might need to supply cmake the path to your VTK build directory

## To use MergePolyData, at the command prompt:
MergePolyData.exe -in file1.vtk file2.vtk -out file12.vtk
### To also export as ABAQUS .inp file
MergePolyData.exe -in file1.vtk file2.vtk -out file12.vtk -abaqus

## To use ImageCapture, at the command prompt:
ImageCapture.exe FileName(.ply)  ImageName(.jpg/.png)  outputImage(.png)  magnification(optional)
### Example on Linux (with X server)
ImageCapture.exe model.ply  model.jpg  test.png

### Example on AWS (without X server)
xvfb-run ImageCapture.exe model.ply  model.jpg  test.png

## To use Quad2Hex, at the command prompt:
Quad2Hex.exe -in QuadFile(.ply) -out HexFile(.vtk)
### Example
Quad2Hex.exe -in chank.ply -out chank_Hex.vtk

## To use MultipleViewPorts, at the command prompt:
MultipleViewPorts meshFile(.ply) textureFile(.jpg/.png) datFile(.dat) outputImage(.jpg/.png) magnification(optional)
### Example
MultipleViewPorts brain3.ply Br_color3.jpg maxstrain.dat maxstrain.png
