#!/bin/bash

# Add file 1 and the target file size below
repo1=MergePolyData-Abaqus
file1="build/examples/4-Simple-Parts/combined.inp"
target_file_size1=1863

# Add file 2 and the target file size below
repo2=MergePolyData-Quad2Hex
file2="build/examples/convert-Quad2Hex/chank_Hex.vtk"
target_file_size2=2147

# Add file 3 and the target file size below
repo3=MergePolyData-ImageWrite
file3="build/examples/ply-to-png/test.png"
target_file_size3=246898

# Add file 4 and the target file size below
repo4=MergePolyData-MultiViewPorts
file4="build/examples/multipleViewPorts/test.png"
target_file_size4=324329

# You shoud not have to modify below
#
myfilesize1=$(wc -c <"$file1")
echo Acutal File1 Size = "$myfilesize1"
echo Target File1 Size = "$target_file_size1"

myfilesize2=$(wc -c <"$file2")
echo Acutal File2 Size = "$myfilesize2"
echo Target File2 Size = "$target_file_size2"

myfilesize3=$(wc -c <"$file3")
echo Acutal File3 Size = "$myfilesize3"
echo Target File3 Size = "$target_file_size3"

myfilesize4=$(wc -c <"$file4")
echo Acutal File4 Size = "$myfilesize4"
echo Target File4 Size = "$target_file_size4"

if [ $myfilesize1 -ge $target_file_size1 ];then
        echo Passed!
        echo "Passed" >> ~/$repo1.PASSED
        echo "Acutal File Size = "$myfilesize1" " >> ~/$repo1.PASSED
        echo "Target File Size = "$target_file_size1" " >> ~/$repo1.PASSED
else
        echo Failed!
        echo "Failed!" >> ~/$repo1.FAILED
        echo "Acutal File1 Size = "$myfilesize1" " >> ~/$repo1.FAILED
        echo "Target File1 Size = "$target_file_size1" " >> ~/$repo1.FAILED

fi

if [ $myfilesize2 -ge $target_file_size2 ];then
        echo Passed!
        echo "Passed" >> ~/$repo2.PASSED
        echo "Acutal File2 Size = "$myfilesize2" " >> ~/$repo2.PASSED
        echo "Target File2 Size = "$target_file_size2" " >> ~/$repo2.PASSED
else
        echo Failed!
        echo "Failed!" >> ~/$repo2.FAILED
        echo "Acutal File2 Size = "$myfilesize2" " >> ~/$repo2.FAILED
        echo "Target File2 Size = "$target_file_size2" " >> ~/$repo2.FAILED

fi

if [ $myfilesize3 -ge $target_file_size3 ];then
        echo Passed!
        echo "Passed" >> ~/$repo3.PASSED
        echo "Acutal File3 Size = "$myfilesize3" " >> ~/$repo3.PASSED
        echo "Target File3 Size = "$target_file_size3" " >> ~/$repo3.PASSED
else
        echo Failed!
        echo "Failed!" >> ~/$repo3.FAILED
        echo "Acutal File3 Size = "$myfilesize3" " >> ~/$repo3.FAILED
        echo "Target File3 Size = "$target_file_size3" " >> ~/$repo3.FAILED

fi

if [ $myfilesize4 -ge $target_file_size4 ];then
        echo Passed!
        echo "Passed" >> ~/$repo4.PASSED
        echo "Acutal File4 Size = "$myfilesize4" " >> ~/$repo4.PASSED
        echo "Target File4 Size = "$target_file_size4" " >> ~/$repo4.PASSED
else
        echo Failed!
        echo "Failed!" >> ~/$repo4.FAILED
        echo "Acutal File4 Size = "$myfilesize4" " >> ~/$repo4.FAILED
        echo "Target File4 Size = "$target_file_size4" " >> ~/$repo4.FAILED

fi

