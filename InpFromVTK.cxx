#include <fstream>
#include <iostream>

#include <vtkAppendFilter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>
#include <vtkIntArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

int main(int argc, char *argv[]) {
  std::string inputFile;
  std::string outFile;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-in") == 0) {
      if (argc > i + 1) {
        inputFile = argv[i + 1];
      } else {
        std::cout << "Input file name expected  after -in option" << std::endl;
      }
    } else if (strcmp(argv[i], "-out") == 0) {
      if (argc > i + 1) {
        outFile = argv[i + 1];
      } else {
        std::cout << "Output file name expected  after -out option"
                  << std::endl;
      }
    }
  }
  // Validate input file paths
  ifstream fileT(inputFile);
  if (!fileT.good()) { // file couldn't be opened
    std::cout << "Error:" << inputFile << " could not be opened" << std::endl;
    return 1;
  }

  cout << "Reading file " << inputFile << endl;
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
      vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(inputFile.c_str());
  reader->Update();
  vtkSmartPointer<vtkUnstructuredGrid> pUnstructedGrid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  pUnstructedGrid->DeepCopy(reader->GetOutput());

  std::ofstream file;
  file.open(outFile);
  file << "*************************************\n";
  file << "*HEADING\n";
  file << "FemTech\n";
  file << "*************************************\n";
  file << "*NODE, NSET=All\n";
  for (int i = 0; i < pUnstructedGrid->GetNumberOfPoints(); i++) {
    double pnt[3];
    pUnstructedGrid->GetPoint(i, pnt);
    file << i + 1 << ",  " << pnt[0] << ",  " << pnt[1] << ",  " << pnt[2]
         << "\n";
  }
  int nCells = pUnstructedGrid->GetNumberOfCells();
  std::string arrayName = "Part ID";
  vtkSmartPointer<vtkIntArray> partArray = dynamic_cast<vtkIntArray *>(
      pUnstructedGrid->GetCellData()->GetArray(arrayName.c_str()));

  if (partArray) {
    double *range = partArray->GetRange();
    int noPart = (int)range[1] + 1;
    for (int j = 0; j < noPart; j++) {
      file << "*ELEMENT,TYPE=C3D8,ELSET=PART_" << j + 1 << std::endl;
      for (int i = 0; i < nCells; i++) {
        int part = partArray->GetValue(i);
        if (part == j) {
          vtkSmartPointer<vtkIdList> cellPointIds =
              vtkSmartPointer<vtkIdList>::New();
          pUnstructedGrid->GetCellPoints(i, cellPointIds);
          // Only for Element with 8 nodes
          assert(cellPointIds->GetNumberOfIds() == 8);
          file << i + 1 << ", " << cellPointIds->GetId(0) + 1 << ", "
               << cellPointIds->GetId(1) + 1 << ", "
               << cellPointIds->GetId(2) + 1 << ", "
               << cellPointIds->GetId(3) + 1 << ", "
               << cellPointIds->GetId(4) + 1 << ", "
               << cellPointIds->GetId(5) + 1 << ", "
               << cellPointIds->GetId(6) + 1 << ", "
               << cellPointIds->GetId(7) + 1 << std::endl;
        }
      }
    }
  } else {
    std::cout << "The file does not have a CellData array named " << arrayName
              << std::endl;
    return 1;
  }
  file.close();
  std::cout << "Written Abaqus inp file : " << outFile << std::endl;
  return 0;
}
