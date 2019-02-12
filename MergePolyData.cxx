#include <iostream>
#include <fstream>

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkIntArray.h>
#include <vtkAppendFilter.h>
#include <vtkCleanUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>

#define MERGE_TOL 1E-12

bool ValidateInputFiles(std::vector<std::string>& inputFiles, std::string& outFile)
{
	//Validate input file paths
	for (int f=0; f<inputFiles.size(); f++)
	{
		ifstream file(inputFiles[f]);
		if (!file.good()) 
        { // file couldn't be opened
            std::cout << "Error:" << inputFiles[f] <<" could not be opened" << std::endl;
            return false;
        }
    }

	return true;
}

void AddAttributeToCells(vtkSmartPointer<vtkUnstructuredGrid> pGrid, vtkSmartPointer<vtkIdList> pCellCounts)
{
	// Setup index array
	vtkSmartPointer<vtkIntArray> indexArray = vtkSmartPointer<vtkIntArray>::New();
	indexArray->SetNumberOfComponents(1);
	indexArray->SetName("Part ID");
	int partId = 0;
	for (vtkIdType i = 0; i < pGrid->GetNumberOfCells(); i++)
	{
		if (i < pCellCounts->GetId(partId))
		{
			indexArray->InsertNextValue(partId);
		}
		else
		{
			indexArray->InsertNextValue(++partId);
		}
	}
	pGrid->GetCellData()->SetScalars(indexArray);
}

void ExportAsAbaqusFile(vtkSmartPointer<vtkUnstructuredGrid> ipGrid, vtkSmartPointer<vtkIdList> pCellCounts, std::string filePath)
{
	std::ofstream file;
	file.open(filePath);
	file << "*************************************\n";
	file << "*HEADING\n";
	file << "FemTech\n";
	file << "*************************************\n";
	file << "*NODE, NSET=All\n";
	for (int i = 0; i < ipGrid->GetNumberOfPoints(); i++)
	{
		double pnt[3];
		ipGrid->GetPoint(i, pnt);
		file << i + 1 << ",  " << pnt[0] << ",  " << pnt[1] << ",  " << pnt[2]<<"\n";
	}
	int currentCell = 0;
	int nSet = pCellCounts->GetNumberOfIds();
	for (int set = 0; set < nSet; set++)
	{
		file << "*ELEMENT,TYPE=C3D8,ELSET=PART_" << set+1 << std::endl;
		int ncell = pCellCounts->GetId(set);
		for (int i = currentCell; i < pCellCounts->GetId(set); i++)
		{
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
			ipGrid->GetCellPoints(i, cellPointIds);
			// Only for Element with 8 nodes 
			assert( cellPointIds->GetNumberOfIds() == 8); 
			file << i+1 << ", " << cellPointIds->GetId(0) << ", " << cellPointIds->GetId(1) << ", " << cellPointIds->GetId(2) << ", " << cellPointIds->GetId(3) << ", " << cellPointIds->GetId(4) << ", " << cellPointIds->GetId(5) << ", " << cellPointIds->GetId(6) << ", " << cellPointIds->GetId(7) << std::endl;
		}
		currentCell = pCellCounts->GetId(set);
	}
	file.close();
}

int main(int argc, char *argv[])
{
	std::vector<std::string> inputFiles;
	std::string outFile;
	bool isInputs = false;
    bool isOutput = false;
	bool writeinp = false;

    for (int i = 1; i < argc; ++i) 
    {
        if (strcmp(argv[i],"-in") == 0)
        {
            isInputs = true;
            isOutput = false;
        }
        else if (strcmp(argv[i],"-out") == 0)
        {
            isInputs = false;
            isOutput = true;
        }
		else if (strcmp(argv[i], "-abaqus") == 0)
		{
			writeinp = true;
		}
        else
        {
            if (isInputs)
            {
                inputFiles.push_back(argv[i]);
            }
            else if(isOutput)
            {
                outFile = argv[i];
                isInputs = false;
                isOutput = false;
            }
        }
    }

	if (!ValidateInputFiles(inputFiles, outFile))
		return 0;

	int partNumber = 0;
	vtkSmartPointer<vtkIdList> pCellCounts = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkAppendFilter> pAppend = vtkSmartPointer<vtkAppendFilter>::New();
	std::cout << "Number of input files " << inputFiles.size() << "\n";
	for (int f=0; f<inputFiles.size(); f++)
	{
		cout << "Reading file " << inputFiles[f] << endl;
		vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
		reader->SetFileName(inputFiles[f].c_str());
		reader->Update();

		vtkSmartPointer<vtkUnstructuredGrid> pUnstructedGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		pUnstructedGrid->DeepCopy(reader->GetOutput());

		pAppend->AddInputData(pUnstructedGrid);

		int nCells = pUnstructedGrid->GetNumberOfCells();
		int nPrevCount = f > 0 ? pCellCounts->GetId(f-1) : 0;
		pCellCounts->InsertNextId(nCells + nPrevCount);
	}
	pAppend->Update();

	vtkSmartPointer<vtkCleanUnstructuredGrid> pGridCleaner = vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
	pGridCleaner->AddInputData(pAppend->GetOutput());
	pGridCleaner->Update();
	vtkSmartPointer<vtkUnstructuredGrid> pMergedGrid = pGridCleaner->GetOutput();

	AddAttributeToCells(pMergedGrid, pCellCounts);

	vtkSmartPointer<vtkUnstructuredGridWriter> writter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writter->SetInputData(pMergedGrid);
	writter->SetFileName(outFile.c_str());
	writter->Write();
	std::cout << "Written to : " << outFile<< std::endl;

	if (writeinp)
	{
		std::string abaqusFilename = outFile.substr(0, outFile.find_last_of('.')) + ".inp";
		ExportAsAbaqusFile(pMergedGrid, pCellCounts, abaqusFilename);
		std::cout << "Written Abaqus inp file : " << abaqusFilename << std::endl;
	}
}
