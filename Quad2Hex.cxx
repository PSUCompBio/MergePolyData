#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkIntArray.h>
#include <vtkCleanPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPLYReader.h>
#include <vtkMeshQuality.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkMath.h>

double* GetCellNormal(vtkSmartPointer<vtkPolyData> polyData, std::vector<vtkIdType> cellPoints)
{
	double* normal = nullptr;
	vtkIdType nPnts = cellPoints.size();
	if (nPnts > 2)
	{
		vtkIdType id1 = cellPoints[0];
		vtkIdType id2 = cellPoints[1];
		vtkIdType id3 = cellPoints[nPnts - 1];

		double pnt1[3], pnt2[3], pnt3[3], vec1[3], vec2[3];
		polyData->GetPoint(id1, pnt1);
		polyData->GetPoint(id2, pnt2);
		polyData->GetPoint(id3, pnt3);

		vtkMath::Subtract(pnt2, pnt1, vec1);
		vtkMath::Subtract(pnt2, pnt3, vec2);

		normal = new double(3);
		vtkMath::Cross(vec1, vec2, normal);
		vtkMath::Normalize(normal);
	}
	return normal;
}

std::vector<vtkIdType> GetNeighbours(vtkSmartPointer<vtkPolyData> pData, vtkIdType cellId)
{

	std::vector<vtkIdType> neighbors;

	vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
	pData->GetCellPoints(cellId, cellPointIds);
	vtkIdType nPnts = cellPointIds->GetNumberOfIds();
	for (vtkIdType c = 0; c < nPnts; c++)
	{
		vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
		idList->InsertNextId(cellPointIds->GetId(c));
		vtkIdType next = (c + 1) % nPnts;
		idList->InsertNextId(cellPointIds->GetId(next));

		//get the neighbors of the cell
		vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
		pData->GetCellNeighbors(cellId, idList, neighborCellIds);
		vtkIdType nNeighCells = neighborCellIds->GetNumberOfIds();

		for (vtkIdType j = 0; j < nNeighCells; j++)
		{
			neighbors.push_back(neighborCellIds->GetId(j));
		}
	}

	return neighbors;
}

std::vector<vtkIdType> GetCellPointIds(vtkSmartPointer<vtkPolyData> pData, vtkIdType cellId)
{
	vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
	pData->GetCellPoints(cellId, cellPoints);

	std::vector<vtkIdType> cellPointIds;
	for (vtkIdType n = 0; n < cellPoints->GetNumberOfIds(); n++)
	{
		cellPointIds.push_back(cellPoints->GetId(n));
	}

	return cellPointIds;
}

std::vector<vtkIdType> FindCellPointsNotOnCellij(std::vector<vtkIdType> C1, std::vector<vtkIdType> neighbourCell)
{
	std::vector<vtkIdType> cellPnts;
	for (auto it : neighbourCell)
	{
		if (std::find(C1.begin(), C1.end(), it) == C1.end())
		{
			cellPnts.push_back(it);
		}
	}
	return cellPnts;
}

std::vector<vtkIdType> FindCellPointsNotOnCells(std::vector<vtkIdType> C1, std::vector<vtkIdType> C2, std::vector<vtkIdType> neighbourCell)
{
	std::vector<vtkIdType> cellPnts;
	for (auto it : neighbourCell)
	{
		if (std::find(C1.begin(), C1.end(), it) == C1.end())
		{
			if (std::find(C2.begin(), C2.end(), it) == C2.end())
			{
				cellPnts.push_back(it);
			}
		}
	}
	return cellPnts;
}

vtkIdType FindCommonCell(std::vector<std::vector<vtkIdType>> Cells)
{
	int n = Cells.size();
	std::vector<vtkIdType> intersect = Cells[0];
	for (int i = 1; i < n; i++)
	{
		std::vector<vtkIdType> currCells = Cells[i];
		std::sort(currCells.begin(), currCells.end());
		std::sort(intersect.begin(), intersect.end());

		std::vector<vtkIdType> inter;
		std::set_intersection(currCells.begin(), currCells.end(), intersect.begin(), intersect.end(),
			std::back_inserter(inter));
		intersect.clear();
		intersect = inter;
	}

	return intersect[0];
}

vtkSmartPointer<vtkUnstructuredGrid> BuildHexElements(vtkSmartPointer<vtkPolyData> pData, std::vector<vtkSmartPointer<vtkIdList>> idLists)
{
	vtkSmartPointer<vtkPoints> pPoints = pData->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
	ug->SetPoints(pPoints);

	for (auto it : idLists)
	{
		ug->InsertNextCell(VTK_HEXAHEDRON, it.GetPointer());
	}
	return ug;
}

vtkSmartPointer<vtkUnstructuredGrid> BuildHexElements(vtkSmartPointer<vtkPolyData> pData, std::vector<std::pair<vtkIdType, vtkIdType>> cellPairs)
{
	vtkSmartPointer<vtkPoints> pPoints = pData->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
	ug->SetPoints(pPoints);
	vtkIdType cellId = 0;
	vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
	for (auto it : cellPairs)
	{
		std::vector<vtkIdType> c1Points = GetCellPointIds(pData, it.first);
		std::vector<vtkIdType> c2Points = GetCellPointIds(pData, it.second);

		std::vector<vtkIdType> set1{ c1Points[0],c2Points[0],c1Points[1] };
		std::vector<vtkIdType> set2{ c1Points[1],c2Points[1],c2Points[0] };
		double* c1Normal = GetCellNormal(pData, set1);
		double* c2Normal = GetCellNormal(pData, set2);

		vtkIdType c2Start = 0;
		double dotProd = (vtkMath::Dot(c1Normal, c2Normal));
		if (dotProd > 0.0)
		{
			//Try rotating back
			set1 = { c1Points[0],c2Points[3],c1Points[1] };
			set2 = { c1Points[1],c2Points[0],c2Points[3] };
			c1Normal = GetCellNormal(pData, set1);
			c2Normal = GetCellNormal(pData, set2);
			c2Start = 3;

			dotProd = (vtkMath::Dot(c1Normal, c2Normal));
			if (dotProd > 0.0)
			{
				//Try rotating front
				set1 = { c1Points[0],c2Points[1],c1Points[1] };
				set2 = { c1Points[1],c2Points[2],c2Points[1] };
				c1Normal = GetCellNormal(pData, set1);
				c2Normal = GetCellNormal(pData, set2);
				c2Start = 1;

				dotProd = (vtkMath::Dot(c1Normal, c2Normal));
			}
		}

		for (auto p1it : c1Points)
		{
			ids->InsertNextId(p1it);
		}

		for (int p = 0; p < c2Points.size(); p++)
		{
			vtkIdType index = (p + c2Start) % c2Points.size();
			ids->InsertNextId(c2Points[index]);
		}

		ug->InsertNextCell(VTK_HEXAHEDRON, ids.GetPointer());
		ids->Reset();

		//vtkCell* aCell = ug->GetCell(cellId);
		//std::cout << "Cell " << cellId << std::endl;
		//std::cout<<"Scaled Jacobian " << vtkMeshQuality::HexScaledJacobian(aCell) <<"\t"<<"Disortion "<< vtkMeshQuality::HexDistortion(aCell)<<"\t" << "Shape " << vtkMeshQuality::HexShape(aCell) <<"\t" << "Shear " << vtkMeshQuality::HexShear(aCell) << "\t" << "Volume " << vtkMeshQuality::HexVolume(aCell) << std::endl;
		cellId++;
	}

	return ug;
}

vtkIdType FindOppositeCellPointId(std::vector<std::vector<vtkIdType>> allCellPointIds, vtkIdType seedPointId)
{
	std::vector<std::vector<vtkIdType>> aNeighCellsPointIds;
	for(auto it: allCellPointIds)
	{
		if (std::find(it.begin(), it.end(), seedPointId) != it.end())
		{
			aNeighCellsPointIds.push_back(it);
		}
	}
	std::sort(aNeighCellsPointIds[0].begin(), aNeighCellsPointIds[0].end());
	std::sort(aNeighCellsPointIds[1].begin(), aNeighCellsPointIds[1].end());

	std::vector<vtkIdType> inter;
	std::set_intersection(aNeighCellsPointIds[0].begin(), aNeighCellsPointIds[0].end(), aNeighCellsPointIds[1].begin(), aNeighCellsPointIds[1].end(),
	std::back_inserter(inter));
	inter.erase(std::remove(inter.begin(), inter.end(), seedPointId), inter.end());

	if (inter.size() > 1)
	{
	std::cout<<"Error : Couldnt find the Opposite cell point\n";
	}

	return inter[0];
}

int main(int argc, char *argv[])
{
	bool isInput = false;
    bool isOutput = false;
	std::string inputFilename;
	std::string outputFilename;
	for (int i = 1; i < argc; ++i) 
    {
		if (strcmp(argv[i],"-in") == 0)
        {
			isInput = true;
		}
		else if(strcmp(argv[i],"-out") == 0)
		{
			isInput = false;
			isOutput = true;
		}
		else
		{
			if (isInput)
			{
				inputFilename = argv[i];
			}
			else if (isOutput)
			{
				outputFilename = argv[i];
			}
		}
	}
	
	vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
	reader->SetFileName(inputFilename.c_str());
	reader->Update();
	vtkSmartPointer<vtkPolyData> pData = reader->GetOutput();
	vtkSmartPointer<vtkCellData> pCellData = pData->GetCellData();

	vtkIdType nCells = pData->GetNumberOfCells();
	std::vector< vtkIdType> cellIds;
	for (vtkIdType c = 0; c < nCells; c++)
	{
		cellIds.push_back(c);
	}

	std::cout<<"Number of cells "<<nCells<<"\n";
	std::vector<std::pair<vtkIdType, vtkIdType>> facePair;
	std::vector<vtkSmartPointer<vtkIdList>> hexIdLists;
	std::vector<vtkIdType> selectedCells;
	for (vtkIdType c = 0; c < cellIds.size(); c++)
	{
		vtkIdType cellId = cellIds[c];

		if (std::find(selectedCells.begin(), selectedCells.end(), cellId) != selectedCells.end())
			continue;

		std::vector<vtkIdType> CijPointIds = GetCellPointIds(pData, cellId);
		double* cellNormal = GetCellNormal(pData, CijPointIds);

		std::vector<vtkIdType> neighbors = GetNeighbours(pData, cellId);
		std::map<vtkIdType, std::vector<vtkIdType>> mapCellPnts;
		std::map<double, vtkIdType> sortedNeighbours;
		for (int n = 0; n < neighbors.size(); n++)
		{
			std::vector<vtkIdType> neighPointIds = GetCellPointIds(pData, neighbors[n]);
			double* normal = GetCellNormal(pData, neighPointIds);
			double dotProd = abs(vtkMath::Dot(cellNormal, normal));
			sortedNeighbours[dotProd] = neighbors[n];
			mapCellPnts[neighbors[n]] = neighPointIds;
		}

		for (auto it : sortedNeighbours)
		{
			vtkIdType C1 = it.second;
			std::vector<vtkIdType> C1PointIds = GetCellPointIds(pData, C1);
			std::vector<vtkIdType> filteredNeigh = neighbors;
			filteredNeigh.erase(std::remove(filteredNeigh.begin(), filteredNeigh.end(), C1), filteredNeigh.end());

			//Points not on Cij
			std::vector<vtkIdType> pnts1 = FindCellPointsNotOnCellij(CijPointIds, C1PointIds);

			// Find the cell which shares any one of the pnt
			vtkIdType c2 = -1;
			for (auto pit : pnts1)
			{
				for (auto nit : filteredNeigh)
				{
					std::vector<vtkIdType> neighbourCellsPnts = mapCellPnts[nit];
					if (std::find(neighbourCellsPnts.begin(), neighbourCellsPnts.end(), pit) != neighbourCellsPnts.end())
					{
						c2 = nit;
						break;
					}
				}
				if (c2 > -1)
					break;
			}
			if (c2 == -1)
				continue;

			std::vector<vtkIdType> C2PointIds = GetCellPointIds(pData, c2);
			filteredNeigh = neighbors;
			filteredNeigh.erase(std::remove(filteredNeigh.begin(), filteredNeigh.end(), c2), filteredNeigh.end());

			//Points not on Cij and C1 
			std::vector<vtkIdType> pnts2 = FindCellPointsNotOnCells(CijPointIds, C1PointIds, C2PointIds);

			// Find the cell which shares any one of the pnt
			vtkIdType c3 = -1;
			for (auto pit : pnts2)
			{
				for (auto nit : filteredNeigh)
				{
					std::vector<vtkIdType> neighbourCellsPnts = mapCellPnts[nit];
					if (std::find(neighbourCellsPnts.begin(), neighbourCellsPnts.end(), pit) != neighbourCellsPnts.end())
					{
						c3 = nit;
						break;
					}
				}
				if (c3 > -1)
					break;
			}
			if (c3 == -1)
				continue;

			std::vector<vtkIdType> C3PointIds = GetCellPointIds(pData, c3);
			filteredNeigh = neighbors;
			filteredNeigh.erase(std::remove(filteredNeigh.begin(), filteredNeigh.end(), c3), filteredNeigh.end());

			//Points not on Cij and C2
			std::vector<vtkIdType> pnts3 = FindCellPointsNotOnCells(CijPointIds, C2PointIds, C3PointIds);

			// Find the cell which shares any one of the pnt
			vtkIdType c4 = -1;
			for (auto pit : pnts3)
			{
				for (auto nit : filteredNeigh)
				{
					std::vector<vtkIdType> neighbourCellsPnts = mapCellPnts[nit];
					if (std::find(neighbourCellsPnts.begin(), neighbourCellsPnts.end(), pit) != neighbourCellsPnts.end())
					{
						c4 = nit;
						break;
					}
				}
				if (c4 > -1)
					break;
			}
			if (c4 == -1)
				continue;

			std::vector<vtkIdType> C4PointIds = GetCellPointIds(pData, c4);
			filteredNeigh = neighbors;
			filteredNeigh.erase(std::remove(filteredNeigh.begin(), filteredNeigh.end(), c4), filteredNeigh.end());
			//Points not on Cij and C3
			std::vector<vtkIdType> pnts4 = FindCellPointsNotOnCells(CijPointIds, C3PointIds, C4PointIds);

			std::vector<vtkIdType> allCells = { C1, c2, c3, c4 };
			std::vector<std::vector<vtkIdType>> allCellNeighs;
			for (auto it : allCells)
			{
				std::vector<vtkIdType> aNeighs = GetNeighbours(pData, it);
				aNeighs.erase(std::remove(aNeighs.begin(), aNeighs.end(), cellId), aNeighs.end());
				allCellNeighs.push_back(aNeighs);
			}

			vtkIdType foundCell = FindCommonCell(allCellNeighs);

			std::vector<vtkIdType> cell1Ids = GetCellPointIds(pData, cellId);
			// std::vector<vtkIdType> cell2Ids = GetCellPointIds(pData, foundCell);
			std::vector<std::vector<vtkIdType>> allCellPointIds = { C1PointIds, C2PointIds, C3PointIds, C4PointIds };

			vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
			for(int i=0; i<cell1Ids.size(); i++)
			{
				idList->InsertNextId(cell1Ids[i]);
			}
			for(int i=0; i<cell1Ids.size(); i++)
			{
				vtkIdType oppPointId = FindOppositeCellPointId(allCellPointIds, cell1Ids[i]);
				idList->InsertNextId(oppPointId);
			}
			
			 //remove cells to prevent duplication
			selectedCells.insert(selectedCells.end(), allCells.begin(), allCells.end());
			selectedCells.push_back(cellId);
			selectedCells.push_back(foundCell);

			facePair.push_back(std::pair<vtkIdType, vtkIdType>(cellId, foundCell));
			hexIdLists.push_back(idList);
			break;
		}
	}

	// vtkSmartPointer<vtkUnstructuredGrid> ug = BuildHexElements(pData, facePair);
	vtkSmartPointer<vtkUnstructuredGrid> ug = BuildHexElements(pData, hexIdLists);

	// Write to out file
	vtkSmartPointer<vtkUnstructuredGridWriter> writter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writter->SetInputData(ug);
	writter->SetFileName(outputFilename.c_str());
	writter->Write();
	std::cout << "Written to : " << outputFilename << std::endl;
}
