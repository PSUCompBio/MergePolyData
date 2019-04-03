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
#include <vtkCleanUnstructuredGrid.h>
#include <vtkCleanUnstructuredGridCells.h>
#include <vtkSmartPointer.h>
#include <vtkPLYReader.h>
#include <vtkMeshQuality.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkMath.h>

#include <vtkPolygon.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

struct Hexa
{
	Hexa(std::vector<vtkIdType> cell1Ids, std::vector<vtkIdType> cell2Ids)
	{
		cellIds = cell1Ids;
		cellIds.insert(cellIds.end(), cell2Ids.begin(), cell2Ids.end());
	}

	std::vector<vtkIdType> cellIds;

	/*inline bool operator==(Hexa aHex)
	{
		std::vector<vtkIdType> cells1 = cellIds;
		std::vector<vtkIdType> cells2 = aHex.cellIds;

		std::sort(cells1.begin(), cells1.end());
		std::sort(cells2.begin(), cells2.end());

		std::vector<vtkIdType> diff;
		std::set_difference(cells1.begin(), cells1.end(), cells2.begin(), cells2.end(), std::back_inserter(diff));

		return diff.size() > 0 ? false : true;
	}*/
};

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

vtkSmartPointer<vtkActor> GetCellActor(vtkSmartPointer<vtkPolyData> polyData, std::vector<vtkIdType> cellPoints)
{
	vtkIdType nPnts = cellPoints.size();
	if (nPnts == 4)
	{

		vtkIdType id1 = cellPoints[0];
		vtkIdType id2 = cellPoints[1];
		vtkIdType id3 = cellPoints[2];
		vtkIdType id4 = cellPoints[3];

		double pnt1[3], pnt2[3], pnt3[3], pnt4[3];
		polyData->GetPoint(id1, pnt1);
		polyData->GetPoint(id2, pnt2);
		polyData->GetPoint(id3, pnt3);
		polyData->GetPoint(id4, pnt4);

		// Setup four points
		vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();
		points->InsertNextPoint(pnt1);
		points->InsertNextPoint(pnt2);
		points->InsertNextPoint(pnt3);
		points->InsertNextPoint(pnt4);

		// Create the polygon
		vtkSmartPointer<vtkPolygon> polygon =
			vtkSmartPointer<vtkPolygon>::New();
		polygon->GetPointIds()->SetNumberOfIds(4); //make a quad
		polygon->GetPointIds()->SetId(0, 0);
		polygon->GetPointIds()->SetId(1, 1);
		polygon->GetPointIds()->SetId(2, 2);
		polygon->GetPointIds()->SetId(3, 3);

		// Add the polygon to a list of polygons
		vtkSmartPointer<vtkCellArray> polygons =
			vtkSmartPointer<vtkCellArray>::New();
		polygons->InsertNextCell(polygon);

		// Create a PolyData
		vtkSmartPointer<vtkPolyData> polygonPolyData =
			vtkSmartPointer<vtkPolyData>::New();
		polygonPolyData->SetPoints(points);
		polygonPolyData->SetPolys(polygons);

		// Create a mapper and actor
		vtkSmartPointer<vtkPolyDataMapper> mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(polygonPolyData);
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.0, 0.0, 1.0);
		actor->GetProperty()->SetLineWidth(3.0);
		return actor;
	}
	return nullptr;
}

void Visualize(vtkSmartPointer<vtkPolyData> pData, std::vector<vtkSmartPointer<vtkActor>> actors)
{
	// Visualize
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(pData);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	actor->GetProperty()->SetOpacity(0.5);
	renderer->AddActor(actor);
	for (auto it : actors)
	{
		renderer->AddActor(it);
	}
	renderer->SetBackground(0.1804, 0.5451, 0.3412); // Sea green

	renderWindow->Render();
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview

	renderWindowInteractor->SetInteractorStyle(style);
	renderWindowInteractor->Start();
}

void VisualizeUnstructedGrid(vtkSmartPointer<vtkUnstructuredGrid> ug)
{
	// Create a mapper and actor
	vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputData(ug);

	vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
	actor1->SetMapper(mapper);

	// Create a renderer, render window, and interactor
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actor to the scene
	actor1->GetProperty()->SetEdgeVisibility(true);
	renderer->AddActor(actor1);
	renderer->SetBackground(.3, .6, .3); // Background color green

	// Render and interact
	renderWindow->Render();
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview

	renderWindowInteractor->SetInteractorStyle(style);
	renderWindowInteractor->Start();
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

int GetConnectedCells(vtkSmartPointer<vtkPolyData> pData, vtkIdType cellId, std::vector<vtkIdType> faces)
{
	std::vector<std::vector<vtkIdType>> allcellNeighbours;
	for (auto it : faces)
	{
		std::vector<vtkIdType> cellNeighbours = GetNeighbours(pData, it);
		cellNeighbours.push_back(it);
		cellNeighbours.erase(std::remove(cellNeighbours.begin(), cellNeighbours.end(), cellId), cellNeighbours.end());
		allcellNeighbours.push_back(cellNeighbours);
	}

	int nNeighbours = allcellNeighbours.size();
	for (int n = 0; n < nNeighbours; n++)
	{
		int next = (n + 1) % nNeighbours;
		std::vector<vtkIdType> cellPoints1 = allcellNeighbours[n];
		std::vector<vtkIdType> cellPoints2 = allcellNeighbours[next];
		std::sort(cellPoints1.begin(), cellPoints1.end());
		std::sort(cellPoints2.begin(), cellPoints2.end());

		std::vector<vtkIdType> intersect;
		std::set_intersection(cellPoints1.begin(), cellPoints1.end(), cellPoints2.begin(), cellPoints2.end(),
			std::back_inserter(intersect));
		intersect.size();
	}
	return 0;
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

static int DumpQualityStats(vtkMeshQuality* iq, const char *arrayname)
{
	cout << "  cardinality: "
		<< iq->GetOutput()->GetFieldData()->GetArray(arrayname)->GetComponent(0, 4)
		<< "  , range: "
		<< iq->GetOutput()->GetFieldData()->GetArray(arrayname)->GetComponent(0, 0)
		<< "  -  "
		<< iq->GetOutput()->GetFieldData()->GetArray(arrayname)->GetComponent(0, 2)
		<< endl;

	cout << "  average: " << iq->GetOutput()->GetFieldData()->GetArray(arrayname)->GetComponent(0, 1)
		<< "  , standard deviation: "
		<< sqrt(fabs(iq->GetOutput()->GetFieldData()->GetArray(arrayname)->GetComponent(0, 3)))
		<< endl;

	return 0;
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
		//cout << "Cell " << cellId << "\t" << dotProd << endl;
		if (dotProd > 0.0)
		{
			//Try rotating back
			set1 = { c1Points[0],c2Points[3],c1Points[1] };
			set2 = { c1Points[1],c2Points[0],c2Points[3] };
			c1Normal = GetCellNormal(pData, set1);
			c2Normal = GetCellNormal(pData, set2);
			c2Start = 3;

			dotProd = (vtkMath::Dot(c1Normal, c2Normal));
			//cout << "Again Cell " << cellId << "\t" << dotProd << endl;
			if (dotProd > 0.0)
			{
				//Try rotating front
				set1 = { c1Points[0],c2Points[1],c1Points[1] };
				set2 = { c1Points[1],c2Points[2],c2Points[1] };
				c1Normal = GetCellNormal(pData, set1);
				c2Normal = GetCellNormal(pData, set2);
				c2Start = 1;

				dotProd = (vtkMath::Dot(c1Normal, c2Normal));
				//cout << "Again again Cell " << cellId << "\t" << dotProd << endl;
				//if (dotProd > 0.0)
				//	cout << "Catch it\n";
			}
		}

		//double dotProd = vtkMath::Dot(c1Normal, c2Normal);
		//cout << "Dot Product " << dotProd << endl;

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

	//vtkSmartPointer<vtkCleanUnstructuredGridCells> cleaner = vtkSmartPointer<vtkCleanUnstructuredGridCells>::New();
	//cleaner->SetInputData(ug);
	//cleaner->Update();
	//ug = cleaner->GetOutput();

	//for (int c = 0; c < ug->GetNumberOfCells(); c++)
	//{
	//	vtkCell* aCell = ug->GetCell(c);
	//	std::cout << "Cell " << c << std::endl;
	//	std::cout<<"Scaled Jacobian " << vtkMeshQuality::HexScaledJacobian(aCell) <<"\t"<<"Disortion "<< vtkMeshQuality::HexDistortion(aCell)<<"\t" << "Shape " << vtkMeshQuality::HexShape(aCell) <<"\t" << "Shear " << vtkMeshQuality::HexShear(aCell) << "\t" << "Volume " << vtkMeshQuality::HexVolume(aCell) << std::endl;
	//}

	/*vtkMeshQuality* iq = vtkMeshQuality::New();
	iq->SetInputData(ug);
	iq->SetHexQualityMeasureToDistortion();
	iq->SaveCellQualityOn();
	iq->Update();

	std::cout << " Distortion:" << endl;
	DumpQualityStats(iq, "Mesh Hexahedron Quality");
	std::cout << endl;*/

	return ug;
}

int main(int argc, char *argv[])
{
	//std::string inputFilename = argv[1];
	//std::string inputFilename = "E://projects//current//RubenCraft//files//quad2hex//chank.ply";
	std::string inputFilename = "E://projects//current//RubenCraft//files//quad2hex//Layer.ply";
	vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
	reader->SetFileName(inputFilename.c_str());
	reader->Update();
	vtkSmartPointer<vtkPolyData> pData = reader->GetOutput();

	//vtkSmartPointer<vtkCleanPolyData> polyCleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	//polyCleaner->SetInputConnection(reader->GetOutputPort());
	//polyCleaner->Update();

	//vtkSmartPointer<vtkPolyData> pData = polyCleaner->GetOutput();
	vtkSmartPointer<vtkCellData> pCellData = pData->GetCellData();

	vtkIdType nCells = pData->GetNumberOfCells();
	std::vector< vtkIdType> cellIds;
	for (vtkIdType c = 0; c < nCells; c++)
	{
		cellIds.push_back(c);
	}

	std::vector<std::pair<vtkIdType, vtkIdType>> facePair;
	std::set<Hexa> hexElements;
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

		//GetConnectedCells(pData, cellId, neighbors);

		for (auto it : sortedNeighbours)
		{
			//1
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
			std::vector<vtkIdType> cell2Ids = GetCellPointIds(pData, foundCell);

			//Hexa aHex(cell1Ids, cell2Ids);
			//hexElements.insert(aHex);

			 //remove cells to prevent duplication
			if (allCells.size() != 4)
			{
				std::cout << "Error" << endl;
			}

			selectedCells.insert(selectedCells.end(), allCells.begin(), allCells.end());
			selectedCells.push_back(cellId);
			selectedCells.push_back(foundCell);

			//std::vector<vtkSmartPointer<vtkActor>> actors;
			//auto act = GetCellActor(pData, CijPointIds);
			//act->GetProperty()->SetColor(1, 0, 0);
			//act->GetProperty()->SetLineWidth(3.0);
			//act->GetProperty()->EdgeVisibilityOn();
			//actors.push_back(act);

			//auto foundAct = GetCellActor(pData, GetCellPointIds(pData, foundCell));
			//foundAct->GetProperty()->SetColor(0, 0, 1);
			//foundAct->GetProperty()->SetLineWidth(3.0);
			//foundAct->GetProperty()->EdgeVisibilityOn();
			//actors.push_back(foundAct);
			//Visualize(pData, actors);

			facePair.push_back(std::pair<vtkIdType, vtkIdType>(cellId, foundCell));
			break;
		}
	}

	vtkSmartPointer<vtkUnstructuredGrid> ug = BuildHexElements(pData, facePair);
	vtkSmartPointer<vtkCleanUnstructuredGrid> gridCleaner = vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
	gridCleaner->AddInputData(ug);
	gridCleaner->Update();

	ug = gridCleaner->GetOutput();

	//std::string outFile = "E://projects//current//RubenCraft//files//quad2hex//chank_Hex.vtk";
	std::string outFile = "E://projects//current//RubenCraft//files//quad2hex//layer_Hex.vtk";
	vtkSmartPointer<vtkUnstructuredGridWriter> writter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writter->SetInputData(ug);
	writter->SetFileName(outFile.c_str());
	writter->Write();
	std::cout << "Written to : " << outFile << std::endl;

	VisualizeUnstructedGrid(ug);
}
