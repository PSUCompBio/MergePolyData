#ifndef vtkCleanUnstructuredGrid_h
#define vtkCleanUnstructuredGrid_h

#include "vtkUnstructuredGridAlgorithm.h"

class vtkPointLocator;

class vtkCleanUnstructuredGrid : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkCleanUnstructuredGrid* New();

  vtkTypeMacro(vtkCleanUnstructuredGrid, vtkUnstructuredGridAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkCleanUnstructuredGrid();
  ~vtkCleanUnstructuredGrid() override;

  vtkPointLocator* Locator;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillInputPortInformation(int port, vtkInformation* info) override;

private:
  vtkCleanUnstructuredGrid(const vtkCleanUnstructuredGrid&) = delete;
  void operator=(const vtkCleanUnstructuredGrid&) = delete;
};
#endif
