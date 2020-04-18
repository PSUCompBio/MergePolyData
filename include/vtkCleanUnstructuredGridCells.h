#ifndef vtkCleanUnstructuredGridCells_h
#define vtkCleanUnstructuredGridCells_h

#include "vtkUnstructuredGridAlgorithm.h"

class  vtkCleanUnstructuredGridCells
  : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkCleanUnstructuredGridCells* New();

  vtkTypeMacro(vtkCleanUnstructuredGridCells, vtkUnstructuredGridAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkCleanUnstructuredGridCells();
  ~vtkCleanUnstructuredGridCells();

  virtual int RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  virtual int FillInputPortInformation(int port, vtkInformation* info) override;

private:
  vtkCleanUnstructuredGridCells(const vtkCleanUnstructuredGridCells&) = delete;
  void operator=(const vtkCleanUnstructuredGridCells&) = delete;
};

#endif
