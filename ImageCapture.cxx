#include "vtkPLYReader.h"

#include "vtkActor.h"
#include "vtkPNGReader.h"
#include <vtkJPEGReader.h>
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRegressionTestImage.h"
#include "vtkTestUtilities.h"
#include "vtkTexture.h"

int main ( int argc, char *argv[] )
{
  // Read file name.
  if (argc < 2)
  {
    return EXIT_FAILURE;
  }
  std::string fn = "/media/vishnu/Local Disk/projects/current/RubenCraft/imageSave/";
  std::string plyName = fn + argv[1];
  std::string imageName = fn + argv[2];
  std::cout<<"ply name : "<<plyName<<std::endl;
  std::cout<<"image name : "<<imageName<<std::endl;
  const char* fname = vtkTestUtilities::ExpandDataFileName(
    argc, argv, plyName.c_str());
  const char* fnameImg = vtkTestUtilities::ExpandDataFileName(
    argc, argv, imageName.c_str());

  // Test if the reader thinks it can open the file.
//   if (0 == vtkPLYReader::CanReadFile(fname))
//   {
//     std::cout << "The PLY reader can not read the input file." << std::endl;
//     return EXIT_FAILURE;
//   }

  // Create the reader.
  vtkPLYReader* reader = vtkPLYReader::New();
  reader->SetFileName(plyName.c_str());
  reader->Update();
//   delete [] fname;

//   vtkPNGReader* readerImg = vtkPNGReader::New();
//   if (0 == readerImg->CanReadFile(fnameImg))
//   {
//      std::cout << "The PNG reader can not read the input file." << std::endl;
//      return EXIT_FAILURE;
//   }
//   readerImg->SetFileName(imageName.c_str());
//   readerImg->Update();
//   delete[] fnameImg;

  vtkSmartPointer<vtkJPEGReader> jpegReader = vtkSmartPointer<vtkJPEGReader>::New();
  jpegReader->SetFileName(imageName.c_str());
  jpegReader->Update();

  // Create the texture.
  vtkTexture* texture = vtkTexture::New();
//   texture->SetInputConnection(readerImg->GetOutputPort());
  texture->SetInputConnection(jpegReader->GetOutputPort());
  texture->InterpolateOn();

  // Create a mapper.
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  mapper->SetInputConnection(reader->GetOutputPort());
  mapper->ScalarVisibilityOn();

  // Create the actor.
  vtkActor* actor = vtkActor::New();
  actor->SetMapper(mapper);
  actor->SetTexture(texture);


  // Basic visualisation.
  vtkRenderWindow* renWin = vtkRenderWindow::New();
  vtkRenderer* ren = vtkRenderer::New();
  renWin->AddRenderer(ren);
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);

  ren->AddActor(actor);
  ren->SetBackground(0.1,0.2,0.3);
  renWin->SetSize(400,400);

  // interact with data
  renWin->Render();

  int retVal = vtkRegressionTestImage( renWin );

//   if ( retVal == vtkRegressionTester::DO_INTERACTOR)
  {
    iren->Start();
  }

  actor->Delete();
  mapper->Delete();
  reader->Delete();
//   readerImg->Delete();
//   jpegReader->Delete();
//   texture->Delete();
  renWin->Delete();
  ren->Delete();
  iren->Delete();

  return !retVal;
}