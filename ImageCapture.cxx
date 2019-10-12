#include <vector>
#include <cstring>

#include "vtkAutoInit.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPLYReader.h>
#include <vtkJPEGReader.h>
#include <vtkProperty.h>
#include <vtkImageReader2Factory.h>
#include <vtkWindowToImageFilter.h>
#include <vtkCamera.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkTexture.h>

#include "include/tinyply.h"
using namespace tinyply;

struct uint3 { uint32_t v1, v2, v3; };
struct float2 { float u, v; };
struct float6 { float u1, v1, u2, v2, u3, v3; };

std::vector<float2> readTexCoordsFromPLYFile(const std::string & filepath)
{

  std::vector<uint3> faces;
  std::map<int, float2> vertTexMap;
	try
	{
		std::ifstream ss(filepath, std::ios::binary);
		if (ss.fail()) throw std::runtime_error("failed to open " + filepath);

		PlyFile file;
		file.parse_header(ss);

		std::shared_ptr<PlyData> pTexcoords, pFaces;

    try { pFaces = file.request_properties_from_element("face", { "vertex_indices" }, 3); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { pTexcoords = file.request_properties_from_element("face", { "texcoord" }, 6); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		file.read(ss);
		if (pTexcoords)
    {
			const size_t numTexsBytes = pTexcoords->buffer.size_bytes();
			std::vector<float6> texCoords(pTexcoords->count);
			std::memcpy(texCoords.data(), pTexcoords->buffer.get(), numTexsBytes);

      if (pFaces)
      {
        const size_t numFacesBytes = pFaces->buffer.size_bytes();
		  	faces.resize(pFaces->count);
			  std::memcpy(faces.data(), pFaces->buffer.get(), numFacesBytes);

        for(int f =0; f<faces.size(); f++)
        {
          uint3 verts = faces[f];
          float6 texs = texCoords[f];
          vertTexMap[verts.v1] = float2{texs.u1, texs.v1};
          vertTexMap[verts.v2] = float2{texs.u2, texs.v2};
          vertTexMap[verts.v3] = float2{texs.u3, texs.v3};
        }
      }
		}
	}
	catch (const std::exception & e)
	{
		std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
	}

  std::vector<float2> texscoords;
  for(auto mit:vertTexMap)
  {
    texscoords.push_back(mit.second);
  }
  faces.clear();
  vertTexMap.clear();
  std::cout<<"Size of Textures "<<texscoords.size()<<std::endl;
  return texscoords;
}

int main(int argc, char *argv[]) {
  if(argc < 4 || argc > 5)
  {
    std::cout << "Usage: " << argv[0] << "  FileName(.ply)" <<"   ImageName(.jpg/.png)"<<"   outputImage(.png)"<<"   magnification(optional)"<<std::endl;
    return EXIT_FAILURE;
  }

  //Read PLY file
  std::string inputFileName = argv[1];
  std::string inputImageName = argv[2];
  std::string outputImageName = argv[3];
  float magnification = 2.0;
  if (argc == 5)
    magnification = std::max(4.0,atof(argv[4]));

  vtkSmartPointer<vtkPLYReader> plyReader = vtkSmartPointer<vtkPLYReader>::New();
  plyReader->SetFileName(inputFileName.c_str());
  plyReader->Update();

  vtkSmartPointer<vtkPolyData> polyData = plyReader->GetOutput();

  vtkSmartPointer<vtkDataArray> texArray = polyData->GetPointData()->GetTCoords();
  if (texArray == nullptr)
  {
    std::vector<float2> texCoords = readTexCoordsFromPLYFile(inputFileName);
    vtkSmartPointer<vtkFloatArray> textureCoordinates = vtkSmartPointer<vtkFloatArray>::New();
    textureCoordinates->SetNumberOfComponents(2);
    textureCoordinates->SetName("TextureCoordinates");

    float tuple[2] = {0.0, 0.0};
    for(auto tIt : texCoords)
    {
      tuple[0] = tIt.u; tuple[1] = tIt.v;
      textureCoordinates->InsertNextTuple(tuple);
    }
    texCoords.clear();

    polyData->GetPointData()->SetTCoords(textureCoordinates);
  }

  //read image
  vtkSmartPointer<vtkImageReader2Factory> readerFactory =	vtkSmartPointer<vtkImageReader2Factory>::New();
  vtkImageReader2 * imageReader = readerFactory->CreateImageReader2(inputImageName.c_str());
  vtkSmartPointer<vtkImageReader2> smartImageReader;
  smartImageReader.TakeReference(imageReader);
  imageReader->SetFileName(inputImageName.c_str());
  imageReader->Update();
  vtkSmartPointer<vtkTexture> texture = vtkSmartPointer<vtkTexture>::New();
  texture->SetInputConnection(imageReader->GetOutputPort());
  texture->Update();

  vtkSmartPointer<vtkPolyDataMapper> mapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(plyReader->GetOutputPort());
  mapper->ScalarVisibilityOn();

  vtkSmartPointer<vtkActor> actor =  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetTexture("normalTex", texture);

  vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetOffScreenRendering(1);
	renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->SetBackground(1, 1, 1);

  vtkSmartPointer<vtkCamera> pCamera = renderer->GetActiveCamera();
  pCamera->Azimuth(30);
  pCamera->Elevation(15);
  renderer->ResetCamera();
  pCamera->Dolly(1.4);
  pCamera->SetFocalPoint(0,0,0);

  renderer->ResetCameraClippingRange();
  float size = 300*magnification;
  renderWindow->SetSize(size, size);
	renderWindow->Render();

  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(renderWindow);
  windowToImageFilter->Update();

  vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName(outputImageName.c_str());
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();
}
