#include <vector>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "json/json.h"

#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL);
VTK_MODULE_INIT(vtkInteractionStyle);

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPLYReader.h>
#include <vtkOBJReader.h>
#include <vtkJPEGReader.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGridReader.h>
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

Json::Value getConfig(const char* inputFile);

std::vector<float2> readTexCoordsFromPLYFile(const std::string & filepath)
{
	try
	{
		std::ifstream ss(filepath, std::ios::binary);
		if (ss.fail()) throw std::runtime_error("failed to open " + filepath);

		PlyFile file;
		file.parse_header(ss);

		std::shared_ptr<PlyData> pTexcoords, pFaces;

        try { pTexcoords = file.request_properties_from_element("vertex", { "s", "t" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		file.read(ss);
		if (pTexcoords)
        {
			const size_t numTexsBytes = pTexcoords->buffer.size_bytes();
            const size_t nTexs = pTexcoords->count;
			std::vector<float2> texCoords(nTexs);
			std::memcpy(texCoords.data(), pTexcoords->buffer.get(), numTexsBytes);

            std::cout<<"Size of Textures "<<texCoords.size()<<std::endl;
            return texCoords;
        }
	}
	catch (const std::exception & e)
	{
		std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
	}
}

int main(int argc, char *argv[]) {
    if(argc < 5 || argc > 6)
    {
        std::cout << "Usage: " << argv[0] << "  meshFile(.ply)" <<"   textureFile(.jpg/.png)"<<"   datFile(.dat)"<<"   outputImage(.jpg/.png)"<<"   magnification(optional)"<<std::endl;
        return EXIT_FAILURE;
    }

    std::string inputFileName = argv[1];
    std::string inputImageName = argv[2];
    // std::string inputDatName = argv[3];
    std::string outputImageName = argv[4];
    float magnification = 2.0;
    if (argc == 6)
        magnification = std::max(4.0,atof(argv[5]));

    //mesh file
    std::string ext = inputFileName.substr(inputFileName.find_last_of(".") + 1);
    double bounds[6];

    vtkSmartPointer<vtkActor> pActor;
    if(ext.compare("ply") == 0 || ext.compare("PLY") == 0)
    {
        std::cout<<"opening ply file \n";
        vtkSmartPointer<vtkPLYReader> pReader = vtkSmartPointer<vtkPLYReader>::New();
        pReader->SetFileName(inputFileName.c_str());
        pReader->Update();
        vtkSmartPointer<vtkPolyData> polyData = pReader->GetOutput();
        polyData->GetBounds(bounds);

        vtkSmartPointer<vtkDataArray> texArray = polyData->GetPointData()->GetTCoords();
        if(texArray == nullptr)
        {
            std::vector<float2> texCoords = readTexCoordsFromPLYFile(inputFileName);
            if(texCoords.size())
            {
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
        }
        else
        {
            std::cout<<"Avilable texture !\n";
        }

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);
        mapper->ScalarVisibilityOn();

        pActor = vtkSmartPointer<vtkActor>::New();
        pActor->SetMapper(mapper);
        pActor->GetProperty()->SetOpacity(0.25);

        texArray = polyData->GetPointData()->GetTCoords();
        if(texArray != nullptr)
        {
            vtkSmartPointer<vtkImageReader2Factory> readerFactory =	vtkSmartPointer<vtkImageReader2Factory>::New();
            vtkImageReader2 * imageReader = readerFactory->CreateImageReader2(inputImageName.c_str());
            vtkSmartPointer<vtkImageReader2> smartImageReader;
            smartImageReader.TakeReference(imageReader);
            imageReader->SetFileName(inputImageName.c_str());
            imageReader->Update();
            vtkSmartPointer<vtkTexture> texture = vtkSmartPointer<vtkTexture>::New();
            texture->SetInputConnection(imageReader->GetOutputPort());
            texture->Update();

            pActor->GetProperty()->SetTexture("normalTex", texture);
        }
    }
    else
    {
        std::cout<<"Not a ply file \n";
    }

    //dat file
    Json::Value outputJson = getConfig(argv[3]);
    std::vector<vtkSmartPointer<vtkActor>> maxSpheres;
    std::vector<vtkSmartPointer<vtkActor>> minSpheres;
    float sphereRadius = (bounds[1] - bounds[0])/25.0;

    float MaxX, MaxY, MaxZ, MinX, MinY, MinZ;
    MaxX = outputJson["principal-max-strain"]["location"][0].asDouble();
    MaxY = outputJson["principal-max-strain"]["location"][1].asDouble();
    MaxZ = outputJson["principal-max-strain"]["location"][2].asDouble();
    MinX = outputJson["principal-min-strain"]["location"][0].asDouble();
    MinY = outputJson["principal-min-strain"]["location"][1].asDouble();
    MinZ = outputJson["principal-min-strain"]["location"][2].asDouble();

    vtkSmartPointer<vtkSphereSource> maxSource = vtkSmartPointer<vtkSphereSource>::New();
    maxSource->SetCenter(MaxX, MaxY, MaxZ);
    maxSource->SetThetaResolution(64);
    maxSource->SetPhiResolution(64);
    maxSource->SetRadius(sphereRadius);
    maxSource->Update();
    vtkSmartPointer<vtkPolyDataMapper> maxMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    maxMapper->SetInputData(maxSource->GetOutput());
    vtkSmartPointer<vtkActor> maxActor = vtkSmartPointer<vtkActor>::New();
    maxActor->SetMapper(maxMapper);
    maxActor->GetProperty()->SetColor(1.0,0.0,0.0);
    maxSpheres.push_back(maxActor);

    vtkSmartPointer<vtkSphereSource> minSource = vtkSmartPointer<vtkSphereSource>::New();
    minSource->SetCenter(MinX, MinY, MinZ);
    minSource->SetThetaResolution(64);
    minSource->SetPhiResolution(64);
    minSource->SetRadius(sphereRadius);
    minSource->Update();
    vtkSmartPointer<vtkPolyDataMapper> minMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    minMapper->SetInputData(minSource->GetOutput());
    vtkSmartPointer<vtkActor> minActor = vtkSmartPointer<vtkActor>::New();
    minActor->SetMapper(minMapper);
    minActor->GetProperty()->SetColor(0.0,0.0,1.0);
    minSpheres.push_back(minActor);

    // Define viewport ranges
    double xmins[4] = {0,.5,0,.5};
    double xmaxs[4] = {0.5,1,0.5,1};
    double ymins[4] = {0,0,.5,.5};
    double ymaxs[4]= {0.5,0.5,1,1};
    double camera_views[4][3] = {{-1.0,0.0,0.0},{-1.0,0.70,0.7},{0.0,0.0,-1.0},{0.0,-1.0,0.0}};
    double camera_viewUp[4][3] = {{0.0,1.0,0.0},{0.0,1.0,0.0},{0.0,1.0,0.0},{0.0,-1.0,0.0}};
    vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();

    for(unsigned i = 0; i < 4; i++)
    {
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        renderWindow->AddRenderer(renderer);
        renderer->SetViewport(xmins[i],ymins[i],xmaxs[i],ymaxs[i]);
        renderer->AddActor(pActor);
        for(auto maxIt:maxSpheres)
        {
            renderer->AddActor(maxIt);
        }
        for(auto minIt:minSpheres)
        {
            renderer->AddActor(minIt);
        }
        renderer->SetBackground(1, 1, 1);

        vtkSmartPointer<vtkCamera> pCamera = renderer->GetActiveCamera();
        pCamera->SetPosition(camera_views[i]);
        pCamera->SetViewUp(camera_viewUp[i]);
        renderer->ResetCamera();
    }

    float size = 400*magnification;
    renderWindow->SetSize(size, size);
    renderWindow->SetOffScreenRendering(1);
    renderWindow->Render();

    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(outputImageName.c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();

    return EXIT_SUCCESS;
}

Json::Value getConfig(const char* inputFile) {
  Json::Value root;
  std::ifstream ifs;
  ifs.open(inputFile);
  if (!ifs.is_open()) {
    std::cout<<"ERROR : Failed to open configuration file" << std::endl;
    exit(1);
  }

  Json::CharReaderBuilder builder;
  JSONCPP_STRING errs;
  if (!parseFromStream(builder, ifs, &root, &errs)) {
    std::cout<<"ERROR : " << errs << std::endl;
    exit(1);
  }
  return root;
}
