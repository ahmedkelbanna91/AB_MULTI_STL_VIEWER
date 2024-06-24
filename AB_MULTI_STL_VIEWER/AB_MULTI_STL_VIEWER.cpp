#include <future> 
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include "rang.hpp"

#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAutoInit.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkIdList.h>
#include <vtkAnnotatedCubeActor.h>
#include <vtkSphereSource.h>
#include <vtkAxesActor.h>
#include <vtkDiskSource.h>
#include <vtkCaptionActor2D.h>
#include <vtkSmartPointer.h>
#include <vtkTextProperty.h>
#include <vtkHardwarePicker.h>
#include <vtkTextActor.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkProperty2D.h>
#include <vtkSliderWidget.h>
#include <vtkCommand.h>
#include <vtkPlaneSource.h>
#include <vtkPlaneWidget.h>
#include <vtkTexturedButtonRepresentation2D.h>
#include <vtkButtonWidget.h>
#include <vtkImageData.h>
#include <vtkFreeTypeTools.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkWin32OpenGLRenderWindow.h>
#include <vtkSplineWidget.h>
#include <vtkSplineWidget2.h>
#include <vtkSplineRepresentation.h>
#include <vtkKochanekSpline.h>
#include <vtkButtonWidget.h>
#include <vtkCoordinate.h>
#include <vtkImageData.h>
#include <vtkNamedColors.h>
#include <vtkCallbackCommand.h>
#include <vtkCameraOrientationWidget.h>
#include <vtkDoubleArray.h>
#include <vtkCameraOrientationRepresentation.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkFloatArray.h>
#include <fstream>
#include <sstream>

#define M_PI 3.14159265358979323846

//<< Red << << ColorEnd <<
auto ColorEnd = [](std::ostream& os) -> std::ostream& { return os << rang::fg::reset; };
auto Red = [](std::ostream& os) -> std::ostream& { return os << rang::fg::red; };
auto Green = [](std::ostream& os) -> std::ostream& { return os << rang::fg::green; };
auto Yellow = [](std::ostream& os) -> std::ostream& { return os << rang::fg::yellow; };
auto Blue = [](std::ostream& os) -> std::ostream& { return os << rang::fg::blue; };
auto Magenta = [](std::ostream& os) -> std::ostream& { return os << rang::fg::magenta; };
auto Cyan = [](std::ostream& os) -> std::ostream& { return os << rang::fg::cyan; };
auto Gray = [](std::ostream& os) -> std::ostream& { return os << rang::fg::gray; };

bool DEBUG = true;
namespace fs = std::filesystem;


VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
//class C_InteractorStyle;

double bkR = 0.129, bkG = 0.129, bkB = 0.141,
stR = 0.7, stG = 0.5, stB = 0.3;


class C_InteractorStyle : public vtkInteractorStyleTrackballCamera {
public:

	static C_InteractorStyle* New() {
		return new C_InteractorStyle;
	}
	vtkTypeMacro(C_InteractorStyle, vtkInteractorStyleTrackballCamera);

	virtual void OnMiddleButtonDown() override {
		this->StartPan();
	}

	virtual void OnMiddleButtonUp() override {
		this->EndPan();
	}

	virtual void OnRightButtonDown() override {
		this->StartRotate();
	}

	virtual void OnRightButtonUp() override {
		this->EndRotate();
	}

	virtual void Pan() override {
		if (this->CurrentRenderer == nullptr || this->Interactor == nullptr) return;
		vtkRenderWindowInteractor* rwi = this->Interactor;
		vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
		if (!camera) return;
		int* lastPos = rwi->GetLastEventPosition();
		int* newPos = rwi->GetEventPosition();
		double dx = newPos[0] - lastPos[0];
		double dy = newPos[1] - lastPos[1];
		double scale = 0.05; // Adjust this scale to control the sensitivity of panning
		dx *= scale;
		dy *= scale;
		double right[3], up[3];
		camera->GetViewUp(up);
		this->CurrentRenderer->GetActiveCamera()->OrthogonalizeViewUp();
		vtkMath::Cross(camera->GetDirectionOfProjection(), up, right);
		vtkMath::Normalize(right);
		double cameraPosition[3], cameraFocalPoint[3];
		camera->GetPosition(cameraPosition);
		camera->GetFocalPoint(cameraFocalPoint);
		for (int i = 0; i < 3; i++) {
			cameraPosition[i] += dx * right[i] + dy * up[i];
			cameraFocalPoint[i] += dx * right[i] + dy * up[i];
		}
		camera->SetPosition(cameraPosition);
		camera->SetFocalPoint(cameraFocalPoint);
		this->CurrentRenderer->ResetCameraClippingRange();
		rwi->Render();
	}

	virtual void Dolly(double amount) override {
		if (this->CurrentRenderer == nullptr || this->Interactor == nullptr) return;
		vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
		if (!camera) return;
		const double DollyScaleFactor = 0.1, minDis = 50, maxDis = 300;
		amount = 1.0 + (amount - 1.0) * DollyScaleFactor;

		double newDistance = (camera->GetDistance()) * (1.0 / amount);  // Adjust the interpretation of amount

		if (newDistance < minDis) camera->SetDistance(minDis);
		else if (newDistance > maxDis) camera->SetDistance(maxDis);
		else camera->Dolly(amount);

		//if (DEBUG) std::cout << Yellow << "      Current Zoom Level: " << ColorEnd << camera->GetDistance() << std::endl;
		this->CurrentRenderer->ResetCameraClippingRange();
		this->Interactor->Render();
	}


	void SetFileNames(const std::vector<std::string>& files) {
		this->FileNames = files;
	}

	virtual void OnKeyPress() override {
		vtkRenderWindowInteractor* rwi = this->Interactor;
		std::string key = rwi->GetKeySym();

		std::cout << "        " << Yellow << "Pressed: " << ColorEnd << key << std::endl;

		if (key == "Up") {
			CurrentIndex = (CurrentIndex + 1) % FileNames.size();
			this->LoadModel(CurrentIndex);
			this->setup();
		}
		else if (key == "Down") {
			CurrentIndex = (CurrentIndex > 0) ? CurrentIndex - 1 : FileNames.size() - 1;
			this->LoadModel(CurrentIndex);
			this->setup();
		}
		else if (key == "Enter" || key == "Return") {
			this->LoadAllModel();
			this->setup();
		}
		else {
			// Call parent class keypress if not Up/Down arrow key
			vtkInteractorStyleTrackballCamera::OnKeyPress();
			return;
		}
		rwi->Render();
	}

	void LoadModel(size_t index) {
		if (!this->Renderer) {
			std::cerr << "Renderer is not set!" << std::endl;
			return;
		}
		fs::path STLfilepath(FileNames[index] + ".stl");
		fs::path PTSfilepath(FileNames[index] + ".pts");
		this->Renderer->RemoveAllViewProps();

		vtkNew<vtkSTLReader> reader;
		reader->SetFileName(STLfilepath.string().c_str());
		reader->Update();
		vtkNew<vtkPolyDataMapper> Modelmapper;
		Modelmapper->SetInputConnection(reader->GetOutputPort());
		vtkNew<vtkActor> Modelactor;
		Modelactor->SetMapper(Modelmapper);
		double bounds[6];
		Modelactor->GetBounds(bounds);
		double Height_ = bounds[5] - bounds[4];

		Modelactor->GetProperty()->SetDiffuse(0.9);
		Modelactor->GetProperty()->SetDiffuseColor(stR, stG, stB); // White color
		Modelactor->GetProperty()->SetSpecular(0.3); // Reduced specular reflection
		Modelactor->GetProperty()->SetSpecularPower(70.0); // Adjust specular power
		Modelactor->GetProperty()->SetInterpolationToPhong();
		this->Renderer->AddActor(Modelactor);
		//this->Renderer->ResetCamera();

		// Model Name Text
		vtkNew<vtkTextActor> ModelTextActor;
		ModelTextActor->GetTextProperty()->SetFontSize(20);
		ModelTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // 0.7, 0.5, 0.3
		ModelTextActor->SetDisplayPosition(10, 10);
		ModelTextActor->SetInput(("ID: " + STLfilepath.filename().string()).c_str());
		this->Renderer->AddActor(ModelTextActor);

		// Model Height Text
		vtkNew<vtkTextActor> HeightTextActor;
		HeightTextActor->GetTextProperty()->SetFontSize(20);
		HeightTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // 0.7, 0.5, 0.3
		HeightTextActor->SetDisplayPosition(10, 35);
		HeightTextActor->SetInput(("Height: " + std::format("{:.2f}", Height_) + " mm").c_str());
		this->Renderer->AddActor(HeightTextActor);


		// Read PTS file
		if (fs::exists(PTSfilepath)) {
			std::ifstream ptsFile(PTSfilepath);
			vtkSmartPointer<vtkPoints> PTSpoints = vtkSmartPointer<vtkPoints>::New();
			std::string line;
			while (std::getline(ptsFile, line)) {
				std::istringstream iss(line);
				double x, y, z;
				iss >> x >> y >> z;
				PTSpoints->InsertNextPoint(x, y, z);
			}

			vtkSmartPointer<vtkPolyData> PTSpointPolyData = vtkSmartPointer<vtkPolyData>::New();
			PTSpointPolyData->SetPoints(PTSpoints);
			vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
			vertexFilter->SetInputData(PTSpointPolyData);
			vertexFilter->Update();
			vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
			sphereSource->SetRadius(0.75);  // Set the radius of the spheres
			vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
			glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
			glyph3D->SetInputConnection(vertexFilter->GetOutputPort());
			glyph3D->ScalingOff();  // Turn off scaling to keep spheres uniform
			vtkSmartPointer<vtkPolyDataMapper> glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			glyphMapper->SetInputConnection(glyph3D->GetOutputPort());
			vtkSmartPointer<vtkActor> glyphActor = vtkSmartPointer<vtkActor>::New();
			glyphActor->SetMapper(glyphMapper);
			this->Renderer->AddActor(glyphActor);
			std::cout << "     >> " << Yellow << "Loaded model: " << ColorEnd << STLfilepath.filename().string()
				<< " & " << PTSfilepath.filename().string() << std::endl;
		}
		else {
			std::cout << "     >> " << Yellow << "Loaded model: " << ColorEnd << STLfilepath.filename().string() << std::endl;
		}
	}

	void LoadAllModel() {
		if (!this->Renderer) {
			std::cerr << "Renderer is not set!" << std::endl;
			return;
		}
		this->Renderer->RemoveAllViewProps();
		for (const std::string& FileName : FileNames) {
			fs::path filepath(FileName + ".stl");
			vtkNew<vtkSTLReader> reader;
			reader->SetFileName(filepath.string().c_str());
			reader->Update();
			vtkNew<vtkPolyDataMapper> mapper;
			mapper->SetInputConnection(reader->GetOutputPort());
			vtkNew<vtkActor> actor;
			actor->SetMapper(mapper);

			actor->GetProperty()->SetDiffuse(0.9);
			actor->GetProperty()->SetDiffuseColor(stR, stG, stB); // White color
			actor->GetProperty()->SetSpecular(0.3); // Reduced specular reflection
			actor->GetProperty()->SetSpecularPower(70.0); // Adjust specular power
			actor->GetProperty()->SetInterpolationToPhong();
			this->Renderer->AddActor(actor);

			// Model Name Text
			vtkNew<vtkTextActor> ModelTextActor;
			ModelTextActor->GetTextProperty()->SetFontSize(20);
			ModelTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // 0.7, 0.5, 0.3
			ModelTextActor->SetDisplayPosition(10, 10);
			ModelTextActor->SetInput("All Models");
			this->Renderer->AddActor(ModelTextActor);

			// Model Height Text
			vtkNew<vtkTextActor> HeightTextActor;
			HeightTextActor->GetTextProperty()->SetFontSize(20);
			HeightTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // 0.7, 0.5, 0.3
			HeightTextActor->SetDisplayPosition(10, 35);
			HeightTextActor->SetInput("");
			this->Renderer->AddActor(HeightTextActor);

			std::cout << "     >> " << Yellow << "Loaded model: " << ColorEnd << filepath.filename().string() << std::endl;
		}
	}

	void setup() {
		if (!this->Renderer) {
			std::cerr << "Renderer is not set!" << std::endl;
			return;
		}
		// Center XYZ axes
		vtkNew<vtkAxesActor> XYZaxes;
		XYZaxes->SetTotalLength(6, 6, 6);
		XYZaxes->SetXAxisLabelText("");
		XYZaxes->SetYAxisLabelText("");
		XYZaxes->SetZAxisLabelText("");
		this->Renderer->AddActor(XYZaxes);

		// Center sphere
		vtkNew<vtkSphereSource> sphereSource;
		sphereSource->SetCenter(0.0, 0.0, 0.0);
		sphereSource->SetRadius(0.5);
		sphereSource->Update();
		vtkNew<vtkPolyDataMapper> sphereMapper;
		sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
		vtkNew<vtkActor> sphereActor;
		sphereActor->SetMapper(sphereMapper);
		sphereActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
		this->Renderer->AddActor(sphereActor);

		// Create a center plane
		vtkNew<vtkPlaneSource> planeSource;
		planeSource->SetXResolution(10);
		planeSource->SetYResolution(10);
		planeSource->SetOrigin(-50, -50, 0);
		planeSource->SetPoint1(50, -50, 0);
		planeSource->SetPoint2(-50, 50, 0);
		vtkNew<vtkPolyDataMapper> planeMapper;
		planeMapper->SetInputConnection(planeSource->GetOutputPort());
		vtkNew<vtkActor> planeActor;
		planeActor->SetMapper(planeMapper);
		planeActor->GetProperty()->SetRepresentationToWireframe();
		planeActor->GetProperty()->SetColor(0.3, 0.3, 0.3);
		this->Renderer->AddActor(planeActor);
	}

	void SetRenderer(vtkRenderer* renderer) {
		this->Renderer = renderer;
	}

private:
	std::vector<std::string> FileNames;
	size_t CurrentIndex = 0;
	vtkRenderer* Renderer = nullptr;
};
//vtkStandardNewMacro(C_InteractorStyle);



int main(int argc, char* argv[]) {
	// Get the current directory path
	std::string folderPath = fs::current_path().string();
	std::vector<std::string> STLFiles;

	// List all STL files in the current directory
	for (const auto& entry : fs::directory_iterator(folderPath)) {
		std::string extension = entry.path().extension().string();
		std::transform(extension.begin(), extension.end(), extension.begin(),
			[](unsigned char c) { return std::tolower(c); });

		if (extension == ".stl") {
			//STLFiles.push_back(entry.path().string());
			fs::path entryPath = entry.path();
			STLFiles.push_back(entryPath.replace_extension().string());
		}
	}


	if (STLFiles.empty()) {
		std::cout << std::endl;
		std::cout << Red << "\n\n        No STL files in input folder\n" << ColorEnd << std::endl;
		std::cout << "        Press enter to continue...";
		std::cin.get();
		return EXIT_SUCCESS;
	}

	if (DEBUG) std::cout << Yellow << "        Preparing Mesh Viewer." << ColorEnd << std::endl;


	vtkNew<vtkRenderer> Renderer;
	Renderer->SetBackground(bkR, bkG, bkB);

	// Window and interactor setup
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(700, 600);
	renderWindow->SetWindowName("Viewer tool");
	renderWindow->AddRenderer(Renderer);

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Custom Interaction Style
	vtkNew<C_InteractorStyle> style;
	style->SetDefaultRenderer(Renderer);
	style->SetRenderer(Renderer);
	style->SetFileNames(STLFiles);
	style->LoadAllModel();
	style->setup();
	renderWindowInteractor->SetInteractorStyle(style);

	// Camera setups
	vtkNew<vtkCameraOrientationWidget> camOrientManipulator;
	camOrientManipulator->SetParentRenderer(Renderer);
	camOrientManipulator->On();

	Renderer->ResetCamera();
	Renderer->GetActiveCamera()->SetPosition(0, 0, 160);
	Renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
	Renderer->GetActiveCamera()->SetViewUp(0, 1, 0);
	Renderer->ResetCameraClippingRange();

	renderWindow->Render();

	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();
	if (DEBUG) std::cout << Yellow << "        Viewer Exited." << ColorEnd << std::endl;

	return EXIT_SUCCESS;
}