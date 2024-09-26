#include <algorithm>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <map>
#include <regex>
#include <vector>

#include <vtkActor.h>
#include <vtkAssembly.h>
#include <vtkAutoInit.h>
#include <vtkAxesActor.h>
#include <vtkButtonWidget.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCameraOrientationWidget.h>
#include <vtkCaptionActor2D.h>
#include <vtkCommand.h>
#include <vtkConeSource.h>
#include <vtkGlyph3D.h>
#include <vtkImageLuminance.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkJPEGWriter.h>
#include <vtkLine.h>
#include <vtkMassProperties.h>
#include <vtkNew.h>
#include <vtkPlaneSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSTLReader.h>
#include <vtkSphereSource.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTexturedButtonRepresentation2D.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkWin32OpenGLRenderWindow.h>
#include <vtkWin32RenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkVectorText.h>
#include <vtkFollower.h>
#include <vtkTransform.h>
#include <vtkLinearExtrusionFilter.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <Eigen/Dense>
#include <rang.hpp>
#include <opencv2/opencv.hpp>
#include <FreeImage.h>
#include <minizip/unzip.h>
#include <hpdf.h>


using namespace std;
using namespace filesystem;
using namespace cv;
using namespace rang;

#define END		fg::reset
#define RED		fg::red	
#define GREEN	fg::green	
#define YELLOW	fg::yellow
#define BLUE	fg::blue	
#define MAGENTA fg::magenta
#define CYAN	fg::cyan	
#define GRAY	fg::gray	
#define BOLD	style::bold
#define NORM	style::reset


regex filenameRegex("([0-9]+)([a-zA-Z]{2,3})([0-9]*)");
regex SVGRegex("<line[^>]*x1=\"([0-9.\\-]+)\"[^>]*y1=\"([0-9.\\-]+)\"[^>]*x2=\"([0-9.\\-]+)\"[^>]*y2=\"([0-9.\\-]+)\"[^>]*text=\"([^\"]+)\"[^>]*height_start=\"([0-9.\\-]+)\"[^>]*height_end=\"([0-9.\\-]+)\"");
regex PTSRegex(R"(\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)(?:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))?(?:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))?(?:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))?)");
regex NCRegex(R"(N\d+\s*B(-?\d*\.?\d*)\s*C(-?\d*\.?\d*)\s*X(-?\d*\.?\d*)\s*Y(-?\d*\.?\d*)\s*Z(-?\d*\.?\d*))");


struct C_Actors {
	vtkSmartPointer<vtkProp> Actor;  // Can be vtkActor or vtkAssembly
	vtkSmartPointer<vtkTextActor> textActor;
};

map<string, C_Actors> actorCache;

#ifdef _MSC_VER
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
#endif


typedef CGAL::Simple_cartesian<double>							Cgal_Kernel;
typedef CGAL::Surface_mesh<Cgal_Kernel::Point_3>				Cgal_Mesh;
typedef Cgal_Kernel::Point_3									Cgal_Point;
typedef Cgal_Mesh::Vertex_index									Cgal_Vertex;

#define M_PI				3.14159265358979323846

bool DEBUG = false;

double
main_BK[3] = { 0.129, 0.129, 0.141 },
main_STLCOLOR[3] = { 0.7, 0.5, 0.3 },
HOLDERCOLOR[3] = { 0.6, 0.6, 0.6 },
AllSTLCOLOR[3] = { 0.8, 0.8, 0.9 },
main_PTSCOLOR[3] = { 0.9, 0.9, 0.9 },
main_NCCOLOR[3] = { 0.196, 0.51, 0.965 };

double svgTextHeight = 1.5;
double svgScale[3] = { 1.5, 1.5, 1.0 };

#include "STL_files.h"
#include "export_pdf.h"
#include "export_mp4.h"

class C_InteractorStyle : public vtkInteractorStyleTrackballCamera {

private:
	bool loadSTL = true, loadHolder = true, ValidPTSnormals = false;
	size_t loadNC = 1, loadPTS = 1, CurrentIndex = 0;
	double NCsphereR = 0.75, threshold = 0.5;
	vector<string> FileNamesPath;
	vtkNew<vtkPoints> PTSpoints;
	vtkNew<vtkPoints> NCpoints;
	vtkSmartPointer<vtkSTLReader> ModelReader = vtkSmartPointer<vtkSTLReader>::New();


public:
	static C_InteractorStyle* New();
	vtkTypeMacro(C_InteractorStyle, vtkInteractorStyleTrackballCamera);

	virtual void OnLeftButtonDown() override {
		this->StartPan();
	}

	virtual void OnLeftButtonUp() override {
		this->EndPan();
	}

	virtual void OnRightButtonDown() override {
		this->StartRotate();
	}

	virtual void OnRightButtonUp() override {
		this->EndRotate();
	}

	virtual void OnMiddleButtonDown() override {
		this->GetDefaultRenderer()->ResetCamera();
		this->GetDefaultRenderer()->GetActiveCamera()->SetPosition(0, 0, 160);
		this->GetDefaultRenderer()->GetActiveCamera()->SetFocalPoint(0, 0, 0);
		this->GetDefaultRenderer()->GetActiveCamera()->SetViewUp(0, 1, 0);
		this->GetDefaultRenderer()->ResetCameraClippingRange();
		this->GetInteractor()->Render();
	}
	virtual void OnMiddleButtonUp() override {}

	//virtual void Pan() override {
	//	if (this->GetDefaultRenderer() == nullptr || this->GetInteractor() == nullptr) return;
	//	//vtkRenderWindowInteractor* rwi = this->GetInteractor();
	//	vtkCamera* camera = this->GetDefaultRenderer()->GetActiveCamera();
	//	if (!camera) return;
	//	int* lastPos = this->GetInteractor()->GetLastEventPosition();
	//	int* newPos = this->GetInteractor()->GetEventPosition();
	//	double dx = newPos[0] - lastPos[0];
	//	double dy = newPos[1] - lastPos[1];
	//	double scale = 0.1; // Adjust this scale to control the sensitivity of panning
	//	dx *= scale;
	//	dy *= scale;
	//	double right[3], up[3];
	//	camera->GetViewUp(up);
	//	this->GetDefaultRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
	//	vtkMath::Cross(camera->GetDirectionOfProjection(), up, right);
	//	vtkMath::Normalize(right);
	//	double cameraPosition[3], cameraFocalPoint[3];
	//	camera->GetPosition(cameraPosition);
	//	camera->GetFocalPoint(cameraFocalPoint);
	//	for (int i = 0; i < 3; i++) {
	//		cameraPosition[i] += dx * right[i] + dy * up[i];
	//		cameraFocalPoint[i] += dx * right[i] + dy * up[i];
	//	}
	//	camera->SetPosition(cameraPosition);
	//	camera->SetFocalPoint(cameraFocalPoint);
	//	this->GetDefaultRenderer()->ResetCameraClippingRange();
	//	this->GetInteractor()->Render();
	//}

	virtual void Dolly(double amount) override {
		if (this->GetDefaultRenderer() == nullptr || this->GetInteractor() == nullptr) return;
		vtkCamera* camera = this->GetDefaultRenderer()->GetActiveCamera();
		if (!camera) return;
		const double DollyScaleFactor = 0.3, minDis = 30, maxDis = 600;
		amount = 1.0 + (amount - 1.0) * DollyScaleFactor;

		double newDistance = (camera->GetDistance()) * (1.0 / amount);  // Adjust the interpretation of amount

		if (newDistance < minDis) camera->SetDistance(minDis);
		else if (newDistance > maxDis) camera->SetDistance(maxDis);
		else camera->Dolly(amount);

		this->GetDefaultRenderer()->ResetCameraClippingRange();
		this->GetInteractor()->Render();
	}

	void SetFileNamesPath(const vector<string>& files) {
		this->FileNamesPath = files;
	}

	virtual void OnKeyPress() override {
		string key = this->GetInteractor()->GetKeySym();

		if (DEBUG) cout << "        " << YELLOW << "Pressed: " << END << key << endl;

		if (key == "Up" || key == "Down" || key == "Right" || key == "Left" || key == "h" || key == "H" ||
			key == "t" || key == "T" || key == "n" || key == "N" || key == "m" || key == "M") {

			this->GetDefaultRenderer()->RemoveAllViewProps();

			if (key == "Up") {
				CurrentIndex = (CurrentIndex + 1) % FileNamesPath.size();
				NCsphereR = 0.75;
			}
			else if (key == "Down") {
				CurrentIndex = (CurrentIndex > 0) ? CurrentIndex - 1 : FileNamesPath.size() - 1;
				NCsphereR = 0.75;
			}
			else if (key == "Right") {
				NCsphereR = min(NCsphereR + 0.05, 1.0);
			}
			else if (key == "Left") {
				NCsphereR = max(NCsphereR - 0.05, 0.1);
			}
			else if (key == "m" || key == "M") {
				loadSTL = !loadSTL;
			}
			else if (key == "t" || key == "T") {
				loadPTS = loadPTS + 1;
				if ((ValidPTSnormals && loadPTS >= 3) || (!ValidPTSnormals && loadPTS >= 2)) {
					loadPTS = 0;
					PTSpoints->Reset();
				}
			}
			else if (key == "n" || key == "N") {
				loadNC = loadNC + 1;
				if (loadNC >= 3) {
					loadNC = 0;
					NCpoints->Reset();
				}
			}
			else if (key == "h" || key == "H") {
				loadHolder = !loadHolder;
			}

			if (loadHolder) this->LoadFixtureHolder();
			if (loadSTL) this->LoadSTLfile(CurrentIndex, ModelReader);
			if (loadPTS >= 1) this->LoadPTSModel(CurrentIndex, 0.65, loadPTS, ModelReader);
			if (loadNC >= 1) this->LoadNCfile(CurrentIndex, NCsphereR, loadNC);
			this->LoadLaserMarkingZip(CurrentIndex);
			this->ShowFilesCount(CurrentIndex);
			this->setup();
		}
		else if (key == "Enter" || key == "Return") {
			this->GetDefaultRenderer()->RemoveAllViewProps();
			this->LoadAllSTLfiles();
			this->setup();

		}
		else if (key == "Escape") {
			this->GetInteractor()->GetRenderWindow()->Finalize();
			this->GetInteractor()->TerminateApp();
			return;
		}
		else {
			vtkInteractorStyleTrackballCamera::OnKeyPress();
			return;
		}
		this->GetInteractor()->Render();
	}

	void ShowFilesCount(size_t index) {
		// Files Count Text
		vtkNew<vtkTextActor> FilesCountTextActor;
		FilesCountTextActor->GetTextProperty()->SetFontSize(23);
		FilesCountTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
		FilesCountTextActor->SetDisplayPosition(10, 310);
		FilesCountTextActor->SetInput(std::format("Files: {}/{}", index + 1, FileNamesPath.size()).c_str());
		this->GetDefaultRenderer()->AddActor(FilesCountTextActor);
	}

	void LoadFixtureHolder() {
		Cgal_Mesh mesh;
		//istringstream iss(string(reinterpret_cast<const char*>(_3pins_holder_Data), sizeof(_3pins_holder_Data)), ios::binary);
		istringstream iss(string(reinterpret_cast<const char*>(_House_Holder_Data), sizeof(_House_Holder_Data)), ios::binary);
		if (CGAL::IO::read_STL(iss, mesh)) {
			vtkNew<vtkPoints> VTKpoints;
			vtkNew<vtkCellArray> polygons;
			map<Cgal_Vertex, vtkIdType> vertexIdMap;
			vtkIdType vtkId = 0;

			for (auto v : mesh.vertices()) {
				const Cgal_Point& p = mesh.point(v);
				VTKpoints->InsertNextPoint(p.x(), p.y(), p.z());
				vertexIdMap[v] = vtkId++;
			}

			for (auto f : mesh.faces()) {
				vtkNew<vtkIdList> polygon;
				CGAL::Vertex_around_face_iterator<Cgal_Mesh> vbegin, vend;
				boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(f), mesh);
				for (; vbegin != vend; ++vbegin) {
					polygon->InsertNextId(vertexIdMap[*vbegin]);
				}
				polygons->InsertNextCell(polygon);
			}

			vtkNew<vtkPolyData> polyData;
			polyData->SetPoints(VTKpoints);
			polyData->SetPolys(polygons);
			vtkNew<vtkPolyDataMapper> staticMapper;
			staticMapper->SetInputData(polyData);
			vtkNew<vtkActor> staticActor;
			staticActor->SetMapper(staticMapper);
			staticActor->GetProperty()->SetColor(HOLDERCOLOR);
			this->GetDefaultRenderer()->AddActor(staticActor);
		}
	}

	void LoadLaserMarkingZip(size_t index) {
		string ZipfileName = FileNamesPath[index];
		string svgFileName = path(ZipfileName + ".svg").filename().string();
		smatch casematch;
		string CaseID;

		if (regex_search(ZipfileName, casematch, filenameRegex)) {
			CaseID = casematch[1];
		}
		else {
			cerr << "Error: Could not extract the zip prefix from filename: " << ZipfileName << endl;
			return;
		}

		path LaserMarkingZipPath(CaseID + ".zip");

		if (!exists(LaserMarkingZipPath)) return;

		unzFile zipfile = unzOpen(LaserMarkingZipPath.string().c_str());
		if (zipfile == nullptr) {
			cerr << "Cannot open zip file: " << LaserMarkingZipPath << endl;
			return;
		}

		if (unzLocateFile(zipfile, svgFileName.c_str(), 0) != UNZ_OK) {
			cerr << "Error: SVG file not found in the zip: " << svgFileName << endl;
			unzClose(zipfile);
			return;
		}

		if (unzOpenCurrentFile(zipfile) != UNZ_OK) {
			cerr << "Error: Cannot open the SVG file inside the zip: " << svgFileName << endl;
			unzClose(zipfile);
			return;
		}

		char buffer[4096];
		string svgContent;
		int bytesRead;
		while ((bytesRead = unzReadCurrentFile(zipfile, buffer, sizeof(buffer))) > 0) {
			svgContent.append(buffer, bytesRead);
		}

		if (bytesRead < 0) {
			cerr << "Error: Failed to read the SVG file inside the zip" << endl;
			unzCloseCurrentFile(zipfile);
			unzClose(zipfile);
			return;
		}

		unzCloseCurrentFile(zipfile);
		unzClose(zipfile);

		smatch svgmatch;
		if (regex_search(svgContent, svgmatch, SVGRegex)) {
			double x1 = stod(svgmatch[1]);
			double y1 = stod(svgmatch[2]);
			double z1 = stod(svgmatch[6]);

			double x2 = stod(svgmatch[3]);
			double y2 = stod(svgmatch[4]);
			double z2 = stod(svgmatch[7]);

			string text = svgmatch[5];

			vtkNew<vtkPoints> SVGpoints;
			SVGpoints->InsertNextPoint(x1, y1, z1);
			SVGpoints->InsertNextPoint(x2, y2, z2);

			if (SVGpoints->GetNumberOfPoints() > 0) {
				vtkNew<vtkPolyData> SVGpointPolyData;
				SVGpointPolyData->SetPoints(SVGpoints);
				vtkNew<vtkVertexGlyphFilter> SVGvertexFilter;
				SVGvertexFilter->SetInputData(SVGpointPolyData);
				SVGvertexFilter->Update();
				vtkNew<vtkSphereSource> SVGsphereSource;
				SVGsphereSource->SetRadius(0.5);
				vtkNew<vtkGlyph3D> SVGglyph3D;
				SVGglyph3D->SetSourceConnection(SVGsphereSource->GetOutputPort());
				SVGglyph3D->SetInputConnection(SVGvertexFilter->GetOutputPort());
				SVGglyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> SVGglyphMapper;
				SVGglyphMapper->SetInputConnection(SVGglyph3D->GetOutputPort());
				vtkNew<vtkActor> SVGglyphActor;
				SVGglyphActor->SetMapper(SVGglyphMapper);
				SVGglyphActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
				this->GetDefaultRenderer()->AddActor(SVGglyphActor);


				double midpoint[3] = { (x1 + x2) / 2.0, (y1 + y2) / 2.0, (z1 + z2) / 2.0 };
				double direction[3] = { x2 - x1, y2 - y1, z2 - z1 };
				double length = sqrt(pow(direction[0], 2) + pow(direction[1], 2) + pow(direction[2], 2));
				direction[0] /= length;
				direction[1] /= length;
				direction[2] /= length;

				vtkNew<vtkVectorText> vectorText;
				vectorText->SetText(text.c_str());
				vtkNew<vtkLinearExtrusionFilter> extrudeFilter;
				extrudeFilter->SetInputConnection(vectorText->GetOutputPort());
				extrudeFilter->SetExtrusionTypeToNormalExtrusion();
				extrudeFilter->SetScaleFactor(svgTextHeight);  // Adjust thickness here
				extrudeFilter->Update();
				vtkNew<vtkPolyDataMapper> textMapper;
				textMapper->SetInputConnection(extrudeFilter->GetOutputPort());
				vtkNew<vtkActor> textActor;
				textActor->SetMapper(textMapper);
				textActor->SetScale(svgScale);
				textActor->GetProperty()->SetColor(0.0, 1.0, 0.0);

				vtkNew<vtkTransform> transform;
				transform->Identity();
				transform->Translate(x1, y1, z1 + 0.5);
				transform->RotateZ(atan2(direction[1], direction[0]) * 180.0 / vtkMath::Pi());
				transform->RotateY((acos(direction[2]) * 180.0 / vtkMath::Pi()) - 90.0);
				transform->RotateX(15);
				textActor->SetUserTransform(transform);
				this->GetDefaultRenderer()->AddActor(textActor);

				// Model SVG Text
				vtkNew<vtkTextActor> ZIPTextActor;
				ZIPTextActor->GetTextProperty()->SetFontSize(23);
				ZIPTextActor->GetTextProperty()->SetColor(0.0, 1.0, 0.0);
				ZIPTextActor->SetPosition(10, 280); // 20, 100
				ZIPTextActor->GetTextProperty()->SetFontFamilyToArial();
				ZIPTextActor->GetTextProperty()->SetFontSize(23);
				ZIPTextActor->SetInput(std::format("Laser Marking: {}", LaserMarkingZipPath.filename().string()).c_str());
				this->GetDefaultRenderer()->AddActor(ZIPTextActor);
			}
		}
		else {
			cerr << "Error: Could not find the expected <line> tag with the required attributes." << endl;
		}
	}

	void LoadNCfile(size_t index, double NCsphereR, size_t loadNCNormals) {
		if (!this->GetDefaultRenderer()) {
			cerr << "Renderer is not set!" << endl;
			return;
		}

		path NCfilepath(FileNamesPath[index] + ".nc");
		if (exists(NCfilepath)) {
			vtkNew<vtkPoints> NCnormalspoints;
			vtkNew<vtkCellArray> NCnormals;
			int NCcount = 0;
			double ToolOffset = -92.0,
				ToolNormalsLength = 4.0;
			ifstream NCfile(NCfilepath);
			string NCline;

			//ofstream outTXTFile("Log.txt");
			//if (!outTXTFile) {
			//	cerr << "Unable to open file for writing";
			//	return;
			//}

			NCpoints->Reset();
			while (getline(NCfile, NCline)) {
				smatch NCmatch;
				if (regex_search(NCline, NCmatch, NCRegex)) {
					double B = stod(NCmatch[1].str());
					double C = stod(NCmatch[2].str());
					double X = stod(NCmatch[3].str());
					double Y = stod(NCmatch[4].str());
					double Z = stod(NCmatch[5].str());

					double B_rad = vtkMath::RadiansFromDegrees(B);
					double C_rad = vtkMath::RadiansFromDegrees(C);

					Eigen::Matrix4d RB, RC, Txyz, TFinal, Tn;

					RB <<
						cos(B_rad), 0, sin(B_rad), 0,
						0, 1, 0, 0,
						-sin(B_rad), 0, cos(B_rad), 0,
						0, 0, 0, 1;

					RC <<
						cos(C_rad), -sin(C_rad), 0, 0,
						sin(C_rad), cos(C_rad), 0, 0,
						0, 0, 1, 0,
						0, 0, 0, 1;

					Txyz <<
						1, 0, 0, X,
						0, 1, 0, Y,
						0, 0, 1, Z,
						0, 0, 0, 1;

					TFinal = RC * Txyz * RB;
					Tn = RC * RB;
					Eigen::Vector4d toolTipPos = TFinal * Eigen::Vector4d(0, 0, ToolOffset, 1);

					NCpoints->InsertNextPoint(toolTipPos(0), toolTipPos(1), toolTipPos(2));

					NCcount++;
					Eigen::Vector3d TNormal = Tn.block<3, 1>(0, 2);

					//outTXTFile << setw(3) << setfill('0') << NCcount
					//	<< " - G  " << B << ", " << C << ", " << X << ", " << Y << ", " << Z << "   T "
					//	<< toolTipPos(0) << ", " << toolTipPos(1) << ", " << toolTipPos(2) << "   N "
					//	<< TNormal(0) << ", " << TNormal(1) << ", " << TNormal(2) << endl;

					if (loadNCNormals == 2) {
						Eigen::Vector3d toolNormal = Tn.block<3, 1>(0, 2) * ToolNormalsLength;
						vtkIdType Pid[2];
						Pid[0] = NCnormalspoints->InsertNextPoint(toolTipPos(0), toolTipPos(1), toolTipPos(2));
						double endpoint[3] = { toolTipPos(0) + toolNormal[0], toolTipPos(1) + toolNormal[1], toolTipPos(2) + toolNormal[2] };
						Pid[1] = NCnormalspoints->InsertNextPoint(endpoint);
						vtkNew<vtkLine> Normalline;
						Normalline->GetPointIds()->SetId(0, Pid[0]);
						Normalline->GetPointIds()->SetId(1, Pid[1]);
						NCnormals->InsertNextCell(Normalline);
					}
				}
			}

			//outTXTFile.close();

			if (NCpoints->GetNumberOfPoints() > 0) {
				path PTSfilepath(FileNamesPath[index] + ".pts");
				if (exists(PTSfilepath) && PTSpoints->GetNumberOfPoints() > 0) {
					const double tolerance = 0.2; // Tolerance for alignment
					double totalAverageDiff = 0;
					double finalAverageTolerance = 0;
					double MaxDiff = 0;
					bool pointsAreAligned = true;
					int count = 0;

					for (vtkIdType i = 0; i < NCpoints->GetNumberOfPoints(); ++i) {
						double* ncPoint = NCpoints->GetPoint(i);
						double minDistance = numeric_limits<double>::max();

						for (vtkIdType j = 0; j < PTSpoints->GetNumberOfPoints(); ++j) {
							double* ptsPoint = PTSpoints->GetPoint(j);
							double distance = 0.0;

							for (int k = 0; k < 3; ++k) {
								distance += pow(ncPoint[k] - ptsPoint[k], 2);
							}
							distance = sqrt(distance); // Euclidean distance

							if (distance < minDistance) minDistance = distance;
						}

						if (minDistance > tolerance) pointsAreAligned = false;

						totalAverageDiff += minDistance;
						MaxDiff = max(MaxDiff, minDistance);
						count++;
					}

					if (count > 0) finalAverageTolerance = totalAverageDiff / count;

					vtkNew<vtkTextActor> AlignmentTextActor;
					AlignmentTextActor->GetTextProperty()->SetFontSize(23);
					AlignmentTextActor->SetDisplayPosition(10, 165);

					if (pointsAreAligned) {
						double newColor[3] = { 0.4, 0.4, 0.9 };
						memcpy(main_NCCOLOR, newColor, 3 * sizeof(double));

						AlignmentTextActor->SetInput(std::format("Good - Aligned\nPTS ({:d})\nNC ({:d})\nA{:.2f} : M{:.2f} : T{:.2f}mm",
							PTSpoints->GetNumberOfPoints(), NCpoints->GetNumberOfPoints(),
							finalAverageTolerance, MaxDiff, tolerance).c_str());

						if (DEBUG) {
							cout << "        PTS(" << PTSpoints->GetNumberOfPoints() << ") & NC("
								<< NCpoints->GetNumberOfPoints() << ") are aligned with average difference of ("
								<< finalAverageTolerance << ") within tolerance (" << tolerance << ")" << endl;
						}
					}
					else {
						double newColor[3] = { 1.0, 0.2, 0.3 };
						memcpy(main_NCCOLOR, newColor, 3 * sizeof(double));

						AlignmentTextActor->SetInput(std::format("Bad - Not Aligned\nPTS ({:d})\nNC ({:d})\nA{:.2f} : M{:.2f} : T{:.2f}mm",
							PTSpoints->GetNumberOfPoints(), NCpoints->GetNumberOfPoints(),
							finalAverageTolerance, MaxDiff, tolerance).c_str());

						if (DEBUG) {
							cout << "        PTS(" << PTSpoints->GetNumberOfPoints() << ") & NC("
								<< NCpoints->GetNumberOfPoints() << ") are not aligned with average difference of ("
								<< finalAverageTolerance << ")" << endl;
						}
					}

					AlignmentTextActor->GetTextProperty()->SetColor(main_NCCOLOR);
					this->GetDefaultRenderer()->AddActor(AlignmentTextActor);
				}

				// Model NC Text
				vtkNew<vtkTextActor> NCTextActor;
				NCTextActor->GetTextProperty()->SetFontSize(23);
				NCTextActor->GetTextProperty()->SetColor(main_NCCOLOR);
				NCTextActor->SetDisplayPosition(10, 100);
				NCTextActor->SetInput(std::format("Tool Diameter: {:.2f} mm\nNC: {} - {} {}",
					NCsphereR * 2.0, NCfilepath.filename().string(), NCpoints->GetNumberOfPoints(),
					loadNCNormals == 2 ? ", Normals" : "").c_str());
				this->GetDefaultRenderer()->AddActor(NCTextActor);

				// NC Tool tip location
				vtkNew<vtkPolyData> NCpointPolyData;
				NCpointPolyData->SetPoints(NCpoints);
				vtkNew<vtkVertexGlyphFilter> NCvertexFilter;
				NCvertexFilter->SetInputData(NCpointPolyData);
				NCvertexFilter->Update();
				vtkNew<vtkSphereSource> NCsphereSource;
				NCsphereSource->SetRadius(NCsphereR);  // Set the radius of the spheres
				vtkNew<vtkGlyph3D> NCglyph3D;
				NCglyph3D->SetSourceConnection(NCsphereSource->GetOutputPort());
				NCglyph3D->SetInputConnection(NCvertexFilter->GetOutputPort());
				NCglyph3D->ScalingOff();  // Turn off scaling to keep spheres uniform
				vtkNew<vtkPolyDataMapper> NCglyphMapper;
				NCglyphMapper->SetInputConnection(NCglyph3D->GetOutputPort());
				vtkNew<vtkActor> NCglyphActor;
				NCglyphActor->SetMapper(NCglyphMapper);
				NCglyphActor->GetProperty()->SetColor(main_NCCOLOR);
				this->GetDefaultRenderer()->AddActor(NCglyphActor);

				// Visualize first point
				double firstPoint[3], secondPoint[3], dir[3], length;
				NCpoints->GetPoint(0, firstPoint);
				NCpoints->GetPoint(1, secondPoint);
				dir[0] = secondPoint[0] - firstPoint[0];
				dir[1] = secondPoint[1] - firstPoint[1];
				dir[2] = secondPoint[2] - firstPoint[2];
				length = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
				dir[0] /= length;
				dir[1] /= length;
				dir[2] /= length;

				vtkNew<vtkConeSource> firstConeSource;
				firstConeSource->SetCenter(firstPoint);
				firstConeSource->SetRadius(1.0);
				firstConeSource->SetHeight(2.0);
				firstConeSource->SetResolution(20);
				firstConeSource->SetDirection(dir);  // Pointing towards the second point
				vtkNew<vtkPolyDataMapper> firstMapper;
				firstMapper->SetInputConnection(firstConeSource->GetOutputPort());
				vtkNew<vtkActor> firstActor;
				firstActor->SetMapper(firstMapper);
				firstActor->GetProperty()->SetColor(0.2, 0.9, 0.2);
				this->GetDefaultRenderer()->AddActor(firstActor);

				if (loadNCNormals == 2) {
					vtkNew<vtkPolyData> NCnormalsPolyData;
					NCnormalsPolyData->SetPoints(NCnormalspoints);
					NCnormalsPolyData->SetLines(NCnormals);
					vtkNew<vtkPolyDataMapper> NCnormalsMapper;
					NCnormalsMapper->SetInputData(NCnormalsPolyData);
					vtkNew<vtkActor> NCnormalsActor;
					NCnormalsActor->SetMapper(NCnormalsMapper);
					NCnormalsActor->GetProperty()->SetLineWidth(2.0);
					NCnormalsActor->GetProperty()->SetColor(main_NCCOLOR);
					this->GetDefaultRenderer()->AddActor(NCnormalsActor);
				}

				if (DEBUG) cout << "     >> " << YELLOW << "Loaded NC: " << END << NCfilepath.filename().string() << endl;
			}
		}
	}

	void LoadPTSModel(size_t index, double NCsphereR, size_t loadPTSNormals, vtkSTLReader* ModelReader) {
		if (!this->GetDefaultRenderer()) {
			cerr << "Renderer is not set!" << endl;
			return;
		}

		path PTSfilepath(FileNamesPath[index] + ".pts");
		if (exists(PTSfilepath)) {
			ifstream ptsFile(PTSfilepath);
			string PTSline;
			vector<string> PTSlines;
			vtkNew<vtkPoints> PTSnormalspoints;
			vtkNew<vtkCellArray> PTSnormals;
			double PTSNormalsLength = 4.0;

			PTSpoints->Reset();
			while (getline(ptsFile, PTSline)) {
				PTSlines.push_back(PTSline);
			}

			// Insert points and normals in reverse order
			for (auto it = PTSlines.rbegin(); it != PTSlines.rend(); ++it) {
				smatch PTSmatch;
				if (regex_match(*it, PTSmatch, PTSRegex)) {
					if (PTSmatch.size() >= 4) {
						double x = stod(PTSmatch[1].str());
						double y = stod(PTSmatch[2].str());
						double z = stod(PTSmatch[3].str());
						PTSpoints->InsertNextPoint(x, y, z);

						// Check if normals are present
						if (PTSmatch.size() == 7 && PTSmatch[4].matched && PTSmatch[5].matched && PTSmatch[6].matched) {
							double u = stod(PTSmatch[4].str());
							double v = stod(PTSmatch[5].str());
							double w = stod(PTSmatch[6].str());
							Eigen::Vector3d PTSNorm(u, v, w);

							PTSNorm.normalize();
							PTSNorm *= PTSNormalsLength;

							vtkIdType Pid[2];
							Pid[0] = PTSnormalspoints->InsertNextPoint(x, y, z);
							Pid[1] = PTSnormalspoints->InsertNextPoint(x + PTSNorm[0], y + PTSNorm[1], z + PTSNorm[2]);
							vtkNew<vtkLine> Normalline;
							Normalline->GetPointIds()->SetId(0, Pid[0]);
							Normalline->GetPointIds()->SetId(1, Pid[1]);
							PTSnormals->InsertNextCell(Normalline);

							ValidPTSnormals = true;
						}
						else {
							ValidPTSnormals = false;
						}
					}
				}
			}

			if (PTSpoints->GetNumberOfPoints() > 0) {
				vtkNew<vtkPoints> inLoopPoints;
				vtkNew<vtkPoints> outLoopPoints;

				for (vtkIdType i = 1; i < PTSpoints->GetNumberOfPoints() - 1; ++i) {
					double point[3];
					double prevPoint[3];
					double nextPoint[3];
					PTSpoints->GetPoint(i, point);
					PTSpoints->GetPoint(i - 1, prevPoint);
					PTSpoints->GetPoint(i + 1, nextPoint);

					// check if a point is far from its neighbors
					double distanceToPrev = sqrt(
						pow(point[0] - prevPoint[0], 2) +
						pow(point[1] - prevPoint[1], 2) +
						pow(point[2] - prevPoint[2], 2));

					double distanceToNext = sqrt(
						pow(point[0] - nextPoint[0], 2) +
						pow(point[1] - nextPoint[1], 2) +
						pow(point[2] - nextPoint[2], 2));

					if (!(distanceToPrev > threshold) || !(distanceToNext > threshold)) {
						inLoopPoints->InsertNextPoint(point);
					}
					else {
						outLoopPoints->InsertNextPoint(point);
					}
				}

				if (ValidPTSnormals && loadPTSNormals == 2) {
					vtkNew<vtkPolyData> PTSnormalsPolyData;
					PTSnormalsPolyData->SetPoints(PTSnormalspoints);
					PTSnormalsPolyData->SetLines(PTSnormals);
					vtkNew<vtkPolyDataMapper> PTSnormalsMapper;
					PTSnormalsMapper->SetInputData(PTSnormalsPolyData);
					vtkNew<vtkActor> PTSnormalsActor;
					PTSnormalsActor->SetMapper(PTSnormalsMapper);
					PTSnormalsActor->GetProperty()->SetLineWidth(2.0);
					PTSnormalsActor->GetProperty()->SetColor(main_PTSCOLOR);
					this->GetDefaultRenderer()->AddActor(PTSnormalsActor);
				}

				// Load the STL file if available
				path STLfilepath(FileNamesPath[index] + ".stl");
				if (exists(STLfilepath) && ModelReader && ModelReader->GetOutput() != nullptr) {

					vtkNew<vtkImplicitPolyDataDistance> implicitDistance;
					implicitDistance->SetInput(ModelReader->GetOutput());

					vtkNew<vtkPoints> surfacePoints;
					vtkNew<vtkPoints> nonSurfacePoints;

					for (vtkIdType i = 0; i < inLoopPoints->GetNumberOfPoints(); ++i) {
						double point[3];
						inLoopPoints->GetPoint(i, point);
						double distance = implicitDistance->EvaluateFunction(point);
						if (abs(distance) <= 0.5)
							surfacePoints->InsertNextPoint(point);
						else
							nonSurfacePoints->InsertNextPoint(point);
					}

					// Visualize in-loop non-surface points (red color)
					vtkNew<vtkPolyData> nonSurfacePolyData;
					nonSurfacePolyData->SetPoints(nonSurfacePoints);
					vtkNew<vtkVertexGlyphFilter> nonSurfaceVertexFilter;
					nonSurfaceVertexFilter->SetInputData(nonSurfacePolyData);
					nonSurfaceVertexFilter->Update();
					vtkNew<vtkSphereSource> nonSurfaceSphereSource;
					nonSurfaceSphereSource->SetRadius(NCsphereR);
					vtkNew<vtkGlyph3D> nonSurfaceGlyph3D;
					nonSurfaceGlyph3D->SetSourceConnection(nonSurfaceSphereSource->GetOutputPort());
					nonSurfaceGlyph3D->SetInputConnection(nonSurfaceVertexFilter->GetOutputPort());
					nonSurfaceGlyph3D->ScalingOff();
					vtkNew<vtkPolyDataMapper> nonSurfaceGlyphMapper;
					nonSurfaceGlyphMapper->SetInputConnection(nonSurfaceGlyph3D->GetOutputPort());
					vtkNew<vtkActor> nonSurfaceGlyphActor;
					nonSurfaceGlyphActor->SetMapper(nonSurfaceGlyphMapper);
					nonSurfaceGlyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  // REDcolor
					this->GetDefaultRenderer()->AddActor(nonSurfaceGlyphActor);

					// Visualize in-loop surface points (original color)
					vtkNew<vtkPolyData> surfacePolyData;
					surfacePolyData->SetPoints(surfacePoints);
					vtkNew<vtkVertexGlyphFilter> surfaceVertexFilter;
					surfaceVertexFilter->SetInputData(surfacePolyData);
					surfaceVertexFilter->Update();
					vtkNew<vtkSphereSource> surfaceSphereSource;
					surfaceSphereSource->SetRadius(NCsphereR);
					vtkNew<vtkGlyph3D> surfaceGlyph3D;
					surfaceGlyph3D->SetSourceConnection(surfaceSphereSource->GetOutputPort());
					surfaceGlyph3D->SetInputConnection(surfaceVertexFilter->GetOutputPort());
					surfaceGlyph3D->ScalingOff();
					vtkNew<vtkPolyDataMapper> surfaceGlyphMapper;
					surfaceGlyphMapper->SetInputConnection(surfaceGlyph3D->GetOutputPort());
					vtkNew<vtkActor> surfaceGlyphActor;
					surfaceGlyphActor->SetMapper(surfaceGlyphMapper);
					surfaceGlyphActor->GetProperty()->SetColor(main_PTSCOLOR);  // Original color
					this->GetDefaultRenderer()->AddActor(surfaceGlyphActor);
				}
				else {
					// Visualize in-loop points
					vtkNew<vtkPolyData> inLoopPolyData;
					inLoopPolyData->SetPoints(inLoopPoints);
					vtkNew<vtkVertexGlyphFilter> inLoopVertexFilter;
					inLoopVertexFilter->SetInputData(inLoopPolyData);
					inLoopVertexFilter->Update();
					vtkNew<vtkSphereSource> inLoopSphereSource;
					inLoopSphereSource->SetRadius(NCsphereR);
					vtkNew<vtkGlyph3D> inLoopGlyph3D;
					inLoopGlyph3D->SetSourceConnection(inLoopSphereSource->GetOutputPort());
					inLoopGlyph3D->SetInputConnection(inLoopVertexFilter->GetOutputPort());
					inLoopGlyph3D->ScalingOff();
					vtkNew<vtkPolyDataMapper> inLoopGlyphMapper;
					inLoopGlyphMapper->SetInputConnection(inLoopGlyph3D->GetOutputPort());
					vtkNew<vtkActor> inLoopGlyphActor;
					inLoopGlyphActor->SetMapper(inLoopGlyphMapper);
					inLoopGlyphActor->GetProperty()->SetColor(main_PTSCOLOR);
					this->GetDefaultRenderer()->AddActor(inLoopGlyphActor);
				}

				// Visualize out-of-loop points
				vtkNew<vtkPolyData> outLoopPolyData;
				outLoopPolyData->SetPoints(outLoopPoints);
				vtkNew<vtkVertexGlyphFilter> outLoopVertexFilter;
				outLoopVertexFilter->SetInputData(outLoopPolyData);
				outLoopVertexFilter->Update();
				vtkNew<vtkSphereSource> outLoopSphereSource;
				outLoopSphereSource->SetRadius(NCsphereR);
				vtkNew<vtkGlyph3D> outLoopGlyph3D;
				outLoopGlyph3D->SetSourceConnection(outLoopSphereSource->GetOutputPort());
				outLoopGlyph3D->SetInputConnection(outLoopVertexFilter->GetOutputPort());
				outLoopGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> outLoopGlyphMapper;
				outLoopGlyphMapper->SetInputConnection(outLoopGlyph3D->GetOutputPort());
				vtkNew<vtkActor> outLoopGlyphActor;
				outLoopGlyphActor->SetMapper(outLoopGlyphMapper);
				outLoopGlyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
				this->GetDefaultRenderer()->AddActor(outLoopGlyphActor);

				// Visualize first point
				double firstPoint[3], secondPoint[3], dir[3], length;
				PTSpoints->GetPoint(0, firstPoint);
				PTSpoints->GetPoint(1, secondPoint);
				dir[0] = secondPoint[0] - firstPoint[0];
				dir[1] = secondPoint[1] - firstPoint[1];
				dir[2] = secondPoint[2] - firstPoint[2];
				length = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
				dir[0] /= length;
				dir[1] /= length;
				dir[2] /= length;

				vtkNew<vtkConeSource> firstConeSource;
				firstConeSource->SetCenter(firstPoint);
				firstConeSource->SetRadius(1.0);
				firstConeSource->SetHeight(2.0);
				firstConeSource->SetResolution(20);
				firstConeSource->SetDirection(dir);  // Pointing towards the second point
				vtkNew<vtkPolyDataMapper> firstMapper;
				firstMapper->SetInputConnection(firstConeSource->GetOutputPort());
				vtkNew<vtkActor> firstActor;
				firstActor->SetMapper(firstMapper);
				firstActor->GetProperty()->SetColor(0.2, 0.8, 0.2);
				this->GetDefaultRenderer()->AddActor(firstActor);

				// Model PTS Text
				vtkNew<vtkTextActor> PTSTextActor;
				PTSTextActor->GetTextProperty()->SetFontSize(23);
				PTSTextActor->GetTextProperty()->SetColor(main_PTSCOLOR);
				PTSTextActor->SetDisplayPosition(10, 77);
				PTSTextActor->SetInput(std::format("PTS{}: {} {} - {}", ValidPTSnormals ? "N" : "", PTSfilepath.filename().string(), loadPTSNormals == 2 ? "- Normals Enabled" : "", PTSpoints->GetNumberOfPoints()).c_str());
				this->GetDefaultRenderer()->AddActor(PTSTextActor);

				if (DEBUG) cout << "     >> " << YELLOW << "Loaded PTS: " << END << PTSfilepath.filename().string() << endl;
				if (DEBUG) cout << "        Initial size of PTS Points: " << PTSpoints->GetNumberOfPoints() << endl;
			}
		}
	}

	void LoadSTLfile(size_t index, vtkSmartPointer<vtkSTLReader>& ModelReader) {
		if (!this->GetDefaultRenderer()) {
			cerr << "Renderer is not set!" << endl;
			return;
		}
		path STLfilepath(FileNamesPath[index] + ".stl");
		if (!exists(STLfilepath)) {
			ModelReader = nullptr;
			return;
		}

		ModelReader->SetFileName(STLfilepath.string().c_str());
		ModelReader->Update();
		vtkNew<vtkPolyDataMapper> Modelmapper;
		Modelmapper->SetInputConnection(ModelReader->GetOutputPort());

		vtkNew<vtkActor> Modelactor;
		Modelactor->SetMapper(Modelmapper);
		double bounds[6];
		Modelactor->GetBounds(bounds);
		double Height_ = bounds[5] - bounds[4];

		vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
		massProperties->SetInputConnection(ModelReader->GetOutputPort());
		massProperties->Update();
		double volumeMilliliters = massProperties->GetVolume() / 1000.0;

		Modelactor->GetProperty()->SetDiffuse(0.9);
		Modelactor->GetProperty()->SetDiffuseColor(main_STLCOLOR);
		Modelactor->GetProperty()->SetSpecular(0.3); // Reduced specular reflection
		Modelactor->GetProperty()->SetSpecularPower(70.0); // Adjust specular power
		Modelactor->GetProperty()->SetInterpolationToPhong();
		this->GetDefaultRenderer()->AddActor(Modelactor);

		// Model Text
		vtkNew<vtkTextActor> ModelTextActor;
		ModelTextActor->GetTextProperty()->SetFontSize(23);
		ModelTextActor->GetTextProperty()->SetColor(main_STLCOLOR);
		ModelTextActor->SetDisplayPosition(10, 10);
		ModelTextActor->SetInput(std::format("ID: {}\nHeight: {:.2f} mm\nVolume: {:.2f} mL",
			STLfilepath.filename().string(), Height_, volumeMilliliters).c_str());
		this->GetDefaultRenderer()->AddActor(ModelTextActor);

		if (DEBUG) cout << "     >> " << YELLOW << "Loaded model: " << END << STLfilepath.filename().string() << endl;
	}

	void LoadAllSTLfiles() {
		if (!this->GetDefaultRenderer()) {
			cerr << "Renderer is not set!" << endl;
			return;
		}
		double totalVolumeMilliliters = 0.0;
		double maxHeight = 0.0;
		int ModelsCount = 0;

		for (const string& FileName : FileNamesPath) {
			path filepath(FileName + ".stl");
			vtkNew<vtkSTLReader> ModelsReader;
			ModelsReader->SetFileName(filepath.string().c_str());
			ModelsReader->Update();
			vtkNew<vtkPolyDataMapper> mapper;
			mapper->SetInputConnection(ModelsReader->GetOutputPort());
			vtkNew<vtkActor> actor;
			actor->SetMapper(mapper);

			actor->GetProperty()->SetDiffuse(0.9);
			actor->GetProperty()->SetDiffuseColor(AllSTLCOLOR); // White color
			actor->GetProperty()->SetSpecular(0.3); // Reduced specular reflection
			actor->GetProperty()->SetSpecularPower(70.0); // Adjust specular power
			actor->GetProperty()->SetInterpolationToPhong();
			this->GetDefaultRenderer()->AddActor(actor);

			double bounds[6];
			actor->GetBounds(bounds);
			double Height_ = bounds[5] - bounds[4];
			maxHeight = max(maxHeight, Height_);

			ModelsCount++;

			vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
			massProperties->SetInputConnection(ModelsReader->GetOutputPort());
			massProperties->Update();
			double volumeMilliliters = massProperties->GetVolume() / 1000.0;
			totalVolumeMilliliters += volumeMilliliters;

			if (DEBUG) cout << "     >> " << YELLOW << "Loaded model: " << END << filepath.filename().string() << endl;
		}

		// Model Count Text
		vtkNew<vtkTextActor> AllModelsTextActor;
		AllModelsTextActor->GetTextProperty()->SetFontSize(23);
		AllModelsTextActor->GetTextProperty()->SetColor(AllSTLCOLOR);
		AllModelsTextActor->SetDisplayPosition(10, 10);
		AllModelsTextActor->SetInput(std::format("Total Count: {:d} Models\nMax Height: {:.2f} mm\nTotal Volume: {:.2f} mL",
			ModelsCount, maxHeight, totalVolumeMilliliters).c_str());
		this->GetDefaultRenderer()->AddActor(AllModelsTextActor);
	}

	void setup() {
		if (!this->GetDefaultRenderer()) {
			cerr << "Renderer is not set!" << endl;
			return;
		}
		// Center XYZ axes
		vtkNew<vtkAxesActor> XYZaxes;
		XYZaxes->SetTotalLength(10, 10, 10);
		XYZaxes->SetXAxisLabelText("");
		XYZaxes->SetYAxisLabelText("");
		XYZaxes->SetZAxisLabelText("");
		this->GetDefaultRenderer()->AddActor(XYZaxes);

		// Center sphere
		vtkNew<vtkSphereSource> sphereSource;
		sphereSource->SetCenter(0.0, 0.0, 0.0);
		sphereSource->SetRadius(1.0);
		sphereSource->Update();
		vtkNew<vtkPolyDataMapper> sphereMapper;
		sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
		vtkNew<vtkActor> sphereActor;
		sphereActor->SetMapper(sphereMapper);
		//sphereActor->GetProperty()->SetColor(0.7, 0.5, 0.3);
		sphereActor->GetProperty()->SetColor(0.9, 0.9, 0.9);
		this->GetDefaultRenderer()->AddActor(sphereActor);

		// Create a center plane
		vtkNew<vtkPlaneSource> planeSource;
		planeSource->SetXResolution(20);
		planeSource->SetYResolution(20);
		planeSource->SetOrigin(-100, -100, 0);
		planeSource->SetPoint1(100, -100, 0);
		planeSource->SetPoint2(-100, 100, 0);
		vtkNew<vtkPolyDataMapper> planeMapper;
		planeMapper->SetInputConnection(planeSource->GetOutputPort());
		vtkNew<vtkActor> planeActor;
		planeActor->SetMapper(planeMapper);
		planeActor->GetProperty()->SetRepresentationToWireframe();
		planeActor->GetProperty()->SetColor(0.3, 0.3, 0.3);
		this->GetDefaultRenderer()->AddActor(planeActor);
	}
};
vtkStandardNewMacro(C_InteractorStyle);


struct ButtonCallbackData {
	vtkTextActor* textActor;
	string textData;
	double buttonMidX;
	double buttonMidY;
};

void CreateButtonImage(vtkImageData* image, int width, int height, const unsigned char* color) {
	image->SetDimensions(width, height, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			copy(color, color + 3, pixel);  // Copy color data into the pixel
		}
	}
}

void CenterTextActor(vtkTextActor* textActor, vtkRenderer* renderer, double buttonMidX, double buttonMidY) {
	double bbox[4];
	textActor->GetBoundingBox(renderer, bbox);
	double textWidth = static_cast<double>(bbox[1] - bbox[0]);
	double textHeight = static_cast<double>(bbox[3] - bbox[2]);
	double textPosX = buttonMidX - (textWidth / 2.0);
	double textPosY = buttonMidY - (textHeight / 2.0);
	textActor->SetPosition(textPosX, textPosY - 2);
}

void ResetTextCallback(vtkObject* caller, long unsigned int, void* clientData, void*) {
	ButtonCallbackData* data = static_cast<ButtonCallbackData*>(clientData);
	vtkTextActor* textActor = data->textActor;
	if (textActor) {
		textActor->SetInput(data->textData.c_str());
		vtkRenderWindowInteractor* interactor = static_cast<vtkRenderWindowInteractor*>(caller);
		vtkRenderer* renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

		CenterTextActor(textActor, renderer, data->buttonMidX, data->buttonMidY);

		interactor->GetRenderWindow()->Render();  // Re-render the window to update the text
	}
}

void SavePDF_ButtonCallback(vtkObject* caller, long unsigned int, void* clientData, void*) {
	ButtonCallbackData* data = static_cast<ButtonCallbackData*>(clientData);
	vtkTextActor* textActor = data->textActor;
	if (!textActor) return;

	const char* currentText = textActor->GetInput();
	if (save_PDF()) {
		textActor->SetInput("SAVED");
	}
	else {
		textActor->SetInput("FAILED");
	}

	vtkButtonWidget* buttonWidget = static_cast<vtkButtonWidget*>(caller);
	vtkTexturedButtonRepresentation2D* buttonRep = static_cast<vtkTexturedButtonRepresentation2D*>(buttonWidget->GetRepresentation());
	vtkRenderer* renderer = buttonRep->GetRenderer();

	ButtonCallbackData* callbackData = new ButtonCallbackData();
	callbackData->textActor = textActor;
	callbackData->textData = "PDF";
	callbackData->buttonMidX = data->buttonMidX;
	callbackData->buttonMidY = data->buttonMidY;
	vtkRenderWindowInteractor* interactor = buttonWidget->GetInteractor();
	interactor->CreateRepeatingTimer(5000);
	vtkNew<vtkCallbackCommand> resetTextCallback;
	resetTextCallback->SetCallback(ResetTextCallback);
	resetTextCallback->SetClientData(callbackData);
	interactor->AddObserver(vtkCommand::TimerEvent, resetTextCallback);

	// Use the passed buttonMidX and buttonMidY for centering the text
	CenterTextActor(textActor, renderer, data->buttonMidX, data->buttonMidY);
	renderer->GetRenderWindow()->Render();
}

void SaveMP4_ButtonCallback(vtkObject* caller, long unsigned int, void* clientData, void*) {
	ButtonCallbackData* data = static_cast<ButtonCallbackData*>(clientData);
	vtkTextActor* textActor = data->textActor;
	if (!textActor) return;

	const char* currentText = textActor->GetInput();
	if (SaveVideo()) {
		textActor->SetInput("SAVED");
	}
	else {
		textActor->SetInput("FAILED");
	}

	vtkButtonWidget* buttonWidget = static_cast<vtkButtonWidget*>(caller);
	vtkTexturedButtonRepresentation2D* buttonRep = static_cast<vtkTexturedButtonRepresentation2D*>(buttonWidget->GetRepresentation());
	vtkRenderer* renderer = buttonRep->GetRenderer();

	ButtonCallbackData* callbackData = new ButtonCallbackData();
	callbackData->textActor = textActor;
	callbackData->textData = "MP4/GIF";
	callbackData->buttonMidX = data->buttonMidX;
	callbackData->buttonMidY = data->buttonMidY;
	vtkRenderWindowInteractor* interactor = buttonWidget->GetInteractor();
	interactor->CreateRepeatingTimer(5000);
	vtkNew<vtkCallbackCommand> resetTextCallback;
	resetTextCallback->SetCallback(ResetTextCallback);
	resetTextCallback->SetClientData(callbackData);
	interactor->AddObserver(vtkCommand::TimerEvent, resetTextCallback);

	// Use the passed buttonMidX and buttonMidY for centering the text
	CenterTextActor(textActor, renderer, data->buttonMidX, data->buttonMidY);
	renderer->GetRenderWindow()->Render();
}


int main(int argc, char* argv[]) {
	// Get the current directory path
	string folderPath = current_path().string();
	set<string> uniqueFiles; // To store unique file paths without extension

	for (const auto& entry : directory_iterator(folderPath)) {
		if (entry.is_regular_file()) { // Check if it's a regular file
			string extension = entry.path().extension().string();
			transform(extension.begin(), extension.end(), extension.begin(),
				[](unsigned char c) { return tolower(c); });

			if (extension == ".stl" || extension == ".pts" || extension == ".nc") {
				path entryPath = entry.path();
				entryPath.replace_extension();
				uniqueFiles.insert(entryPath.string());
			}
		}
	}

	vector<string> Files(uniqueFiles.begin(), uniqueFiles.end());

	if (Files.empty()) {
		cout << endl;
		cout << RED << "\n\n        No STL or PTS or NC files in input folder\n" << END << endl;
		cout << "        Press enter to continue...";
		cin.get();
		return EXIT_SUCCESS;
	}

	if (DEBUG) cout << YELLOW << "        Preparing Mesh Viewer." << END << endl;

	cout << endl;
	cout << GREEN << "        Keyboard Shortcuts:" << END << endl;
	cout << YELLOW << "          DOWN/UP    key  ->>  Opens Next/Previous file" << endl;
	cout << "          ENTER      key  ->>  Opens All STL files together" << endl;
	cout << "          M/m        Key  ->>  Enable/Disable STL View" << endl;
	cout << "          T/t        Key  ->>  Enable/Disable PTS (PTSN) View" << endl;
	cout << "          N/n        Key  ->>  Enable/Disable NC (Normals) View" << endl;
	cout << "          H/h        Key  ->>  Enable/Disable Fixture Holder View" << endl;
	cout << "          LEFT/RIGHT key  ->>  Decrease/Increase Tool Diameter" << endl;
	cout << "          S/W        Key  ->>  Solid/WireFrame STL View" << END << endl;


	vtkNew<vtkRenderer> Renderer;
	Renderer->SetBackground(main_BK);
	Renderer->SetUseShadows(true);
	Renderer->SetUseDepthPeeling(1);        // Enable depth peeling
	Renderer->SetMaximumNumberOfPeels(100);
	Renderer->SetOcclusionRatio(0.7);

	// Window and interactor setup
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(1024, 768);
	renderWindow->SetWindowName("AB Viewer tool (STL, PTS, NC)");
	renderWindow->AddRenderer(Renderer);
	renderWindow->SetAlphaBitPlanes(1);  // Enable alpha channel if supported
	renderWindow->SetMultiSamples(0);    // Disable multi-sampling to see depth peeling

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Custom Interaction Style
	vtkNew<C_InteractorStyle> style;
	style->SetDefaultRenderer(Renderer);
	//style->SetRenderer(Renderer);
	style->SetFileNamesPath(Files);
	style->LoadFixtureHolder();
	vtkSmartPointer<vtkSTLReader> ModelReader = vtkSmartPointer<vtkSTLReader>::New();
	style->LoadSTLfile(0, ModelReader);
	style->LoadPTSModel(0, 0.65, 1, ModelReader);
	style->LoadNCfile(0, 0.75, 1);
	style->LoadLaserMarkingZip(0);
	style->ShowFilesCount(0);
	style->setup();
	renderWindowInteractor->SetInteractorStyle(style);


	// save pdf button
	unsigned char TextureColor[3] = { 179, 128, 77 };

	int width = 60;
	int height = 20;
	double buttonMidX = 500;
	double buttonMidY = 35;

	vtkNew<vtkImageData> texture;
	CreateButtonImage(texture, width, height, TextureColor);
	vtkNew<vtkTexturedButtonRepresentation2D> buttonRepresentation;
	buttonRepresentation->SetNumberOfStates(1);
	buttonRepresentation->SetButtonTexture(0, texture);
	vtkNew<vtkButtonWidget> buttonWidget;
	buttonWidget->SetInteractor(renderWindowInteractor);
	buttonWidget->SetRepresentation(buttonRepresentation);
	double bounds[6] = {
		buttonMidX - width, buttonMidX + width,  // Adjust width boundaries
		buttonMidY - height, buttonMidY + height,  // Adjust height boundaries
		0.0, 0.0
	};
	buttonRepresentation->SetPlaceFactor(1);
	buttonRepresentation->PlaceWidget(bounds);
	vtkNew<vtkTextActor> textActor;
	textActor->SetInput("PDF");
	textActor->GetTextProperty()->SetFontSize(20);
	textActor->GetTextProperty()->SetBold(1);
	textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
	CenterTextActor(textActor, Renderer, buttonMidX, buttonMidY);
	ButtonCallbackData* callbackData = new ButtonCallbackData();
	callbackData->textActor = textActor;
	callbackData->textData = "PDF";
	callbackData->buttonMidX = buttonMidX;
	callbackData->buttonMidY = buttonMidY;
	vtkNew<vtkCallbackCommand> SavePDFButtonCallback;
	SavePDFButtonCallback->SetCallback(SavePDF_ButtonCallback);
	SavePDFButtonCallback->SetClientData(callbackData);
	buttonWidget->AddObserver(vtkCommand::StateChangedEvent, SavePDFButtonCallback);
	buttonWidget->On();
	Renderer->AddActor2D(textActor);

	// save mp4 button
	double buttonMidX2 = 650;
	double buttonMidY2 = 35;

	vtkNew<vtkTexturedButtonRepresentation2D> buttonRepresentation2;
	buttonRepresentation2->SetNumberOfStates(1);
	buttonRepresentation2->SetButtonTexture(0, texture);
	vtkNew<vtkButtonWidget> buttonWidget2;
	buttonWidget2->SetInteractor(renderWindowInteractor);
	buttonWidget2->SetRepresentation(buttonRepresentation2);
	double bounds2[6] = {
		buttonMidX2 - width, buttonMidX2 + width,  // Adjust width boundaries
		buttonMidY2 - height, buttonMidY2 + height,  // Adjust height boundaries
		0.0, 0.0
	};
	buttonRepresentation2->SetPlaceFactor(1);
	buttonRepresentation2->PlaceWidget(bounds2);
	vtkNew<vtkTextActor> textActor2;
	textActor2->SetInput("MP4/GIF");
	textActor2->GetTextProperty()->SetFontSize(20);
	textActor2->GetTextProperty()->SetBold(1);
	textActor2->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
	CenterTextActor(textActor2, Renderer, buttonMidX2, buttonMidY2);
	ButtonCallbackData* callbackData2 = new ButtonCallbackData();
	callbackData2->textActor = textActor2;
	callbackData2->textData = "MP4/GIF";
	callbackData2->buttonMidX = buttonMidX2;
	callbackData2->buttonMidY = buttonMidY2;
	vtkNew<vtkCallbackCommand> SaveMP4ButtonCallback;
	SaveMP4ButtonCallback->SetCallback(SaveMP4_ButtonCallback);
	SaveMP4ButtonCallback->SetClientData(callbackData2);
	buttonWidget2->AddObserver(vtkCommand::StateChangedEvent, SaveMP4ButtonCallback);
	buttonWidget2->On();
	Renderer->AddActor2D(textActor2);

	//// observer to the renderer
	//vtkNew<RenderingObserver> renderingObserver;
	//Renderer->AddObserver(vtkCommand::StartEvent, renderingObserver);
	//Renderer->AddObserver(vtkCommand::EndEvent, renderingObserver);

	//Renderer->AutomaticLightCreationOff();
	//// Create a light
	//vtkNew<vtkLight> light;
	//light->SetLightTypeToSceneLight();
	//light->PositionalOn();
	//light->SetPosition(0, 0, 1);
	//light->SetConeAngle(10);
	//light->SetFocalPoint(0, 0, 0);
	//light->SetDiffuseColor(1.0, 0.0, 0.0);
	//light->SetAmbientColor(0.0, 1.0, 0.0);
	//light->SetSpecularColor(0.0, 0.0, 1.0);

	renderWindow->Render();

	//Renderer->AddLight(light);

	// Camera setups
	vtkNew<vtkCameraOrientationWidget> camOrientManipulator;
	camOrientManipulator->SetParentRenderer(Renderer);
	camOrientManipulator->On();

	Renderer->ResetCamera();
	Renderer->GetActiveCamera()->SetPosition(0, 0, 160);
	Renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
	Renderer->GetActiveCamera()->SetViewUp(0, 1, 0);
	//Renderer->GetActiveCamera()->ParallelProjectionOn();
	Renderer->ResetCameraClippingRange();

	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();
	if (DEBUG) cout << YELLOW << "        Viewer Exited." << END << endl;

	return EXIT_SUCCESS;
}