#include <future> 
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <set>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <tuple>
#include <Windows.h>

#include "rang.hpp"
#include "Holders_STL.h"

#include <vtkWin32RenderWindowInteractor.h>
#include <vtkWin32OpenGLRenderWindow.h>
#include <vtkLightCollection.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
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
#include <vtkSphereSource.h>
#include <vtkAxesActor.h>
#include <vtkCaptionActor2D.h>
#include <vtkTextProperty.h>
#include <vtkTextActor.h>
#include <vtkTextActor3D.h>
#include <vtkPlaneSource.h>
#include <vtkCameraOrientationWidget.h>
#include <vtkSTLReader.h>
#include <vtkGlyph3D.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkMassProperties.h>
#include <vtkLine.h>
#include <vtkCommand.h>

#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>


typedef CGAL::Simple_cartesian<double>							Cgal_Kernel;
typedef CGAL::Surface_mesh<Cgal_Kernel::Point_3>				Cgal_Mesh;
typedef Cgal_Kernel::Point_3									Cgal_Point;
typedef Cgal_Mesh::Vertex_index									Cgal_Vertex;

#define M_PI 3.14159265358979323846
#define Min(a,b)            (((a) < (b)) ? (a) : (b))
#define Max(a,b)            (((a) > (b)) ? (a) : (b))

#define deg2rad(degrees)  degrees * M_PI / 180.0


//<< Red << << ColorEnd <<
auto ColorEnd = [](std::ostream& os) -> std::ostream& { return os << rang::fg::reset; };
auto Red = [](std::ostream& os) -> std::ostream& { return os << rang::fg::red; };
auto Green = [](std::ostream& os) -> std::ostream& { return os << rang::fg::green; };
auto Yellow = [](std::ostream& os) -> std::ostream& { return os << rang::fg::yellow; };
auto Blue = [](std::ostream& os) -> std::ostream& { return os << rang::fg::blue; };
auto Magenta = [](std::ostream& os) -> std::ostream& { return os << rang::fg::magenta; };
auto Cyan = [](std::ostream& os) -> std::ostream& { return os << rang::fg::cyan; };
auto Gray = [](std::ostream& os) -> std::ostream& { return os << rang::fg::gray; };

bool DEBUG = false;
namespace fs = std::filesystem;


VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

double bkR = 0.129, bkG = 0.129, bkB = 0.141,
stR = 0.7, stG = 0.5, stB = 0.3,
fR = 0.6, fG = 0.6, fB = 0.6,
stRall = 0.8, stGall = 0.8, stBall = 0.9,
tR = 0.9, tG = 0.9, tB = 0.9,
nR = 0.4, nG = 0.4, nB = 0.9;


std::tuple<bool, double, double> arePointsAligned(double* pt1, double* pt2, double tolerance) {
	double totalDiff = 0.0;
	double Diff = 0.0;
	bool aligned = true;

	for (int i = 0; i < 3; ++i) {
		Diff = std::fabs(pt1[i] - pt2[i]);
		totalDiff += Diff;
		if (Diff > tolerance) {
			aligned = false;
		}
	}

	double averageDiff = totalDiff / 3;  // Calculate the average difference
	return { aligned, averageDiff, Diff };
}


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

	//virtual void Pan() override {
	//	if (this->CurrentRenderer == nullptr || this->Interactor == nullptr) return;
	//	vtkRenderWindowInteractor* rwi = this->Interactor;
	//	vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
	//	if (!camera) return;
	//	int* lastPos = rwi->GetLastEventPosition();
	//	int* newPos = rwi->GetEventPosition();
	//	double dx = newPos[0] - lastPos[0];
	//	double dy = newPos[1] - lastPos[1];
	//	double scale = 0.1; // Adjust this scale to control the sensitivity of panning
	//	dx *= scale;
	//	dy *= scale;
	//	double right[3], up[3];
	//	camera->GetViewUp(up);
	//	this->CurrentRenderer->GetActiveCamera()->OrthogonalizeViewUp();
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
	//	this->CurrentRenderer->ResetCameraClippingRange();
	//	rwi->Render();
	//}

	virtual void Dolly(double amount) override {
		if (this->CurrentRenderer == nullptr || this->Interactor == nullptr) return;
		vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
		if (!camera) return;
		const double DollyScaleFactor = 0.3, minDis = 50, maxDis = 600;
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

		if (DEBUG) std::cout << "        " << Yellow << "Pressed: " << ColorEnd << key << std::endl;

		if (key == "Up" || key == "Down" || key == "Right" || key == "Left" || key == "h" || key == "H" ||
			key == "t" || key == "T" || key == "n" || key == "N" || key == "m" || key == "M") {

			this->Renderer->RemoveAllViewProps();

			if (key == "Up") {
				CurrentIndex = (CurrentIndex + 1) % FileNames.size();
				NCsphereR = 0.75;
			}
			else if (key == "Down") {
				CurrentIndex = (CurrentIndex > 0) ? CurrentIndex - 1 : FileNames.size() - 1;
				NCsphereR = 0.75;
			}
			else if (key == "Right") {
				NCsphereR = Min(NCsphereR + 0.05, 1.0);
			}
			else if (key == "Left") {
				NCsphereR = Max(NCsphereR - 0.05, 0.1);
			}
			else if (key == "m" || key == "M") {
				loadSTL = !loadSTL;
			}
			else if (key == "t" || key == "T") {
				loadPTS = !loadPTS;
			}
			else if (key == "n" || key == "N") {
				loadNC = loadNC + 1;
				if (loadNC >= 3) loadNC = 0;
			}
			else if (key == "h" || key == "H") {
				loadHolder = !loadHolder;
			}

			if (loadHolder) LoadFixtureHolder();
			if (loadSTL) this->LoadSTLfile(CurrentIndex);
			if (loadPTS) this->LoadPTSModel(CurrentIndex);
			if (loadNC >= 1) this->LoadNCfile(CurrentIndex, NCsphereR, loadNC);
			this->setup();
		}
		else if (key == "Enter" || key == "Return") {
			this->Renderer->RemoveAllViewProps();
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
		rwi->Render();
	}


	void LoadFixtureHolder() {
		Cgal_Mesh mesh;
		//std::istringstream iss(std::string(reinterpret_cast<const char*>(_3pins_holder_Data), sizeof(_3pins_holder_Data)), std::ios::binary);
		std::istringstream iss(std::string(reinterpret_cast<const char*>(_House_Holder_Data), sizeof(_House_Holder_Data)), std::ios::binary);
		if (CGAL::IO::read_STL(iss, mesh)) {
			vtkNew<vtkPoints> VTKpoints;
			vtkNew<vtkCellArray> polygons;
			std::map<Cgal_Vertex, vtkIdType> vertexIdMap;
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
			staticActor->GetProperty()->SetColor(fR, fG, fB);
			this->Renderer->AddActor(staticActor);
		}
	}

	void LoadNCfile(size_t index, double NCsphereR, size_t loadNormals) {
		if (!this->Renderer) {
			std::cerr << "Renderer is not set!" << std::endl;
			return;
		}

		fs::path NCfilepath(FileNames[index] + ".nc");
		if (fs::exists(NCfilepath)) {
			vtkNew<vtkPoints> NCnormalspoints;
			vtkNew<vtkCellArray> NCnormals;

			double ToolOffset = -92.0, ToolNormalsLength = 5.0;
			std::ifstream NCfile(NCfilepath);
			std::string NCline;

			std::regex NCRegex(R"(N\d+\s*B(-?\d*\.?\d*)\s*C(-?\d*\.?\d*)\s*X(-?\d*\.?\d*)\s*Y(-?\d*\.?\d*)\s*Z(-?\d*\.?\d*))");

			NCpoints->Reset();
			while (std::getline(NCfile, NCline)) {
				std::smatch NCmatch;
				if (std::regex_search(NCline, NCmatch, NCRegex)) {
					double B = std::stod(NCmatch[1].str());
					double C = std::stod(NCmatch[2].str());
					double X = std::stod(NCmatch[3].str());
					double Y = std::stod(NCmatch[4].str());
					double Z = std::stod(NCmatch[5].str());

					double B_rad = deg2rad(B);
					double C_rad = deg2rad(C);

					Eigen::Matrix4d RB, RC, Txyz;

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

					Eigen::Matrix4d TFinal = RC * Txyz * RB;
					Eigen::Vector4d toolTipPos = TFinal * Eigen::Vector4d(0, 0, ToolOffset, 1);

					NCpoints->InsertNextPoint(toolTipPos(0), toolTipPos(1), toolTipPos(2));

					//std::cout << "X: " << toolTipPos(0) << ", Y: " << toolTipPos(1) << ", Z: " << toolTipPos(2) << std::endl;

					if (loadNormals == 2) {
						Eigen::Matrix4d Tn = RC * RB;
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

			if (NCpoints->GetNumberOfPoints() > 0) {
				const double tolerance = 0.2;
				double totalAverageDiff = 0;
				double finalAverageTolerance = 0;
				double MaxDiff = 0;
				bool pointsAreEqual = true;
				int count = 0;

				for (vtkIdType i = 0; i < PTSpoints->GetNumberOfPoints(); ++i) {
					double* ptsPoint = PTSpoints->GetPoint(i);
					double* ncPoint = NCpoints->GetPoint(i);

					auto [aligned, averageDiff, Diff] = arePointsAligned(ptsPoint, ncPoint, tolerance);
					totalAverageDiff += averageDiff;
					MaxDiff = Max(MaxDiff, Diff);
					count++;

					if (!aligned) {
						pointsAreEqual = false;
					}
				}

				if (count > 0) finalAverageTolerance = totalAverageDiff / count;

				// PTS and NC Alignment Text
				vtkNew<vtkTextActor> AlignmentTextActor;
				AlignmentTextActor->GetTextProperty()->SetFontSize(23);
				AlignmentTextActor->SetDisplayPosition(10, 145);
				if (pointsAreEqual) {
					if (DEBUG) std::cout << "        PTS(" << PTSpoints->GetNumberOfPoints() << ") & NC("
						<< NCpoints->GetNumberOfPoints() << ") are Aligned (" << finalAverageTolerance << ") within tolerance ("
						<< tolerance << ")" << std::endl;

					AlignmentTextActor->GetTextProperty()->SetColor(0.2, 1.0, 0.3);
					AlignmentTextActor->SetInput(std::format("Good - Aligned\nPTS ({:d})\nNC ({:d})\nA{:.2f} : M{:.2f} : T{:.2f}mm", PTSpoints->GetNumberOfPoints(),
						NCpoints->GetNumberOfPoints(), finalAverageTolerance, MaxDiff, tolerance).c_str());
					nR = 0.4, nG = 0.4, nB = 0.9;
				}
				else {
					if (DEBUG) std::cout << "        PTS(" << PTSpoints->GetNumberOfPoints() << ") & NC("
						<< NCpoints->GetNumberOfPoints() << ") are not Aligned (" << finalAverageTolerance << ")" << std::endl;

					AlignmentTextActor->GetTextProperty()->SetColor(1.0, 0.2, 0.3);
					AlignmentTextActor->SetInput(std::format("Bad - Not Aligned\nPTS ({:d})\nNC ({:d})\nA{:.2f} : M{:.2f} : T{:.2f}mm", PTSpoints->GetNumberOfPoints(),
						NCpoints->GetNumberOfPoints(), finalAverageTolerance, MaxDiff, tolerance).c_str());

					nR = 1.0, nG = 0.2, nB = 0.3;
				}
				this->Renderer->AddActor(AlignmentTextActor);

				// Model NC Text
				vtkNew<vtkTextActor> NCTextActor;
				NCTextActor->GetTextProperty()->SetFontSize(23);
				NCTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
				NCTextActor->SetDisplayPosition(10, 100);
				NCTextActor->SetInput(std::format("Tool Diameter: {:.2f} mm\nNC: {:s}",
					NCsphereR * 2.0, NCfilepath.filename().string()).c_str());
				this->Renderer->AddActor(NCTextActor);

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
				NCglyphActor->GetProperty()->SetColor(nR, nG, nB);
				this->Renderer->AddActor(NCglyphActor);

				if (loadNormals == 2) {
					// NC Tool normals
					vtkNew<vtkPolyData> NCnormalsPolyData;
					NCnormalsPolyData->SetPoints(NCnormalspoints);
					NCnormalsPolyData->SetLines(NCnormals);
					vtkNew<vtkPolyDataMapper> NCnormalsMapper;
					NCnormalsMapper->SetInputData(NCnormalsPolyData);
					vtkNew<vtkActor> NCnormalsActor;
					NCnormalsActor->SetMapper(NCnormalsMapper);
					NCnormalsActor->GetProperty()->SetColor(nR, nG, nB);
					this->Renderer->AddActor(NCnormalsActor);
				}

				if (DEBUG) std::cout << "     >> " << Yellow << "Loaded NC: " << ColorEnd << NCfilepath.filename().string() << std::endl;
			}
		}
	}

	// Function to check if a point is significantly far from its neighbors
	bool isPointOutOfLoop(double point[3], double prevPoint[3], double nextPoint[3], double threshold) {
		double distanceToPrev = std::sqrt(std::pow(point[0] - prevPoint[0], 2) +
			std::pow(point[1] - prevPoint[1], 2) +
			std::pow(point[2] - prevPoint[2], 2));

		double distanceToNext = std::sqrt(std::pow(point[0] - nextPoint[0], 2) +
			std::pow(point[1] - nextPoint[1], 2) +
			std::pow(point[2] - nextPoint[2], 2));

		return (distanceToPrev > threshold) || (distanceToNext > threshold);
	}

	void LoadPTSModel(size_t index) {
		if (!this->Renderer) {
			std::cerr << "Renderer is not set!" << std::endl;
			return;
		}

		fs::path PTSfilepath(FileNames[index] + ".pts");
		if (fs::exists(PTSfilepath)) {
			std::ifstream ptsFile(PTSfilepath);
			std::string PTSline;
			std::vector<std::string> PTSlines;
			std::regex PTSRegex(R"(\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))");

			PTSpoints->Reset();
			while (std::getline(ptsFile, PTSline)) {
				PTSlines.push_back(PTSline);
			}

			// Insert points in reverse order
			for (auto it = PTSlines.rbegin(); it != PTSlines.rend(); ++it) {
				std::smatch PTSmatch;
				if (std::regex_match(*it, PTSmatch, PTSRegex)) {
					if (PTSmatch.size() == 4) {
						double x = std::stod(PTSmatch[1].str());
						double y = std::stod(PTSmatch[2].str());
						double z = std::stod(PTSmatch[3].str());
						PTSpoints->InsertNextPoint(x, y, z);
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

					if (isPointOutOfLoop(point, prevPoint, nextPoint, threshold)) {
						outLoopPoints->InsertNextPoint(point);
					}
					else {
						inLoopPoints->InsertNextPoint(point);
					}
				}


				// Visualize in-loop points
				vtkNew<vtkPolyData> inLoopPolyData;
				inLoopPolyData->SetPoints(inLoopPoints);
				vtkNew<vtkVertexGlyphFilter> inLoopVertexFilter;
				inLoopVertexFilter->SetInputData(inLoopPolyData);
				inLoopVertexFilter->Update();
				vtkNew<vtkSphereSource> inLoopSphereSource;
				inLoopSphereSource->SetRadius(0.65);  // Set the radius of the spheres
				vtkNew<vtkGlyph3D> inLoopGlyph3D;
				inLoopGlyph3D->SetSourceConnection(inLoopSphereSource->GetOutputPort());
				inLoopGlyph3D->SetInputConnection(inLoopVertexFilter->GetOutputPort());
				inLoopGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> inLoopGlyphMapper;
				inLoopGlyphMapper->SetInputConnection(inLoopGlyph3D->GetOutputPort());
				vtkNew<vtkActor> inLoopGlyphActor;
				inLoopGlyphActor->SetMapper(inLoopGlyphMapper);
				inLoopGlyphActor->GetProperty()->SetColor(tR, tG, tB); // Blue color for in-loop points
				this->Renderer->AddActor(inLoopGlyphActor);

				// Visualize out-of-loop points
				vtkNew<vtkPolyData> outLoopPolyData;
				outLoopPolyData->SetPoints(outLoopPoints);
				vtkNew<vtkVertexGlyphFilter> outLoopVertexFilter;
				outLoopVertexFilter->SetInputData(outLoopPolyData);
				outLoopVertexFilter->Update();
				vtkNew<vtkSphereSource> outLoopSphereSource;
				outLoopSphereSource->SetRadius(0.65);  // Set the radius of the spheres
				vtkNew<vtkGlyph3D> outLoopGlyph3D;
				outLoopGlyph3D->SetSourceConnection(outLoopSphereSource->GetOutputPort());
				outLoopGlyph3D->SetInputConnection(outLoopVertexFilter->GetOutputPort());
				outLoopGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> outLoopGlyphMapper;
				outLoopGlyphMapper->SetInputConnection(outLoopGlyph3D->GetOutputPort());
				vtkNew<vtkActor> outLoopGlyphActor;
				outLoopGlyphActor->SetMapper(outLoopGlyphMapper);
				outLoopGlyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0); // Red color for out-of-loop points
				this->Renderer->AddActor(outLoopGlyphActor);
				
				// Visualize first point
				double firstPoint[3];
				PTSpoints->GetPoint(0, firstPoint);
				//vtkNew<vtkSphereSource> firstSphereSource;
				//firstSphereSource->SetCenter(firstPoint);
				//firstSphereSource->SetRadius(0.85);
				//vtkNew<vtkPolyDataMapper> firstMapper;
				//firstMapper->SetInputConnection(firstSphereSource->GetOutputPort());
				//vtkNew<vtkActor> firstActor;
				//firstActor->SetMapper(firstMapper);
				//firstActor->GetProperty()->SetColor(0.2, 0.2, 0.2); // Green color
				//this->Renderer->AddActor(firstActor);
				//vtkNew<vtkTextActor3D> textActor;
				//textActor->GetTextProperty()->SetFontSize(50);
				//textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); 
				//textActor->SetPosition(firstPoint[0], firstPoint[1], firstPoint[2]);
				//textActor->SetInput("1");
				//this->Renderer->AddActor(textActor);
				//
				//// Visualize last point
				//double lastPoint[3];
				//PTSpoints->GetPoint(PTSpoints->GetNumberOfPoints() - 1, lastPoint);
				//vtkNew<vtkSphereSource> lastSphereSource;
				//lastSphereSource->SetCenter(lastPoint);
				//lastSphereSource->SetRadius(0.85);
				//vtkNew<vtkPolyDataMapper> lastMapper;
				//lastMapper->SetInputConnection(lastSphereSource->GetOutputPort());
				//vtkNew<vtkActor> lastActor;
				//lastActor->SetMapper(lastMapper);
				//lastActor->GetProperty()->SetColor(0.2, 0.2, 0.2);
				//this->Renderer->AddActor(lastActor);

				// Model PTS Text
				vtkNew<vtkTextActor> PTSTextActor;
				PTSTextActor->GetTextProperty()->SetFontSize(23);
				PTSTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
				PTSTextActor->SetDisplayPosition(10, 77);
				PTSTextActor->SetInput(std::format("PTS: {:s}", PTSfilepath.filename().string()).c_str());
				this->Renderer->AddActor(PTSTextActor);

				if (DEBUG) std::cout << "     >> " << Yellow << "Loaded PTS: " << ColorEnd << PTSfilepath.filename().string() << std::endl;
				if (DEBUG) std::cout << "        Initial size of PTS Points: " << PTSpoints->GetNumberOfPoints() << std::endl;
			}
		}
	}

	void LoadSTLfile(size_t index) {
		if (!this->Renderer) {
			std::cerr << "Renderer is not set!" << std::endl;
			return;
		}
		fs::path STLfilepath(FileNames[index] + ".stl");
		if (fs::exists(STLfilepath)) {
			vtkNew<vtkSTLReader> ModelReader;
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
			Modelactor->GetProperty()->SetDiffuseColor(stR, stG, stB);
			Modelactor->GetProperty()->SetSpecular(0.3); // Reduced specular reflection
			Modelactor->GetProperty()->SetSpecularPower(70.0); // Adjust specular power
			Modelactor->GetProperty()->SetInterpolationToPhong();
			this->Renderer->AddActor(Modelactor);

			// Model Text
			vtkNew<vtkTextActor> ModelTextActor;
			ModelTextActor->GetTextProperty()->SetFontSize(23);
			ModelTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
			ModelTextActor->SetDisplayPosition(10, 10);
			ModelTextActor->SetInput(std::format("ID: {:s}\nHeight: {:.2f} mm\nVolume: {:.2f} mL",
				STLfilepath.filename().string(), Height_, volumeMilliliters).c_str());
			this->Renderer->AddActor(ModelTextActor);

			if (DEBUG) std::cout << "     >> " << Yellow << "Loaded model: " << ColorEnd << STLfilepath.filename().string() << std::endl;
		}
	}

	void LoadAllSTLfiles() {
		if (!this->Renderer) {
			std::cerr << "Renderer is not set!" << std::endl;
			return;
		}
		double totalVolumeMilliliters = 0.0;
		double maxHeight = 0.0;
		int ModelsCount = 0;

		for (const std::string& FileName : FileNames) {
			fs::path filepath(FileName + ".stl");
			vtkNew<vtkSTLReader> ModelsReader;
			ModelsReader->SetFileName(filepath.string().c_str());
			ModelsReader->Update();
			vtkNew<vtkPolyDataMapper> mapper;
			mapper->SetInputConnection(ModelsReader->GetOutputPort());
			vtkNew<vtkActor> actor;
			actor->SetMapper(mapper);

			actor->GetProperty()->SetDiffuse(0.9);
			actor->GetProperty()->SetDiffuseColor(stRall, stGall, stBall); // White color
			actor->GetProperty()->SetSpecular(0.3); // Reduced specular reflection
			actor->GetProperty()->SetSpecularPower(70.0); // Adjust specular power
			actor->GetProperty()->SetInterpolationToPhong();
			this->Renderer->AddActor(actor);

			double bounds[6];
			actor->GetBounds(bounds);
			double Height_ = bounds[5] - bounds[4];
			maxHeight = Max(maxHeight, Height_);

			ModelsCount++;

			vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
			massProperties->SetInputConnection(ModelsReader->GetOutputPort());
			massProperties->Update();
			double volumeMilliliters = massProperties->GetVolume() / 1000.0;
			totalVolumeMilliliters += volumeMilliliters;

			if (DEBUG) std::cout << "     >> " << Yellow << "Loaded model: " << ColorEnd << filepath.filename().string() << std::endl;
		}

		// Model Count Text
		vtkNew<vtkTextActor> AllModelsTextActor;
		AllModelsTextActor->GetTextProperty()->SetFontSize(23);
		AllModelsTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
		AllModelsTextActor->SetDisplayPosition(10, 10);
		AllModelsTextActor->SetInput(std::format("Total Count: {:d} Models\nMax Height: {:.2f} mm\nTotal Volume: {:.2f} mL",
			ModelsCount, maxHeight, totalVolumeMilliliters).c_str());
		this->Renderer->AddActor(AllModelsTextActor);
	}

	void setup() {
		if (!this->Renderer) {
			std::cerr << "Renderer is not set!" << std::endl;
			return;
		}
		// Center XYZ axes
		vtkNew<vtkAxesActor> XYZaxes;
		XYZaxes->SetTotalLength(10, 10, 10);
		XYZaxes->SetXAxisLabelText("");
		XYZaxes->SetYAxisLabelText("");
		XYZaxes->SetZAxisLabelText("");
		this->Renderer->AddActor(XYZaxes);

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
		this->Renderer->AddActor(sphereActor);

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
		this->Renderer->AddActor(planeActor);
	}

	void SetRenderer(vtkRenderer* renderer) {
		this->Renderer = renderer;
	}

private:
	bool loadSTL = true, loadPTS = true, loadHolder = true;
	size_t loadNC = 1;
	std::vector<std::string> FileNames;
	size_t CurrentIndex = 0;
	double NCsphereR = 0.75;
	double threshold = 0.5;
	vtkRenderer* Renderer = nullptr;
	vtkNew<vtkPoints> PTSpoints;
	vtkNew<vtkPoints> NCpoints;
};
//vtkStandardNewMacro(C_InteractorStyle);



int main(int argc, char* argv[]) {

	// Get the current directory path
	std::string folderPath = fs::current_path().string();
	std::set<std::string> uniqueFiles; // To store unique file paths without extension

	for (const auto& entry : fs::directory_iterator(folderPath)) {
		if (entry.is_regular_file()) { // Check if it's a regular file
			std::string extension = entry.path().extension().string();
			std::transform(extension.begin(), extension.end(), extension.begin(),
				[](unsigned char c) { return std::tolower(c); });

			if (extension == ".stl" || extension == ".pts" || extension == ".nc") {
				fs::path entryPath = entry.path();
				entryPath.replace_extension();
				uniqueFiles.insert(entryPath.string());
			}
		}
	}

	std::vector<std::string> Files(uniqueFiles.begin(), uniqueFiles.end());


	if (Files.empty()) {
		std::cout << std::endl;
		std::cout << Red << "\n\n        No STL or PTS or NC files in input folder\n" << ColorEnd << std::endl;
		std::cout << "        Press enter to continue...";
		std::cin.get();
		return EXIT_SUCCESS;
	}

	if (DEBUG) std::cout << Yellow << "        Preparing Mesh Viewer." << ColorEnd << std::endl;

	std::cout << std::endl;
	std::cout << Green << "        Keyboard Shortcuts:" << ColorEnd << std::endl;
	std::cout << Yellow << "          DOWN/UP    key  ->>  Opens Next/Previous file" << std::endl;
	std::cout << "          ENTER      key  ->>  Opens All STL files together" << std::endl;
	std::cout << "          M          Key  ->>  Enable/Disable STL View" << std::endl;
	std::cout << "          T          Key  ->>  Enable/Disable PTS View" << std::endl;
	std::cout << "          N          Key  ->>  Enable/Disable NC View" << std::endl;
	std::cout << "          F          Key  ->>  Enable/Disable Fixture Holder View" << std::endl;
	std::cout << "          LEFT/RIGHT key  ->>  Decrease/Increase Tool Diameter" << std::endl;
	std::cout << "          S/W        Key  ->>  Solid/WireFrame STL View" << ColorEnd << std::endl;


	vtkNew<vtkRenderer> Renderer;
	Renderer->SetBackground(bkR, bkG, bkB);
	Renderer->SetUseShadows(true);
	Renderer->SetUseDepthPeeling(1);        // Enable depth peeling
	Renderer->SetMaximumNumberOfPeels(100);
	Renderer->SetOcclusionRatio(0.1);

	// Window and interactor setup
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(800, 700);
	renderWindow->SetWindowName("AB Viewer tool (STL, PTS, NC)");
	renderWindow->AddRenderer(Renderer);
	renderWindow->SetAlphaBitPlanes(1);  // Enable alpha channel if supported
	renderWindow->SetMultiSamples(0);    // Disable multi-sampling to see depth peeling

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Custom Interaction Style
	vtkNew<C_InteractorStyle> style;
	style->SetDefaultRenderer(Renderer);
	style->SetRenderer(Renderer);
	style->SetFileNames(Files);
	//style->LoadAllSTLfiles();
	style->LoadFixtureHolder();
	style->LoadSTLfile(0);
	style->LoadPTSModel(0);
	style->LoadNCfile(0, 0.75, 1);
	style->setup();
	renderWindowInteractor->SetInteractorStyle(style);

	// Camera setups
	vtkNew<vtkCameraOrientationWidget> camOrientManipulator;
	camOrientManipulator->SetParentRenderer(Renderer);
	camOrientManipulator->On();

	//Renderer->AutomaticLightCreationOff();

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