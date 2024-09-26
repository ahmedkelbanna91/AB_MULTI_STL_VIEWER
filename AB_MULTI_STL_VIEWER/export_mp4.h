#pragma once

double BACKGROUND4VIDEO[3] = { 0.129, 0.129, 0.141 };
//double BACKGROUND4VIDEO[3] = { 1.0, 1.0, 1.0 };
double STLCOLOR4VIDEO[3] = { 0.7, 0.5, 0.3 };	// 1.0, 0.8, 0.4
double PTSCOLOR4VIDEO[3] = { 0.9, 0.8, 0.9 };	// 0.95, 0.95, 0.95
double NCCOLOR4VIDEO[3] = { 0.196, 0.51, 0.965 };	// 0.7, 0.9, 0.9
double SVGCOLOR4VIDEO[3] = { 0.0, 1.0, 0.0 };

double TEXTCOLOR4VIDEO[3] = { 1.0, 1.0, 1.0 };
//double TEXTCOLOR4VIDEO[3] = { 0.0, 0.0, 0.0 };

double ToolRadius4Video = 0.75;
double ToolOffsetLength4Video = -92.0;

bool ShowDetails = true;
int F_width = 1024;
int F_height = 768;
//int F_width = 1680;
//int F_height = 1050;
double fps = 1;
double FarPointThreshold = 0.6;


C_Actors LoadLaserMarkingZip4Video(const string& FilePath, string& CaseID) {
	string svgFileName = path(FilePath + ".svg").filename().string();
	path LaserMarkingZipPath(CaseID + ".zip");

	if (!exists(LaserMarkingZipPath)) return { nullptr, nullptr };

	unzFile zipfile = unzOpen(LaserMarkingZipPath.string().c_str());
	if (zipfile == nullptr) {
		cerr << "Cannot open zip file: " << LaserMarkingZipPath << endl;
		return { nullptr, nullptr };
	}

	if (unzLocateFile(zipfile, svgFileName.c_str(), 0) != UNZ_OK) {
		cerr << "Error: SVG file not found in the zip: " << svgFileName << endl;
		unzClose(zipfile);
		return { nullptr, nullptr };
	}

	if (unzOpenCurrentFile(zipfile) != UNZ_OK) {
		cerr << "Error: Cannot open the SVG file inside the zip: " << svgFileName << endl;
		unzClose(zipfile);
		return { nullptr, nullptr };
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
		return { nullptr, nullptr };
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

		vtkNew<vtkAssembly> SVGassembly;
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
			SVGassembly->AddPart(SVGglyphActor);


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
			textActor->GetProperty()->SetColor(SVGCOLOR4VIDEO);

			vtkNew<vtkTransform> transform;
			transform->Identity();
			transform->Translate(x1, y1, z1 + 0.5);
			transform->RotateZ(atan2(direction[1], direction[0]) * 180.0 / vtkMath::Pi());
			transform->RotateY((acos(direction[2]) * 180.0 / vtkMath::Pi()) - 90.0);
			transform->RotateX(15);
			textActor->SetUserTransform(transform);
			SVGassembly->AddPart(textActor);

			// Model SVG Text
			vtkNew<vtkTextActor> ZIPTextActor;
			ZIPTextActor->GetTextProperty()->SetFontSize(23);
			ZIPTextActor->GetTextProperty()->SetColor(SVGCOLOR4VIDEO);
			ZIPTextActor->SetPosition(20, 10); // 20, 100
			ZIPTextActor->GetTextProperty()->SetFontFamilyToArial();
			ZIPTextActor->GetTextProperty()->SetFontSize(23);
			ZIPTextActor->SetInput(std::format("Laser Marking: {}", LaserMarkingZipPath.filename().string()).c_str());

			cout << CYAN << "          >> ZIP: " << END << LaserMarkingZipPath.filename().string() << endl;
			return { SVGassembly, ZIPTextActor };
		}
		else {
			return { nullptr, nullptr };
		}
	}
	else {
		cerr << "Error: Could not find the expected <line> tag with the required attributes." << endl;
	}
}

C_Actors LoadNC4Video(const string& FilePath, double sphereRadius, double toolOffset, vtkSTLReader* ModelReader) {
	if (!exists(FilePath)) return { nullptr, nullptr };
	ifstream ncFile(FilePath);
	string line;

	vtkNew<vtkPoints> ncPoints;

	while (getline(ncFile, line)) {
		smatch match;
		if (regex_search(line, match, NCRegex)) {
			double B = stod(match[1].str());
			double C = stod(match[2].str());
			double X = stod(match[3].str());
			double Y = stod(match[4].str());
			double Z = stod(match[5].str());

			double B_rad = vtkMath::RadiansFromDegrees(B);
			double C_rad = vtkMath::RadiansFromDegrees(C);

			Eigen::Matrix4d RB, RC, Txyz, TFinal;

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
			Eigen::Vector4d toolTipPos = TFinal * Eigen::Vector4d(0, 0, toolOffset, 1);

			ncPoints->InsertNextPoint(toolTipPos(0), toolTipPos(1), toolTipPos(2));
		}
	}

	if (ncPoints->GetNumberOfPoints() > 0) {
		vtkNew<vtkPoints> inLoopPoints;
		vtkNew<vtkPoints> outLoopPoints;

		for (vtkIdType i = 1; i < ncPoints->GetNumberOfPoints() - 1; ++i) {
			double point[3];
			double prevPoint[3];
			double nextPoint[3];
			ncPoints->GetPoint(i, point);
			ncPoints->GetPoint(i - 1, prevPoint);
			ncPoints->GetPoint(i + 1, nextPoint);

			// check if a point is far from its neighbors
			double distanceToPrev = sqrt(
				pow(point[0] - prevPoint[0], 2) +
				pow(point[1] - prevPoint[1], 2) +
				pow(point[2] - prevPoint[2], 2));

			double distanceToNext = sqrt(
				pow(point[0] - nextPoint[0], 2) +
				pow(point[1] - nextPoint[1], 2) +
				pow(point[2] - nextPoint[2], 2));

			if (!(distanceToPrev >= FarPointThreshold) || !(distanceToNext >= FarPointThreshold)) {
				inLoopPoints->InsertNextPoint(point);
			}
			else {
				outLoopPoints->InsertNextPoint(point);
			}
		}

		vtkNew<vtkAssembly> NCassembly;

		string STLPath = path(FilePath).replace_extension(".stl").string();
		if (exists(STLPath) && ModelReader && ModelReader->GetOutput() != nullptr) {

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


			if (nonSurfacePoints->GetNumberOfPoints() > 0) {
				// Visualize in-loop non-surface points (red color)
				vtkNew<vtkPolyData> nonSurfacePolyData;
				nonSurfacePolyData->SetPoints(nonSurfacePoints);
				vtkNew<vtkVertexGlyphFilter> nonSurfaceVertexFilter;
				nonSurfaceVertexFilter->SetInputData(nonSurfacePolyData);
				nonSurfaceVertexFilter->Update();
				vtkNew<vtkSphereSource> nonSurfaceSphereSource;
				nonSurfaceSphereSource->SetRadius(sphereRadius);
				vtkNew<vtkGlyph3D> nonSurfaceGlyph3D;
				nonSurfaceGlyph3D->SetSourceConnection(nonSurfaceSphereSource->GetOutputPort());
				nonSurfaceGlyph3D->SetInputConnection(nonSurfaceVertexFilter->GetOutputPort());
				nonSurfaceGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> nonSurfaceGlyphMapper;
				nonSurfaceGlyphMapper->SetInputConnection(nonSurfaceGlyph3D->GetOutputPort());
				vtkNew<vtkActor> nonSurfaceGlyphActor;
				nonSurfaceGlyphActor->SetMapper(nonSurfaceGlyphMapper);
				nonSurfaceGlyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  // REDcolor
				NCassembly->AddPart(nonSurfaceGlyphActor);
			}

			if (surfacePoints->GetNumberOfPoints() > 0) {
				// Visualize in-loop surface points (original color)
				vtkNew<vtkPolyData> surfacePolyData;
				surfacePolyData->SetPoints(surfacePoints);
				vtkNew<vtkVertexGlyphFilter> surfaceVertexFilter;
				surfaceVertexFilter->SetInputData(surfacePolyData);
				surfaceVertexFilter->Update();
				vtkNew<vtkSphereSource> surfaceSphereSource;
				surfaceSphereSource->SetRadius(sphereRadius);
				vtkNew<vtkGlyph3D> surfaceGlyph3D;
				surfaceGlyph3D->SetSourceConnection(surfaceSphereSource->GetOutputPort());
				surfaceGlyph3D->SetInputConnection(surfaceVertexFilter->GetOutputPort());
				surfaceGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> surfaceGlyphMapper;
				surfaceGlyphMapper->SetInputConnection(surfaceGlyph3D->GetOutputPort());
				vtkNew<vtkActor> surfaceGlyphActor;
				surfaceGlyphActor->SetMapper(surfaceGlyphMapper);
				surfaceGlyphActor->GetProperty()->SetColor(NCCOLOR4VIDEO);  // Original color
				NCassembly->AddPart(surfaceGlyphActor);
			}
		}
		else {
			if (inLoopPoints->GetNumberOfPoints() > 0) {
				// Visualize in-loop points
				vtkNew<vtkPolyData> inLoopPolyData;
				inLoopPolyData->SetPoints(inLoopPoints);
				vtkNew<vtkVertexGlyphFilter> inLoopVertexFilter;
				inLoopVertexFilter->SetInputData(inLoopPolyData);
				inLoopVertexFilter->Update();
				vtkNew<vtkSphereSource> inLoopSphereSource;
				inLoopSphereSource->SetRadius(sphereRadius);
				vtkNew<vtkGlyph3D> inLoopGlyph3D;
				inLoopGlyph3D->SetSourceConnection(inLoopSphereSource->GetOutputPort());
				inLoopGlyph3D->SetInputConnection(inLoopVertexFilter->GetOutputPort());
				inLoopGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> inLoopGlyphMapper;
				inLoopGlyphMapper->SetInputConnection(inLoopGlyph3D->GetOutputPort());
				vtkNew<vtkActor> inLoopGlyphActor;
				inLoopGlyphActor->SetMapper(inLoopGlyphMapper);
				inLoopGlyphActor->GetProperty()->SetColor(NCCOLOR4VIDEO);
				NCassembly->AddPart(inLoopGlyphActor);
			}
		}

		if (outLoopPoints->GetNumberOfPoints() > 0) {
			// Visualize out-of-loop points
			vtkNew<vtkPolyData> outLoopPolyData;
			outLoopPolyData->SetPoints(outLoopPoints);
			vtkNew<vtkVertexGlyphFilter> outLoopVertexFilter;
			outLoopVertexFilter->SetInputData(outLoopPolyData);
			outLoopVertexFilter->Update();
			vtkNew<vtkSphereSource> outLoopSphereSource;
			outLoopSphereSource->SetRadius(sphereRadius);
			vtkNew<vtkGlyph3D> outLoopGlyph3D;
			outLoopGlyph3D->SetSourceConnection(outLoopSphereSource->GetOutputPort());
			outLoopGlyph3D->SetInputConnection(outLoopVertexFilter->GetOutputPort());
			outLoopGlyph3D->ScalingOff();
			vtkNew<vtkPolyDataMapper> outLoopGlyphMapper;
			outLoopGlyphMapper->SetInputConnection(outLoopGlyph3D->GetOutputPort());
			vtkNew<vtkActor> outLoopGlyphActor;
			outLoopGlyphActor->SetMapper(outLoopGlyphMapper);
			outLoopGlyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
			NCassembly->AddPart(outLoopGlyphActor);
		}

		// Model NC Text
		vtkNew<vtkTextActor> NCTextActor;
		NCTextActor->GetTextProperty()->SetFontSize(23);
		NCTextActor->GetTextProperty()->SetColor(NCCOLOR4VIDEO);
		NCTextActor->SetPosition(20, 1); // 20, 100
		NCTextActor->GetTextProperty()->SetFontFamilyToArial();
		NCTextActor->GetTextProperty()->SetFontSize(23);
		NCTextActor->SetInput(std::format("Tool Diameter: {:.2f} mm\nNC: {}",
			sphereRadius * 2.0, path(FilePath).filename().string()).c_str());

		cout << CYAN << "          >> NC: " << END << path(FilePath).filename().string() << endl;
		return { NCassembly, NCTextActor };
	}
	else {
		return { nullptr, nullptr };
	}
}

C_Actors LoadPTS4Video(const string& FilePath, double sphereRadius, vtkSTLReader* ModelReader) {
	if (!exists(FilePath)) return { nullptr, nullptr };
	vtkNew<vtkPoints> ptsPoints;

	ifstream ptsFile(FilePath);
	string line;

	while (getline(ptsFile, line)) {
		smatch match;
		if (regex_match(line, match, PTSRegex) && match.size() >= 4) {
			double x = stod(match[1].str());
			double y = stod(match[2].str());
			double z = stod(match[3].str());
			ptsPoints->InsertNextPoint(x, y, z);
		}
	}

	if (ptsPoints->GetNumberOfPoints() > 0) {
		vtkNew<vtkPoints> inLoopPoints;
		vtkNew<vtkPoints> outLoopPoints;

		for (vtkIdType i = 1; i < ptsPoints->GetNumberOfPoints() - 1; ++i) {
			double point[3];
			double prevPoint[3];
			double nextPoint[3];
			ptsPoints->GetPoint(i, point);
			ptsPoints->GetPoint(i - 1, prevPoint);
			ptsPoints->GetPoint(i + 1, nextPoint);

			// check if a point is far from its neighbors
			double distanceToPrev = sqrt(
				pow(point[0] - prevPoint[0], 2) +
				pow(point[1] - prevPoint[1], 2) +
				pow(point[2] - prevPoint[2], 2));

			double distanceToNext = sqrt(
				pow(point[0] - nextPoint[0], 2) +
				pow(point[1] - nextPoint[1], 2) +
				pow(point[2] - nextPoint[2], 2));

			if (!(distanceToPrev >= FarPointThreshold) || !(distanceToNext >= FarPointThreshold)) {
				inLoopPoints->InsertNextPoint(point);
			}
			else {
				outLoopPoints->InsertNextPoint(point);
			}
		}

		vtkNew<vtkAssembly> PTSassembly;

		string STLPath = path(FilePath).replace_extension(".stl").string();

		if (exists(STLPath) && ModelReader && ModelReader->GetOutput() != nullptr) {

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


			if (nonSurfacePoints->GetNumberOfPoints() > 0) {
				// Visualize in-loop non-surface points (red color)
				vtkNew<vtkPolyData> nonSurfacePolyData;
				nonSurfacePolyData->SetPoints(nonSurfacePoints);
				vtkNew<vtkVertexGlyphFilter> nonSurfaceVertexFilter;
				nonSurfaceVertexFilter->SetInputData(nonSurfacePolyData);
				nonSurfaceVertexFilter->Update();
				vtkNew<vtkSphereSource> nonSurfaceSphereSource;
				nonSurfaceSphereSource->SetRadius(sphereRadius);
				vtkNew<vtkGlyph3D> nonSurfaceGlyph3D;
				nonSurfaceGlyph3D->SetSourceConnection(nonSurfaceSphereSource->GetOutputPort());
				nonSurfaceGlyph3D->SetInputConnection(nonSurfaceVertexFilter->GetOutputPort());
				nonSurfaceGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> nonSurfaceGlyphMapper;
				nonSurfaceGlyphMapper->SetInputConnection(nonSurfaceGlyph3D->GetOutputPort());
				vtkNew<vtkActor> nonSurfaceGlyphActor;
				nonSurfaceGlyphActor->SetMapper(nonSurfaceGlyphMapper);
				nonSurfaceGlyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  // REDcolor
				PTSassembly->AddPart(nonSurfaceGlyphActor);
			}

			if (surfacePoints->GetNumberOfPoints() > 0) {
				// Visualize in-loop surface points (original color)
				vtkNew<vtkPolyData> surfacePolyData;
				surfacePolyData->SetPoints(surfacePoints);
				vtkNew<vtkVertexGlyphFilter> surfaceVertexFilter;
				surfaceVertexFilter->SetInputData(surfacePolyData);
				surfaceVertexFilter->Update();
				vtkNew<vtkSphereSource> surfaceSphereSource;
				surfaceSphereSource->SetRadius(sphereRadius);
				vtkNew<vtkGlyph3D> surfaceGlyph3D;
				surfaceGlyph3D->SetSourceConnection(surfaceSphereSource->GetOutputPort());
				surfaceGlyph3D->SetInputConnection(surfaceVertexFilter->GetOutputPort());
				surfaceGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> surfaceGlyphMapper;
				surfaceGlyphMapper->SetInputConnection(surfaceGlyph3D->GetOutputPort());
				vtkNew<vtkActor> surfaceGlyphActor;
				surfaceGlyphActor->SetMapper(surfaceGlyphMapper);
				surfaceGlyphActor->GetProperty()->SetColor(PTSCOLOR4VIDEO);  // Original color
				PTSassembly->AddPart(surfaceGlyphActor);
			}
		}
		else {
			if (inLoopPoints->GetNumberOfPoints() > 0) {
				// Visualize in-loop points
				vtkNew<vtkPolyData> inLoopPolyData;
				inLoopPolyData->SetPoints(inLoopPoints);
				vtkNew<vtkVertexGlyphFilter> inLoopVertexFilter;
				inLoopVertexFilter->SetInputData(inLoopPolyData);
				inLoopVertexFilter->Update();
				vtkNew<vtkSphereSource> inLoopSphereSource;
				inLoopSphereSource->SetRadius(sphereRadius);
				vtkNew<vtkGlyph3D> inLoopGlyph3D;
				inLoopGlyph3D->SetSourceConnection(inLoopSphereSource->GetOutputPort());
				inLoopGlyph3D->SetInputConnection(inLoopVertexFilter->GetOutputPort());
				inLoopGlyph3D->ScalingOff();
				vtkNew<vtkPolyDataMapper> inLoopGlyphMapper;
				inLoopGlyphMapper->SetInputConnection(inLoopGlyph3D->GetOutputPort());
				vtkNew<vtkActor> inLoopGlyphActor;
				inLoopGlyphActor->SetMapper(inLoopGlyphMapper);
				inLoopGlyphActor->GetProperty()->SetColor(PTSCOLOR4VIDEO);
				PTSassembly->AddPart(inLoopGlyphActor);
			}
		}

		if (outLoopPoints->GetNumberOfPoints() > 0) {
			// Visualize out-of-loop points
			vtkNew<vtkPolyData> outLoopPolyData;
			outLoopPolyData->SetPoints(outLoopPoints);
			vtkNew<vtkVertexGlyphFilter> outLoopVertexFilter;
			outLoopVertexFilter->SetInputData(outLoopPolyData);
			outLoopVertexFilter->Update();
			vtkNew<vtkSphereSource> outLoopSphereSource;
			outLoopSphereSource->SetRadius(sphereRadius);
			vtkNew<vtkGlyph3D> outLoopGlyph3D;
			outLoopGlyph3D->SetSourceConnection(outLoopSphereSource->GetOutputPort());
			outLoopGlyph3D->SetInputConnection(outLoopVertexFilter->GetOutputPort());
			outLoopGlyph3D->ScalingOff();
			vtkNew<vtkPolyDataMapper> outLoopGlyphMapper;
			outLoopGlyphMapper->SetInputConnection(outLoopGlyph3D->GetOutputPort());
			vtkNew<vtkActor> outLoopGlyphActor;
			outLoopGlyphActor->SetMapper(outLoopGlyphMapper);
			outLoopGlyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
			PTSassembly->AddPart(outLoopGlyphActor);
		}

		// Model PTS Text
		vtkNew<vtkTextActor> PTSTextActor;
		PTSTextActor->GetTextProperty()->SetFontSize(23);
		PTSTextActor->GetTextProperty()->SetColor(PTSCOLOR4VIDEO);
		PTSTextActor->SetPosition(20, 1); // 20, 77
		PTSTextActor->GetTextProperty()->SetFontFamilyToArial();
		PTSTextActor->GetTextProperty()->SetFontSize(23);
		PTSTextActor->SetInput(std::format("PTS: {} -> {} Points", path(FilePath).filename().string(), ptsPoints->GetNumberOfPoints()).c_str());

		cout << CYAN << "          >> PTS: " << END << path(FilePath).filename().string() << endl;
		return { PTSassembly, PTSTextActor };
	}
	else {
		return { nullptr, nullptr };
	}
}

C_Actors LoadSTL4Video(const string& FilePath, vtkSmartPointer<vtkSTLReader>& stlReader) {
	if (!exists(FilePath)) {
		stlReader = nullptr;
		return { nullptr, nullptr };
	}

	stlReader->SetFileName(FilePath.c_str());
	stlReader->Update();
	vtkNew<vtkPolyDataMapper> stlMapper;
	stlMapper->SetInputConnection(stlReader->GetOutputPort());
	vtkNew<vtkActor> stlActor;
	stlActor->SetMapper(stlMapper);
	stlActor->GetProperty()->SetColor(STLCOLOR4VIDEO);
	double bounds[6];
	stlActor->GetBounds(bounds);
	double Height_ = bounds[5] - bounds[4];
	vtkNew<vtkMassProperties> massProperties;
	massProperties->SetInputConnection(stlReader->GetOutputPort());
	massProperties->Update();
	double volumeMilliliters = massProperties->GetVolume() / 1000.0;

	// Model Text
	vtkNew<vtkTextActor> STLTextActor;
	STLTextActor->GetTextProperty()->SetFontSize(23);
	STLTextActor->GetTextProperty()->SetColor(STLCOLOR4VIDEO);
	STLTextActor->SetPosition(20, 10);
	STLTextActor->GetTextProperty()->SetFontFamilyToArial();
	STLTextActor->GetTextProperty()->SetFontSize(23);
	STLTextActor->SetInput(std::format("ID: {}\nHeight: {:.2f} mm\nVolume: {:.2f} mL",
		path(FilePath).filename().string(), Height_, volumeMilliliters).c_str());

	cout << CYAN << "         >> STL: " << END << path(FilePath).filename().string() << endl;
	return { stlActor, STLTextActor };
}

void AddTextLabel4Video(vtkRenderer* renderer, const string& text, int X, int Y, int fontsize, bool B, bool I) {
	vtkNew<vtkTextActor> textActor;
	textActor->SetInput(text.c_str());
	textActor->SetPosition(X, Y);
	textActor->GetTextProperty()->SetFontFamilyToArial();
	textActor->GetTextProperty()->SetFontSize(fontsize);
	textActor->GetTextProperty()->SetColor(TEXTCOLOR4VIDEO);
	textActor->GetTextProperty()->SetBold(B);
	textActor->GetTextProperty()->SetItalic(I);
	renderer->AddActor2D(textActor);
}

void AddActors2Renderer4Video(vtkRenderer* renderer, vector<string>& FilePaths, string& FilePath,
	string& CaseID, int& index) {

	renderer->RemoveAllViewProps();
	renderer->SetBackground(BACKGROUND4VIDEO);

	C_Actors stlActors = { nullptr, nullptr };
	C_Actors ptsActors = { nullptr, nullptr };
	C_Actors ncActors = { nullptr, nullptr };
	C_Actors svgActors = { nullptr, nullptr };

	string stlPath = FilePath + ".stl";
	string ptsPath = FilePath + ".pts";
	string ncPath = FilePath + ".NC";
	string svgPath = FilePath;

	vtkSmartPointer<vtkSTLReader> ModelReader = vtkSmartPointer<vtkSTLReader>::New();

	vector<string> views = { "Left", "Right", "Front", "Back" };
	size_t numViews = views.size();
	size_t cols = 2;
	size_t rows = (numViews + cols - 1) / cols;
	for (int i = 0; i < numViews; ++i) {
		double x0 = (i % cols) * 0.5;
		double y0 = (rows - 1 - i / cols) * 0.5;
		double x1 = x0 + 0.5;
		double y1 = y0 + 0.5;

		vtkNew<vtkRenderer> subRenderer;
		subRenderer->SetViewport(x0, y0, x1, y1);
		subRenderer->SetBackground(BACKGROUND4VIDEO);
		subRenderer->SetUseShadows(true);
		subRenderer->SetUseDepthPeeling(1);        // Enable depth peeling
		subRenderer->SetMaximumNumberOfPeels(100);
		subRenderer->SetOcclusionRatio(0.7);

		stlActors = actorCache.contains(stlPath) ? actorCache[stlPath] : LoadSTL4Video(stlPath, ModelReader);
		if (!actorCache.contains(stlPath)) actorCache[stlPath] = stlActors;

		ptsActors = actorCache.contains(ptsPath) ? actorCache[ptsPath] : LoadPTS4Video(ptsPath, ToolRadius4Video - 0.1, ModelReader);
		if (!actorCache.contains(ptsPath)) actorCache[ptsPath] = ptsActors;

		ncActors = actorCache.contains(ncPath) ? actorCache[ncPath] : LoadNC4Video(ncPath, ToolRadius4Video, ToolOffsetLength4Video, ModelReader);
		if (!actorCache.contains(ncPath)) actorCache[ncPath] = ncActors;

		svgActors = actorCache.contains(svgPath) ? actorCache[svgPath] : LoadLaserMarkingZip4Video(svgPath, CaseID);
		if (!actorCache.contains(svgPath)) actorCache[svgPath] = svgActors;

		if (ptsActors.Actor) subRenderer->AddActor(ptsActors.Actor);
		if (i == 0 && ptsActors.textActor && !ncActors.textActor) subRenderer->AddActor(ptsActors.textActor);
		else if (i == 1 && ncActors.textActor && ptsActors.textActor) subRenderer->AddActor(ptsActors.textActor);

		if (ncActors.Actor) subRenderer->AddActor(ncActors.Actor);
		if (i == 0 && ncActors.textActor && !ptsActors.textActor) subRenderer->AddActor(ncActors.textActor);
		else if (i == 0 && ncActors.textActor && ptsActors.textActor) subRenderer->AddActor(ncActors.textActor);

		if (stlActors.Actor) subRenderer->AddActor(stlActors.Actor);
		if (i == 2 && stlActors.textActor) subRenderer->AddActor(stlActors.textActor);

		if (svgActors.Actor) subRenderer->AddActor(svgActors.Actor);
		if (i == 3 && svgActors.textActor) subRenderer->AddActor(svgActors.textActor);


		// Get the size of the viewport
		int* windowSize;
		windowSize = renderer->GetRenderWindow()->GetSize();
		int VPWidth = static_cast<int>((x1 - x0) * windowSize[0]);
		int VPHeight = static_cast<int>((y1 - y0) * windowSize[1]);

		string DateTime;
		GetCurrentDateTime(DateTime);
		AddTextLabel4Video(subRenderer, views[i] + " View", VPWidth - 150, VPHeight - 50, 20, true, true);
		if (i == 0) {
			AddTextLabel4Video(subRenderer, CaseID, 20, VPHeight - 50, 26, true, false);
			AddTextLabel4Video(subRenderer, std::format("# {}/{}", index, FilePaths.size()).c_str(), 20, VPHeight - 80, 26, false, true);
		}
		if (i == 1) {
			AddTextLabel4Video(subRenderer, Header + DateTime, VPWidth - 250, VPHeight - 15, 8, false, false);
		}

		vtkNew<vtkCamera> camera;
		camera->SetViewUp(0, 0, 1);

		if (views[i] == "Left") {
			camera->SetPosition(1.0, 0.0, 0.5);
			camera->Azimuth(35);
		}
		else if (views[i] == "Right") {
			camera->SetPosition(-1.0, 0.0, 0.5);
			camera->Azimuth(-35);
		}
		else if (views[i] == "Front") {
			camera->SetPosition(0.0, 1.0, 0.2);
		}
		else if (views[i] == "Back") {
			camera->SetPosition(0.0, -1.0, 1.8);
		}
		subRenderer->SetActiveCamera(camera);
		subRenderer->ResetCamera();
		subRenderer->ResetCameraClippingRange();
		subRenderer->GetActiveCamera()->Zoom(1.3);

		renderer->GetRenderWindow()->AddRenderer(subRenderer);
	}
}

vector<unsigned char> GetImageData4Video(vector<string>& FilePaths, string FilePath, string& CaseID,
	int& index, int& width, int& height) {

	vtkNew<vtkRenderer> mainRenderer;
	mainRenderer->SetUseShadows(true);
	mainRenderer->SetUseDepthPeeling(1);
	mainRenderer->SetMaximumNumberOfPeels(100);
	mainRenderer->SetOcclusionRatio(0.7);
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(width, height);
	renderWindow->SetAlphaBitPlanes(1);
	//renderWindow->SetMultiSamples(4);
	renderWindow->AddRenderer(mainRenderer);

	AddActors2Renderer4Video(mainRenderer, FilePaths, FilePath, CaseID, index); // enabled quad view

	vtkNew<vtkCamera> camera;
	camera->SetPosition(0.0, 1.0, 0.2);
	camera->SetViewUp(0, 0, 1);

	mainRenderer->SetActiveCamera(camera);
	mainRenderer->ResetCamera();
	mainRenderer->ResetCameraClippingRange();
	mainRenderer->GetActiveCamera()->Zoom(1.5);
	//mainRenderer->GetActiveCamera()->Azimuth(rot);
	renderWindow->OffScreenRenderingOn();
	renderWindow->Render();

	vtkNew<vtkWindowToImageFilter> windowToImageFilter;
	windowToImageFilter->SetInput(renderWindow);
	windowToImageFilter->SetInputBufferTypeToRGB();
	windowToImageFilter->Update();

	//string PngPath = FilePath + ".png";
	//vtkNew<vtkPNGWriter> PNGWriter;
	//PNGWriter->SetFileName(PngPath.c_str());
	//PNGWriter->SetInputConnection(windowToImageFilter->GetOutputPort());
	//PNGWriter->Write();
	//if (!PNGWriter->GetErrorCode()) {
	//	cout << "PNG file saved successfully as: " << PngPath << endl;
	//}
	//else {
	//	cerr << "Error: Failed to save PNG file!" << endl;
	//}

	vtkNew<vtkJPEGWriter> JPEGWriter;
	JPEGWriter->WriteToMemoryOn();
	JPEGWriter->SetInputConnection(windowToImageFilter->GetOutputPort());
	//imageWriter->SetQuality(100);
	JPEGWriter->Write();

	vtkUnsignedCharArray* jpegData = JPEGWriter->GetResult();
	if (!jpegData) {
		cerr << RED << "Error: No JPEG data returned." << END << endl;
		return {};
	}

	vector<unsigned char> buffer(jpegData->GetNumberOfTuples());
	memcpy(buffer.data(), jpegData->GetPointer(0), buffer.size());

	return buffer;
}


bool SaveMP4(string& mp4Filename, vector<vector<unsigned char>>& allImageData) {
	cout << YELLOW << "     << Opening MP4 File: " << END << mp4Filename << endl;

	VideoWriter mp4Writer(mp4Filename, VideoWriter::fourcc('X', '2', '6', '4'), fps, Size(F_width, F_height));

	if (!mp4Writer.isOpened()) {
		cerr << RED << "Could not open the video file for writing!" << END << endl;
		return false;
	}

	int ImageCount = 0;
	for (auto& imageData : allImageData) {
		Mat jpgData(1, static_cast<int>(imageData.size()), CV_8UC1, (void*)imageData.data());
		Mat frame = imdecode(jpgData, IMREAD_COLOR);

		if (frame.empty()) {
			cerr << RED << "Error: Failed to decode JPEG data." << END << endl;
			continue;
		}

		if (frame.cols != F_width || frame.rows != F_height) {
			resize(frame, frame, Size(F_width, F_height));
		}

		mp4Writer.write(frame);

		ImageCount++;
		if (ImageCount == 1 || ImageCount == allImageData.size()) {
			for (int i = 0; i < fps; i++) {
				mp4Writer.write(frame);
			}
		}
	}

	mp4Writer.release();
	cout << GREEN << "     << MP4 File Saved." << END << endl;
	return true;
}

bool SaveGIF(string& gifFilename, vector<vector<unsigned char>>& allImageData) {
	cout << YELLOW << "     << Opening GIF File: " << END << gifFilename << endl;

	FreeImage_Initialise();

	FIMULTIBITMAP* multiBitmap = FreeImage_OpenMultiBitmap(FIF_GIF, gifFilename.c_str(), TRUE, FALSE);

	if (!multiBitmap) {
		cerr << "Error: Failed to create multi-page GIF!" << endl;
		FreeImage_DeInitialise();
		return false;
	}

	for (auto& imageData : allImageData) {
		FIMEMORY* memoryStream = FreeImage_OpenMemory(imageData.data(), static_cast<DWORD>(imageData.size()));

		FIBITMAP* bitmap = FreeImage_LoadFromMemory(FIF_JPEG, memoryStream, JPEG_DEFAULT);
		FreeImage_CloseMemory(memoryStream);

		if (!bitmap) {
			cerr << RED << "Error: Failed to load JPEG image from memory!" << END << endl;
			continue;
		}

		//FIBITMAP* indexedBitmap = FreeImage_ColorQuantizeEx(bitmap, FIQ_WUQUANT, 32); // Convert to indexed color
		//FreeImage_Unload(bitmap);

		FIBITMAP* resizedBitmap = FreeImage_Rescale(bitmap, F_width, F_height, FILTER_BICUBIC);
		FreeImage_Unload(bitmap);
		FIBITMAP* indexedBitmap = FreeImage_ColorQuantizeEx(resizedBitmap, FIQ_WUQUANT, 32);
		FreeImage_Unload(resizedBitmap);

		FITAG* tag = FreeImage_CreateTag();
		if (tag) {
			DWORD dwFrameTime = (DWORD)((1000.0f / fps) + 0.5f);  // Frame delay in milliseconds
			FreeImage_SetTagKey(tag, "FrameTime");
			FreeImage_SetTagType(tag, FIDT_LONG);   // Data type is a long integer
			FreeImage_SetTagCount(tag, 1);          // One value
			FreeImage_SetTagLength(tag, 4);         // Length of the value in bytes (DWORD is 4 bytes)
			FreeImage_SetTagValue(tag, &dwFrameTime); // Set the actual value (frame time)
			FreeImage_SetMetadata(FIMD_ANIMATION, indexedBitmap, FreeImage_GetTagKey(tag), tag);
			FreeImage_DeleteTag(tag);  // Clean up the tag after setting the metadata
		}

		FreeImage_AppendPage(multiBitmap, indexedBitmap);
		FreeImage_Unload(indexedBitmap);
	}

	if (!FreeImage_CloseMultiBitmap(multiBitmap)) {
		cerr << "Error: Failed to save the GIF file!" << endl;
		return false;
	}
	else {
		cout << GREEN << "     << GIF File Saved." << END << endl;
	}

	FreeImage_DeInitialise();
}


bool SaveVideo() {
	string Main_path = current_path().string() + "\\";
	bool has_stl_files = false;
	vector<string> filenames;

	for (auto& entry : directory_iterator(Main_path)) {
		if (entry.path().extension() == ".stl" || entry.path().extension() == ".STL") {
			string filename_WE = entry.path().stem().string();
			filenames.push_back(filename_WE);
			has_stl_files = true;
		}
	}

	if (!has_stl_files) {
		cout << RED << "        No STL files in input folder\n" << END << endl;
		cout << "        Press enter to continue...";
		cin.get();
		return EXIT_SUCCESS;
	}

	map<string, vector<string>> Cases;
	for (auto& filename : filenames) {
		smatch matches;
		if (regex_match(filename, matches, filenameRegex) && matches.size() > 1) {
			string CaseID = matches[1].str();
			Cases[CaseID].push_back(filename);
		}
	}

	vector<string> SortedCaseIDs;
	for (auto& pair : Cases)
		SortedCaseIDs.push_back(pair.first);
	sort(SortedCaseIDs.begin(), SortedCaseIDs.end(), [](const string& a, const string& b) {
		return stoi(a) < stoi(b);
		});

	int Count = 0;
	bool success = true;

	for (auto& CaseID : SortedCaseIDs) {
		Count++;
		cout << (Count < 10 ? "   0" : "   ") + to_string(Count) << " > " << GREEN << "Processing CaseID: " << END << CaseID << endl;

		vector<vector<unsigned char>> allImageData;
		int countindex = 0;
		for (string fileID : Cases[CaseID]) {
			countindex++;
			vector<unsigned char> imageData = GetImageData4Video(Cases[CaseID], Main_path + fileID, CaseID, countindex, F_width, F_height); // Quad View
			if (imageData.empty()) {
				cerr << RED << "  Error: JPEG buffer is empty for " << fileID << END << endl;
				success = false;
				continue;
			}
			allImageData.push_back(imageData);
		}

		if (!allImageData.empty()) {
			string mp4Filename = CaseID + "_TL_Views.mp4";
			string gifFilename = CaseID + "_TL_Views.gif";

			bool mp4Success = SaveMP4(mp4Filename, allImageData);
			bool gifSuccess = SaveGIF(gifFilename, allImageData);

			if (!mp4Success || !gifSuccess) {
				success = false;
			}
		}
	}

	return success && Count > 0;
}


