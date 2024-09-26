#pragma once

double BACKGROUND[3] = { 1.0, 1.0, 1.0 };	// 1.0, 1.0, 1.0 

double STLTEXTCOLOR[3] = { 0.0, 0.0, 0.1 };	// 0.0, 0.0, 0.0
double TEXTCOLOR[3] = { 0.0, 0.0, 0.0 };

double STLCOLOR[3] = { 1.0, 1.0, 1.0 };	// 1.0, 0.8, 0.4
double PTSCOLOR[3] = { 0.6, 0.6, 0.6 };	// 0.95, 0.95, 0.95
double NCCOLOR[3] = { 0.5, 0.5, 0.5 };	// 0.7, 0.9, 0.9
double SVGCOLOR[3] = { 0.0, 0.0, 0.0 };

double ToolDiameter = 0.75;
double ToolOffsetLength = -92.0;

string Header = "AB3DVIEW2MP4/GIF/PDF V2 (STL,PTS,NC) - ";

HPDF_Date GetCurrentDateTime(string& DateTime) {
	time_t t = time(nullptr);
	tm tm;
	localtime_s(&tm, &t);

	ostringstream oss;
	oss << put_time(&tm, "%d/%m/%Y, %H:%M:%S");
	DateTime = oss.str();

	long timezoneOffset;
	_get_timezone(&timezoneOffset);
	if (tm.tm_isdst > 0)  timezoneOffset -= 3600;

	int offsetHours = static_cast<int>(timezoneOffset / 3600);
	int offsetMinutes = static_cast<int>((timezoneOffset % 3600) / 60);

	HPDF_Date hpdfDate;
	hpdfDate.year = tm.tm_year + 1900;
	hpdfDate.month = tm.tm_mon + 1;
	hpdfDate.day = tm.tm_mday;
	hpdfDate.hour = tm.tm_hour;
	hpdfDate.minutes = tm.tm_min;
	hpdfDate.seconds = tm.tm_sec;
	hpdfDate.ind = (offsetHours >= 0) ? '-' : '+';
	hpdfDate.off_hour = abs(offsetHours);
	hpdfDate.off_minutes = abs(offsetMinutes);

	return hpdfDate;
}

void ExportToPDF(const string& pdfFile, const vector<unsigned char>& UPPER_ImageData, const vector<unsigned char>& LOWER_ImageData) {
	HPDF_Doc pdf = HPDF_New(NULL, NULL);
	HPDF_SetCompressionMode(pdf, HPDF_COMP_ALL);

	string DateTime;
	HPDF_Date creationDate = GetCurrentDateTime(DateTime);
	HPDF_SetInfoAttr(pdf, HPDF_INFO_TITLE, "AB3DVIEW2MP4/GIF/PDF V2");
	HPDF_SetInfoAttr(pdf, HPDF_INFO_AUTHOR, "Ahmed Elbanna");
	HPDF_SetInfoAttr(pdf, HPDF_INFO_SUBJECT, "STL & PTS/NC quadViews");
	HPDF_SetInfoAttr(pdf, HPDF_INFO_KEYWORDS, "STL, PTS, NC, 3D, PDF");
	HPDF_SetInfoDateAttr(pdf, HPDF_INFO_CREATION_DATE, creationDate);
	//HPDF_SetInfoDateAttr(pdf, HPDF_INFO_MOD_DATE, creationDate);

	HPDF_Page page = HPDF_AddPage(pdf);
	HPDF_Page_SetSize(page, HPDF_PAGE_SIZE_A4, HPDF_PAGE_PORTRAIT);

	float margin = 15.0f;
	float pageWidth = HPDF_Page_GetWidth(page);
	float pageHeight = HPDF_Page_GetHeight(page);
	float imageWidth = pageWidth - 2 * margin;
	float imageHeight = (pageHeight - 2 * margin) / 2;

	if (!UPPER_ImageData.empty()) {
		HPDF_Image UPPER_Image = HPDF_LoadJpegImageFromMem(pdf, UPPER_ImageData.data(), UPPER_ImageData.size());
		HPDF_Page_DrawImage(page, UPPER_Image, margin, pageHeight - margin - imageHeight, imageWidth, imageHeight);
	}

	if (!LOWER_ImageData.empty()) {
		HPDF_Image LOWER_Image = HPDF_LoadJpegImageFromMem(pdf, LOWER_ImageData.data(), LOWER_ImageData.size());
		if (!UPPER_ImageData.empty()) {
			HPDF_Page_DrawImage(page, LOWER_Image, margin, margin, imageWidth, imageHeight);
			HPDF_Page_MoveTo(page, margin, pageHeight - margin - imageHeight);
			HPDF_Page_LineTo(page, pageWidth - margin, pageHeight - margin - imageHeight);
			HPDF_Page_Stroke(page);
		}
		else {
			HPDF_Page_DrawImage(page, LOWER_Image, margin, pageHeight - margin - imageHeight, imageWidth, imageHeight);
		}
	}

	string headerText = Header + DateTime;
	HPDF_Page_SetFontAndSize(page, HPDF_GetFont(pdf, "Helvetica", NULL), 4);
	float textWidth = HPDF_Page_TextWidth(page, headerText.c_str());
	HPDF_Page_BeginText(page);
	HPDF_Page_TextOut(page, (pageWidth - textWidth) / 2, (pageHeight / 2) + (margin / 2), headerText.c_str());
	HPDF_Page_TextOut(page, (pageWidth - textWidth) / 2, margin, headerText.c_str());
	HPDF_Page_EndText(page);

	HPDF_SaveToFile(pdf, pdfFile.c_str());
	HPDF_Free(pdf);
}


C_Actors LoadLaserMarkingZip(const string& FilePath, string& CaseID) {
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
			textActor->GetProperty()->SetColor(SVGCOLOR);

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
			ZIPTextActor->GetTextProperty()->SetFontSize(30);
			ZIPTextActor->GetTextProperty()->SetColor(SVGCOLOR);
			ZIPTextActor->SetPosition(20, 105); // 20, 100
			ZIPTextActor->GetTextProperty()->SetFontFamilyToArial();
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

C_Actors LoadNC(const string& FilePath, double sphereRadius, double toolOffset) {
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
		vtkNew<vtkPolyData> ncPolyData;
		ncPolyData->SetPoints(ncPoints);
		vtkNew<vtkVertexGlyphFilter> ncGlyphFilter;
		ncGlyphFilter->SetInputData(ncPolyData);
		ncGlyphFilter->Update();

		vtkNew<vtkSphereSource> ncSphereSource;
		ncSphereSource->SetRadius(sphereRadius);

		vtkNew<vtkGlyph3D> ncGlyph3D;
		ncGlyph3D->SetSourceConnection(ncSphereSource->GetOutputPort());
		ncGlyph3D->SetInputConnection(ncGlyphFilter->GetOutputPort());
		ncGlyph3D->ScalingOff();

		vtkNew<vtkPolyDataMapper> ncMapper;
		ncMapper->SetInputConnection(ncGlyph3D->GetOutputPort());

		vtkNew<vtkActor> ncActor;
		ncActor->SetMapper(ncMapper);
		ncActor->GetProperty()->SetColor(NCCOLOR);

		// Model NC Text
		vtkNew<vtkTextActor> NCTextActor;
		NCTextActor->GetTextProperty()->SetFontSize(30);
		NCTextActor->GetTextProperty()->SetColor(TEXTCOLOR);
		NCTextActor->SetPosition(20, 1);
		NCTextActor->GetTextProperty()->SetFontFamilyToArial();
		NCTextActor->SetInput(std::format("Tool Diameter: {:.2f} mm\nNC: {} -> {} Points",
			sphereRadius * 2.0, path(FilePath).filename().string(), ncPoints->GetNumberOfPoints()).c_str());

		cout << CYAN << "          >> NC: " << END << path(FilePath).filename().string() << endl;
		return { ncActor, NCTextActor };
	}
	else {
		return { nullptr, nullptr };
	}
}

C_Actors LoadPTS(const string& FilePath, double sphereRadius) {
	if (!exists(FilePath)) return { nullptr, nullptr };
	vtkNew<vtkPoints> ptsPoints;
	vtkNew<vtkPolyData> ptsPolyData;
	vtkNew<vtkVertexGlyphFilter> ptsGlyphFilter;

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
		ptsPolyData->SetPoints(ptsPoints);
		ptsGlyphFilter->SetInputData(ptsPolyData);
		ptsGlyphFilter->Update();

		vtkNew<vtkSphereSource> sphereSource;
		sphereSource->SetRadius(sphereRadius);

		vtkNew<vtkGlyph3D> glyph3D;
		glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
		glyph3D->SetInputConnection(ptsGlyphFilter->GetOutputPort());
		glyph3D->ScalingOff();

		vtkNew<vtkPolyDataMapper> ptsMapper;
		ptsMapper->SetInputConnection(glyph3D->GetOutputPort());

		vtkNew<vtkActor> ptsActor;
		ptsActor->SetMapper(ptsMapper);
		ptsActor->GetProperty()->SetColor(PTSCOLOR);

		// Model PTS Text
		vtkNew<vtkTextActor> PTSTextActor;
		PTSTextActor->GetTextProperty()->SetFontSize(30);
		PTSTextActor->GetTextProperty()->SetColor(TEXTCOLOR);
		PTSTextActor->SetPosition(20, 1);
		PTSTextActor->GetTextProperty()->SetFontFamilyToArial();
		PTSTextActor->SetInput(std::format("PTS: {} -> {} Points", path(FilePath).filename().string(), ptsPoints->GetNumberOfPoints()).c_str());

		cout << CYAN << "          >> PTS: " << END << path(FilePath).filename().string() << endl;
		return { ptsActor, PTSTextActor };
	}
	else {
		return { nullptr, nullptr };
	}
}

C_Actors LoadSTL(const string& FilePath) {
	if (!exists(FilePath)) return { nullptr, nullptr };

	vtkNew<vtkSTLReader> stlReader;
	stlReader->SetFileName(FilePath.c_str());
	stlReader->Update();

	vtkNew<vtkPolyDataMapper> stlMapper;
	stlMapper->SetInputConnection(stlReader->GetOutputPort());

	vtkNew<vtkActor> stlActor;
	stlActor->SetMapper(stlMapper);
	stlActor->GetProperty()->SetColor(STLCOLOR);
	double bounds[6];
	stlActor->GetBounds(bounds);
	double Height_ = bounds[5] - bounds[4];
	vtkNew<vtkMassProperties> massProperties;
	massProperties->SetInputConnection(stlReader->GetOutputPort());
	massProperties->Update();
	double volumeMilliliters = massProperties->GetVolume() / 1000.0;

	// Model Text
	vtkNew<vtkTextActor> STLTextActor;
	STLTextActor->GetTextProperty()->SetFontSize(30);
	STLTextActor->GetTextProperty()->SetColor(TEXTCOLOR);
	STLTextActor->SetPosition(20, 15);
	STLTextActor->GetTextProperty()->SetFontFamilyToArial();
	STLTextActor->SetInput(std::format("ID: {}\nHeight: {:.2f} mm\nVolume: {:.2f} mL",
		path(FilePath).filename().string(), Height_, volumeMilliliters).c_str());

	cout << CYAN << "         >> STL: " << END << path(FilePath).filename().string() << endl;
	return { stlActor, STLTextActor };
}

void AddTextLabel(vtkRenderer* renderer, const string& text, int X, int Y, int fontsize, bool B, bool I) {
	vtkNew<vtkTextActor> textActor;
	textActor->SetInput(text.c_str());
	textActor->SetPosition(X, Y);
	textActor->GetTextProperty()->SetFontFamilyToArial();
	textActor->GetTextProperty()->SetFontSize(fontsize);
	textActor->GetTextProperty()->SetColor(TEXTCOLOR);
	textActor->GetTextProperty()->SetBold(B);
	textActor->GetTextProperty()->SetItalic(I);
	renderer->AddActor2D(textActor);
}

void ArrangeViewsInRenderWindow(vtkRenderer* renderer, const string& FilePath, string& CaseID, string& filter) {
	vector<string> views = { "Left", "Right", "Front", "Back" };
	size_t numViews = views.size();
	size_t cols = 2;
	size_t rows = (numViews + cols - 1) / cols;

	renderer->RemoveAllViewProps();

	C_Actors stlActors = { nullptr, nullptr };
	C_Actors ptsActors = { nullptr, nullptr };
	C_Actors ncActors = { nullptr, nullptr };
	C_Actors svgActors = { nullptr, nullptr };

	string stlPath = FilePath + ".stl",
		ptsPath = FilePath + ".pts",
		ncPath = FilePath + ".NC",
		svgPath = FilePath;

	stlActors = LoadSTL(stlPath);

	if (exists(ptsPath))
		ptsActors = LoadPTS(ptsPath, ToolDiameter);
	else
		ncActors = LoadNC(ncPath, ToolDiameter, ToolOffsetLength);

	svgActors = LoadLaserMarkingZip(svgPath, CaseID);


	for (int i = 0; i < numViews; ++i) {
		double x0 = (i % cols) * 0.5;
		double y0 = (rows - 1 - i / cols) * 0.5;
		double x1 = x0 + 0.5;
		double y1 = y0 + 0.5;

		vtkNew<vtkRenderer> subRenderer;
		subRenderer->SetViewport(x0, y0, x1, y1);
		subRenderer->SetBackground(BACKGROUND);

		if (ncActors.Actor) subRenderer->AddActor(ncActors.Actor);
		if (i == 0 && ncActors.textActor) subRenderer->AddActor(ncActors.textActor);
		if (ptsActors.Actor) subRenderer->AddActor(ptsActors.Actor);
		if (i == 0 && ptsActors.textActor) subRenderer->AddActor(ptsActors.textActor);
		if (stlActors.Actor) subRenderer->AddActor(stlActors.Actor);
		if (i == 2 && stlActors.textActor) subRenderer->AddActor(stlActors.textActor);
		if (svgActors.Actor) subRenderer->AddActor(svgActors.Actor);
		if (i == 2 && svgActors.textActor) subRenderer->AddActor(svgActors.textActor);

		// Get the size of the viewport
		int* windowSize;
		windowSize = renderer->GetRenderWindow()->GetSize();
		int VPWidth = static_cast<int>((x1 - x0) * windowSize[0]);
		int VPHeight = static_cast<int>((y1 - y0) * windowSize[1]);

		if (i == 0) {
			AddTextLabel(subRenderer, CaseID, 20, VPHeight - 70, 50, true, false);
			AddTextLabel(subRenderer, filter == "U" ? "UPPER" : "LOWER", 20, VPHeight - 110, 40, false, true);
		}
		AddTextLabel(subRenderer, views[i] + " View", VPWidth - 200, VPHeight - 50, 35, true, true);

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

void GetImageData(vector<string>& filenames, string& CaseID, string filter, string& FilesPath, vector<unsigned char>& ImageData, bool Grayscale) {
	vector<string> filtered_files;
	copy_if(filenames.begin(), filenames.end(), back_inserter(filtered_files),
		[&filter](const string& name) {
			return name.find(filter) != string::npos;
		});
	sort(filtered_files.begin(), filtered_files.end());

	int lowestStep = INT_MAX;
	string fileWithLowestStep;

	smatch match;
	for (string& fileID : filtered_files) {
		if (regex_search(fileID, match, filenameRegex)) {
			string stepStr = match[3].str();
			if (stepStr.empty()) stepStr = "00";

			int stepNumber = stoi(stepStr);
			if (stepNumber < lowestStep) {
				lowestStep = stepNumber;
				fileWithLowestStep = fileID;
			}
		}
	}

	if (!fileWithLowestStep.empty()) {
		cout << "      > " << CYAN << "First step number : " << YELLOW << fileWithLowestStep << END << endl;

		vtkNew<vtkRenderer> mainRenderer;
		mainRenderer->SetUseShadows(true);
		mainRenderer->SetUseDepthPeeling(1);        // Enable depth peeling
		mainRenderer->SetMaximumNumberOfPeels(100);
		mainRenderer->SetOcclusionRatio(0.7);
		vtkNew<vtkRenderWindow> renderWindow;
		renderWindow->SetSize(1830, 1464);			// 2440, 1952   1830, 1464     1220, 976
		renderWindow->SetAlphaBitPlanes(1);			// Enable alpha channel if supported
		renderWindow->SetMultiSamples(4);			// anti-aliasing
		renderWindow->AddRenderer(mainRenderer);

		ArrangeViewsInRenderWindow(mainRenderer, FilesPath + fileWithLowestStep, CaseID, filter);

		renderWindow->OffScreenRenderingOn();
		renderWindow->Render();

		vtkNew<vtkWindowToImageFilter> windowToImageFilter;
		windowToImageFilter->SetInput(renderWindow);
		windowToImageFilter->Update();

		vtkNew<vtkJPEGWriter> imageWriter;
		imageWriter->WriteToMemoryOn();

		if (Grayscale) {
			vtkNew<vtkImageLuminance> luminanceFilter;
			luminanceFilter->SetInputConnection(windowToImageFilter->GetOutputPort());
			luminanceFilter->Update();
			imageWriter->SetInputConnection(luminanceFilter->GetOutputPort());
		}
		else {
			imageWriter->SetInputConnection(windowToImageFilter->GetOutputPort());
		}

		imageWriter->SetQuality(30);
		imageWriter->Write();

		vtkUnsignedCharArray* JpegData = imageWriter->GetResult();
		vector<unsigned char> buffer(JpegData->GetNumberOfTuples());
		memcpy(buffer.data(), JpegData->GetPointer(0), buffer.size());
		ImageData = buffer;
	}
	else {
		ImageData.clear();
	}
}

bool save_PDF() {
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
		std::cout << RED << "        No STL files in input folder\n" << END << std::endl;
		std::cout << "        Press enter to continue...";
		std::cin.get();
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
	//sort(SortedCaseIDs.begin(), SortedCaseIDs.end());
	sort(SortedCaseIDs.begin(), SortedCaseIDs.end(), [](const string& a, const string& b) {
		return stoi(a) < stoi(b);
		});

	vector<unsigned char> UPPER_ImageData, LOWER_ImageData;

	int Count = 0;
	for (auto& CaseID : SortedCaseIDs) {
		Count++;
		cout << (Count < 10 ? "   0" : "   ") + to_string(Count) << " > " << GREEN << "Processing CaseID: " << END << CaseID << endl;

		GetImageData(Cases[CaseID], CaseID, "U", Main_path, UPPER_ImageData, true);

		GetImageData(Cases[CaseID], CaseID, "L", Main_path, LOWER_ImageData, true);

		string OutPDFname = CaseID + "_TL_Views.pdf";
		ExportToPDF(Main_path + OutPDFname, UPPER_ImageData, LOWER_ImageData);

		cout << YELLOW << "     << File Saved: " << END << OutPDFname << endl;

		//ShellExecuteA(NULL, "open", (Main_path + OutPDFname).c_str(), NULL, NULL, SW_SHOWNORMAL);
	}
	return true;
}
