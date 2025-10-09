#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
#include <limits>
#include <utility>      // For std::pair
#include <cstring>      // For std::memcpy
#include <iomanip>      // For std::setw, std::setfill

// VTK Headers
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkIdList.h>
#include <vtkCellTypes.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkTextProperty.h>
#include <vtkCommand.h>     // 用于动画回调
#include <vtkTextActor.h>     // 用于显示文本
#include <vtkUnstructuredGridReader.h>

// ===============================================
// 1. 自定义数据结构
// ===============================================
struct Cell {
    int vtk_type;
    std::vector<long long> point_ids;
};

struct InternalMesh {
    long long numberOfPoints = 0;
    std::vector<double> points; // 使用 double 保证精度
    long long numberOfCells = 0;
    std::vector<Cell> cells;
    std::map<std::string, std::vector<double>> pointData;
    std::map<std::string, std::vector<double>> cellData;
    std::map<std::string, double> fieldData;

    // A helper method to reset the struct
    void clear() {
        numberOfPoints = 0;
        points.clear();
        numberOfCells = 0;
        cells.clear();
        pointData.clear();
        cellData.clear();
        fieldData.clear();
    }
};

// ===============================================
// 2. VtkReader 类
// ===============================================
class VtkReader {
public:
    bool readFile(const std::string& filePath);
    const InternalMesh& getMeshData() const { return m_meshData; }

private:
    InternalMesh m_meshData;

    std::string toUpper(std::string s);
    std::vector<std::string> Stringsplit(const std::string& str);
    bool parseLine(const std::string& line, std::ifstream& file);
    bool readPoints(std::ifstream& in, int numPoints);
    bool readCells(std::ifstream& in, int numCells);
    bool readCellTypes(std::ifstream& in, int numCellTypes);
    bool readPointData(std::ifstream& in, int numPoints);
    bool readCellData(std::ifstream& in, int numCells);
};

bool VtkReader::readFile(const std::string& filePath) {
    // Clear previous data to make the reader instance reusable
    m_meshData.clear();

    std::ifstream file(filePath);
    if (!file.is_open()) { std::cerr << "Error opening file: " << filePath << std::endl; return false; }
    std::string line;
    while (std::getline(file, line)) {
        if (!parseLine(line, file)) { std::cerr << "Failed to parse line: " << line << std::endl; return false; }
    }
    file.close();
    return true;
}

// MODIFIED: Simplified toUpper function using a basic for loop
std::string VtkReader::toUpper(std::string s) {
    for (size_t i = 0; i < s.length(); ++i) {
        s[i] = std::toupper(s[i]);
    }
    return s;
}

std::vector<std::string> VtkReader::Stringsplit(const std::string& str) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string token;
    while (ss >> token) {
        result.push_back(token);
    }
    return result;
}

bool VtkReader::parseLine(const std::string& line, std::ifstream& file) {
    if (line.empty() || line[0] == '#') { return true; }
    std::vector<std::string> parts = Stringsplit(line);
    if (parts.empty()) { return true; }
    std::string keyword = toUpper(parts[0]);
    try {
        if (keyword == "POINTS") {
            return readPoints(file, std::stoi(parts[1]));
        }
        else if (keyword == "CELLS") {
            return readCells(file, std::stoi(parts[1]));
        }
        else if (keyword == "CELL_TYPES") {
            return readCellTypes(file, std::stoi(parts[1]));
        }
        else if (keyword == "POINT_DATA") {
            return readPointData(file, std::stoi(parts[1]));
        }
        else if (keyword == "CELL_DATA") {
            return readCellData(file, std::stoi(parts[1]));
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error during parsing keyword line: " << line << ". Details: " << e.what() << std::endl;
        return false;
    }
    return true;
}

// MODIFIED: Simplified the reading logic to be more direct.
bool VtkReader::readPoints(std::ifstream& in, int numPoints) {
    m_meshData.numberOfPoints = numPoints;
    m_meshData.points.assign(numPoints * 3, 0.0);
    for (int i = 0; i < numPoints; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (!(in >> m_meshData.points[i * 3 + j])) {
                std::cerr << "Error: Failed to parse coordinate for point " << i << std::endl;
                return false;
            }
        }
    }
    return true;
}

bool VtkReader::readCells(std::ifstream& in, int numCells) {
    m_meshData.numberOfCells = numCells;
    m_meshData.cells.resize(numCells);
    for (int i = 0; i < numCells; ++i) {
        int num_points_in_cell;
        if (!(in >> num_points_in_cell)) {
            std::cerr << "Error: Failed to parse cell structure for cell " << i << std::endl;
            return false;
        }
        m_meshData.cells[i].point_ids.resize(num_points_in_cell);
        for (int j = 0; j < num_points_in_cell; ++j) {
            if (!(in >> m_meshData.cells[i].point_ids[j])) {
                std::cerr << "Error reading point_id for cell " << i << std::endl;
                return false;
            }
        }
    }
    return true;
}

bool VtkReader::readCellTypes(std::ifstream& in, int numCellTypes) {
    if (m_meshData.cells.size() != static_cast<size_t>(numCellTypes)) { std::cerr << "Error: Mismatch between cell count and cell type count." << std::endl; return false; }
    for (int i = 0; i < numCellTypes; ++i) {
        if (!(in >> m_meshData.cells[i].vtk_type)) {
            std::cerr << "Error reading vtk_type for cell " << i << std::endl;
            return false;
        }
    }
    return true;
}

// NOTE: This reader is simplified to handle the "SCALARS" format.
// It does not handle the "FIELD" format from your previous version.
bool VtkReader::readPointData(std::ifstream& in, int numPoints) {
    std::string line;
    // Consume potential empty lines until a keyword is found
    while (std::getline(in, line) && line.find_first_not_of(" \t\r\n") == std::string::npos);

    std::stringstream ss(line);
    std::string keyword, dataName, dataType;
    ss >> keyword >> dataName >> dataType;

    if (toUpper(keyword) == "SCALARS") {
        std::getline(in, line); // Skip "LOOKUP_TABLE default" line
        std::vector<double>& values = m_meshData.pointData[dataName];
        values.resize(numPoints);
        for (int i = 0; i < numPoints; ++i) {
            if (!(in >> values[i])) {
                std::cerr << "Error reading scalar value for point " << i << std::endl;
                return false;
            }
        }
        return true;
    }
    std::cerr << "Unsupported POINT_DATA format. Only SCALARS is supported." << std::endl;
    return false;
}

bool VtkReader::readCellData(std::ifstream& in, int numCells) {
    std::string line;
    while (std::getline(in, line) && line.find_first_not_of(" \t\r\n") == std::string::npos);

    std::stringstream ss(line);
    std::string keyword, dataName, dataType;
    ss >> keyword >> dataName >> dataType;
    if (toUpper(keyword) == "SCALARS") {
        std::getline(in, line); // Skip "LOOKUP_TABLE default" line
        std::vector<double>& values = m_meshData.cellData[dataName];
        values.resize(numCells);
        for (int i = 0; i < numCells; ++i) {
            if (!(in >> values[i])) {
                std::cerr << "Error reading scalar value for cell " << i << std::endl;
                return false;
            }
        }
        return true;
    }
    std::cerr << "Unsupported CELL_DATA format. Only SCALARS is supported." << std::endl;
    return false;
}

// ===============================================
// 3. 辅助函数: 数据转换 (No changes needed here)
// ===============================================
vtkSmartPointer<vtkUnstructuredGrid> convertToVtkGrid(const InternalMesh& mesh) {
    auto vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(mesh.numberOfPoints);
    if (mesh.numberOfPoints > 0) {
        std::memcpy(points->GetVoidPointer(0), mesh.points.data(), mesh.numberOfPoints * 3 * sizeof(double));
    }
    vtkGrid->SetPoints(points);
    if (mesh.numberOfCells > 0) {
        vtkGrid->Allocate(mesh.numberOfCells);
        for (const auto& cell : mesh.cells) {
            vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
            for (long long id : cell.point_ids) {
                idList->InsertNextId(static_cast<vtkIdType>(id));
            }
            vtkGrid->InsertNextCell(cell.vtk_type, idList);
        }
    }
    for (const auto& pair : mesh.pointData) {
        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(pair.first.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(pair.second.size());
        if (!pair.second.empty()) {
            std::memcpy(dataArray->GetVoidPointer(0), pair.second.data(), pair.second.size() * sizeof(double));
        }
        vtkGrid->GetPointData()->AddArray(dataArray);
    }
    for (const auto& pair : mesh.cellData) {
        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(pair.first.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(pair.second.size());
        if (!pair.second.empty()) {
            std::memcpy(dataArray->GetVoidPointer(0), pair.second.data(), pair.second.size() * sizeof(double));
        }
        vtkGrid->GetCellData()->AddArray(dataArray);
    }
    return vtkGrid;
}

// ===============================================
// 4. 动画与可视化相关 (No changes needed here)
// ===============================================
class AnimationTimerCallback : public vtkCommand {
public:
    static AnimationTimerCallback* New() { return new AnimationTimerCallback; }
    void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
        this->CurrentFrameIndex = (this->CurrentFrameIndex + 1) % this->Grids->size();
        vtkDataSet* currentGrid = (*this->Grids)[this->CurrentFrameIndex];
        this->Mapper->SetInputData(currentGrid);
        this->TextActor->SetInput((*this->FileNames)[this->CurrentFrameIndex].c_str());
        this->RenderWindow->Render();
    }
    void SetData(
        std::vector<vtkSmartPointer<vtkDataSet>>* grids,
        const std::vector<std::string>* fileNames,
        vtkDataSetMapper* mapper,
        vtkTextActor* textActor,
        vtkRenderWindow* renderWindow) {
        this->Grids = grids;
        this->FileNames = fileNames;
        this->Mapper = mapper;
        this->TextActor = textActor;
        this->RenderWindow = renderWindow;
    }
private:
    int CurrentFrameIndex = 0;
    std::vector<vtkSmartPointer<vtkDataSet>>* Grids = nullptr;
    const std::vector<std::string>* FileNames = nullptr;
    vtkDataSetMapper* Mapper = nullptr;
    vtkTextActor* TextActor = nullptr;
    vtkRenderWindow* RenderWindow = nullptr;
};

void animateGrids(
    std::vector<vtkSmartPointer<vtkDataSet>>& grids,
    const std::vector<std::string>& fileNames,
    double globalMin, double globalMax,
    const char* scalarBarTitle) {
    if (grids.empty()) {
        std::cout << "No grids to animate." << std::endl;
        return;
    }
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0.1, 0.2, 0.4);
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(1000, 800);
    renderWindow->SetWindowName("VTK Animation Player");
    auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);
    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(grids[0]);
    mapper->SetScalarRange(globalMin, globalMax);
    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->AddRGBPoint(globalMin, 0.2, 0.2, 1.0);
    ctf->AddRGBPoint(0.0, 1.0, 1.0, 1.0);
    ctf->AddRGBPoint(globalMax, 1.0, 0.2, 0.2);
    mapper->SetLookupTable(ctf);
    if (grids[0]->GetPointData()->GetScalars()) {
        mapper->SetScalarModeToUsePointData();
        mapper->SelectColorArray(grids[0]->GetPointData()->GetScalars()->GetName());
    }
    else if (grids[0]->GetCellData()->GetScalars()) {
        mapper->SetScalarModeToUseCellData();
        mapper->SelectColorArray(grids[0]->GetCellData()->GetScalars()->GetName());
    }
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    renderer->AddActor(actor);
    auto scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle(scalarBarTitle);
    scalarBar->SetNumberOfLabels(5);
    renderer->AddActor2D(scalarBar);
    auto textActor = vtkSmartPointer<vtkTextActor>::New();
    textActor->SetInput(fileNames[0].c_str());
    textActor->GetTextProperty()->SetFontSize(16);
    textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    textActor->SetPosition(20, 20);
    renderer->AddActor2D(textActor);
    auto callback = vtkSmartPointer<AnimationTimerCallback>::New();
    callback->SetData(&grids, &fileNames, mapper, textActor, renderWindow);
    interactor->AddObserver(vtkCommand::TimerEvent, callback);
    interactor->Initialize();
    interactor->CreateRepeatingTimer(100);
    std::cout << "\nStarting animation... Close the window to stop." << std::endl;
    renderWindow->Render();
    interactor->Start();
}

std::pair<double, double> calculateGlobalScalarRange(const std::vector<vtkSmartPointer<vtkDataSet>>& grids, std::string& determinedArrayName) {
    double globalMin = std::numeric_limits<double>::max();
    double globalMax = std::numeric_limits<double>::lowest();
    if (grids.empty()) return { 0.0, 0.0 };
    const char* arrayName = nullptr;
    if (grids[0]->GetPointData()->GetNumberOfArrays() > 0) {
        arrayName = grids[0]->GetPointData()->GetArray(0)->GetName();
        grids[0]->GetPointData()->SetActiveScalars(arrayName);
    }
    else if (grids[0]->GetCellData()->GetNumberOfArrays() > 0) {
        arrayName = grids[0]->GetCellData()->GetArray(0)->GetName();
        grids[0]->GetCellData()->SetActiveScalars(arrayName);
    }
    else {
        return { 0.0, 0.0 };
    }
    determinedArrayName = arrayName;
    std::cout << "\n--- Calculating Global Scalar Range from loaded data ('" << determinedArrayName << "') ---" << std::endl;
    for (size_t i = 0; i < grids.size(); ++i) {
        auto& grid = grids[i];
        vtkDataArray* dataArray = grid->GetPointData()->GetArray(arrayName);
        if (!dataArray) dataArray = grid->GetCellData()->GetArray(arrayName);
        if (dataArray) {
            if (grid->GetPointData()->GetArray(arrayName)) grid->GetPointData()->SetActiveScalars(arrayName);
            else if (grid->GetCellData()->GetArray(arrayName)) grid->GetCellData()->SetActiveScalars(arrayName);
            double range[2];
            dataArray->GetRange(range);
            globalMin = std::min(globalMin, range[0]);
            globalMax = std::max(globalMax, range[1]);
        }
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    return { globalMin, globalMax };
}

// MODIFIED: New helper functions to replace <filesystem>
// =======================================================
// This function generates a list of file paths based on a pattern.
std::vector<std::string> generateFilePaths(
    const std::string& basePath,
    const std::string& prefix,
    const std::string& extension,
    int startIndex,
    int endIndex,
    int padding) {

    std::vector<std::string> paths;
    for (int i = startIndex; i <= endIndex; ++i) {
        std::stringstream ss;
        ss << basePath << prefix
            << std::setw(padding) << std::setfill('0') << i
            << extension;
        paths.push_back(ss.str());
    }
    return paths;
}

// This function extracts the filename from a full path.
std::string getShortFileName(const std::string& filePath) {
    size_t lastSlash = filePath.find_last_of("/\\");
    if (lastSlash != std::string::npos) {
        return filePath.substr(lastSlash + 1);
    }
    return filePath;
}
// =======================================================


// ===============================================
// 5. Main Function
// ===============================================
int main() {
    // MODIFIED: Removed <filesystem> and replaced with a pattern-based approach.
    // --- Step 1: Define the list of files to read ---
    // Please configure these parameters to match your file sequence.
    std::string directoryPath = "D:/data/";
    std::string filePrefix = "vtkout";
    std::string fileExtension = ".vtk";
    int startIndex = 1054;
    int endIndex = 1060;
    int numberPadding = 6; // e.g., for "001054", padding is 6.

    std::vector<std::string> vtkFilePaths = generateFilePaths(
        directoryPath, filePrefix, fileExtension, startIndex, endIndex, numberPadding
    );

    if (vtkFilePaths.empty()) {
        std::cout << "No file paths were generated. Please check the parameters in main()." << std::endl;
        return 0;
    }
    std::cout << "Generated " << vtkFilePaths.size() << " VTK file paths to read." << std::endl;

    // --- Step 2: Read all files into memory ---
    std::vector<vtkSmartPointer<vtkDataSet>> allGrids;
    std::vector<std::string> shortFileNames;
    std::cout << "\n--- Loading all files into memory using custom VtkReader ---" << std::endl;

    for (const auto& filePath : vtkFilePaths) {
        // MODIFIED: Use the new helper function to get the short filename
        std::string shortName = getShortFileName(filePath);
        std::cout << "Reading: " << shortName << "..." << std::endl;

        VtkReader myReader;
        if (myReader.readFile(filePath)) {
            const InternalMesh& meshData = myReader.getMeshData();
            vtkSmartPointer<vtkUnstructuredGrid> grid = convertToVtkGrid(meshData);
            if (grid && grid->GetNumberOfPoints() > 0) {
                allGrids.push_back(grid);
                shortFileNames.push_back(shortName);
            }
            else {
                std::cerr << "Failed to convert custom mesh data to VTK grid, or grid is empty: " << filePath << std::endl;
            }
        }
        else {
            std::cerr << "Failed to read file using custom VtkReader: " << filePath << std::endl;
        }
    }

    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Successfully loaded " << allGrids.size() << " datasets." << std::endl;

    // --- Step 3: Calculate global scalar range (no changes) ---
    std::string scalarBarTitle;
    std::pair<double, double> globalRange = calculateGlobalScalarRange(allGrids, scalarBarTitle);
    if (globalRange.first > globalRange.second) {
        std::cout << "\nCould not determine a valid global scalar range. Animation may not color correctly." << std::endl;
        globalRange.first = 0.0;
        globalRange.second = 1.0;
        scalarBarTitle = "N/A";
    }
    else {
        std::cout << "\nGlobal Scalar Range calculated: ["
            << globalRange.first << ", " << globalRange.second << "]" << std::endl;
    }

    // --- Step 4: Call animation function (no changes) ---
    animateGrids(allGrids, shortFileNames, globalRange.first, globalRange.second, scalarBarTitle.c_str());

    std::cout << "\nAnimation window closed. Exiting program." << std::endl;
    return 0;
}
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <algorithm>
//#include <sstream>
//#include <map>
//#include <limits>
//#include <filesystem> // C++17: 用于文件系统操作
//#include <utility>    // 用于 std::pair
//#include <cstring>    // 用于 std::memcpy
//
// VTK Headers
//#include <vtkSmartPointer.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkDataSet.h>
//#include <vtkPoints.h>
//#include <vtkCellArray.h>
//#include <vtkDoubleArray.h>
//#include <vtkPointData.h>
//#include <vtkCellData.h>
//#include <vtkDataSetMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkProperty.h>
//#include <vtkIdList.h>
//#include <vtkCellTypes.h>
//#include <vtkScalarBarActor.h>
//#include <vtkLookupTable.h>
//#include <vtkColorTransferFunction.h>
//#include <vtkTextProperty.h>
//#include <vtkCommand.h>     // 用于动画回调
//#include <vtkTextActor.h>     // 用于显示文本
//#include <vtkUnstructuredGridReader.h>
//
//
// ===============================================
// 1. 自定义数据结构
// ===============================================
//struct Cell {
//    int vtk_type;
//    std::vector<long long> point_ids;
//};
//
//struct InternalMesh {
//    long long numberOfPoints = 0;
//    std::vector<double> points; // 使用 double 保证精度
//    long long numberOfCells = 0;
//    std::vector<Cell> cells;
//    std::map<std::string, std::vector<double>> pointData;
//    std::map<std::string, std::vector<double>> cellData;
//    std::map<std::string, double> fieldData;
//};
//
// ===============================================
// 2. VtkReader 类
// ===============================================
//class VtkReader {
//public:
//    bool readFile(const std::string& filePath);
//    const InternalMesh& getMeshData() const { return m_meshData; }
//
//private:
//    InternalMesh m_meshData;
//
//    std::string toUpper(std::string s);
//    std::vector<std::string> Stringsplit(const std::string& str);
//    bool parseLine(const std::string& line, std::ifstream& file);
//    bool readPoints(std::ifstream& in, int numPoints);
//    bool readCells(std::ifstream& in, int numCells);
//    bool readCellTypes(std::ifstream& in, int numCellTypes);
//    bool readPointData(std::ifstream& in, int numPoints);
//    bool readCellData(std::ifstream& in, int numCells);
//};
//
//bool VtkReader::readFile(const std::string& filePath) {
//    std::ifstream file(filePath);
//    if (!file.is_open()) { std::cerr << "Error opening file: " << filePath << std::endl; return false; }
//    std::string line;
//    while (std::getline(file, line)) {
//        if (!parseLine(line, file)) { std::cerr << "Failed to parse line: " << line << std::endl; return false; }
//    }
//    file.close();
//    return true;
//}
//
//std::string VtkReader::toUpper(std::string s) {
//    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
//    return s;
//}
//
//std::vector<std::string> VtkReader::Stringsplit(const std::string& str) {
//    std::vector<std::string> result;
//    std::stringstream ss(str);
//    std::string token;
//    while (ss >> token) {
//        result.push_back(token);
//    }
//    return result;
//}
//
//bool VtkReader::parseLine(const std::string& line, std::ifstream& file) {
//    if (line.empty() || line[0] == '#') { return true; }
//    std::vector<std::string> parts = Stringsplit(line);
//    if (parts.empty()) { return true; }
//    std::string keyword = toUpper(parts[0]);
//    try {
//        if (keyword == "POINTS") {
//            return readPoints(file, std::stoi(parts[1]));
//        }
//        else if (keyword == "CELLS") {
//            return readCells(file, std::stoi(parts[1]));
//        }
//        else if (keyword == "CELL_TYPES") {
//            return readCellTypes(file, std::stoi(parts[1]));
//        }
//        else if (keyword == "POINT_DATA") {
//            return readPointData(file, std::stoi(parts[1]));
//        }
//        else if (keyword == "CELL_DATA") {
//            return readCellData(file, std::stoi(parts[1]));
//        }
//    }
//    catch (const std::exception& e) {
//        std::cerr << "Error during parsing keyword line: " << line << ". Details: " << e.what() << std::endl;
//        return false;
//    }
//    return true;
//}
//
//bool VtkReader::readPoints(std::ifstream& in, int numPoints) {
//    m_meshData.numberOfPoints = numPoints;
//    m_meshData.points.assign(numPoints * 3, 0.0);
//    for (int i = 0; i < numPoints; ++i) {
//        if (!(in >> m_meshData.points[i * 3] >> m_meshData.points[i * 3 + 1] >> m_meshData.points[i * 3 + 2])) {
//            std::cerr << "Error: Failed to parse point coordinates for point " << i << std::endl;
//            return false;
//        }
//    }
//    return true;
//}
//
//bool VtkReader::readCells(std::ifstream& in, int numCells) {
//    m_meshData.numberOfCells = numCells;
//    m_meshData.cells.resize(numCells);
//    for (int i = 0; i < numCells; ++i) {
//        int num_points_in_cell;
//        if (!(in >> num_points_in_cell)) {
//            std::cerr << "Error: Failed to parse cell structure for cell " << i << std::endl;
//            return false;
//        }
//        m_meshData.cells[i].point_ids.resize(num_points_in_cell);
//        for (int j = 0; j < num_points_in_cell; ++j) {
//            if (!(in >> m_meshData.cells[i].point_ids[j])) {
//                std::cerr << "Error reading point_id for cell " << i << std::endl;
//                return false;
//            }
//        }
//    }
//    return true;
//}
//
//bool VtkReader::readCellTypes(std::ifstream& in, int numCellTypes) {
//    if (m_meshData.cells.size() != numCellTypes) { std::cerr << "Error: Mismatch between cell count and cell type count." << std::endl; return false; }
//    for (int i = 0; i < numCellTypes; ++i) {
//        if (!(in >> m_meshData.cells[i].vtk_type)) {
//            std::cerr << "Error reading vtk_type for cell " << i << std::endl;
//            return false;
//        }
//    }
//    return true;
//}
//
//bool VtkReader::readPointData(std::ifstream& in, int numPoints) {
//     This is a simplified reader for SCALARS data, a common format.
//     The FIELD format is more complex and might require adjustments.
//    std::string line;
//    std::getline(in, line); // Consume the rest of the POINT_DATA line
//
//     Look for SCALARS definition line
//    while (std::getline(in, line)) {
//        if (line.empty()) continue;
//        std::stringstream ss(line);
//        std::string keyword, dataName, dataType;
//        ss >> keyword >> dataName >> dataType;
//        if (toUpper(keyword) == "SCALARS") {
//             Skip lookup table line
//            std::getline(in, line);
//            std::vector<double> values;
//            values.reserve(numPoints);
//            for (int i = 0; i < numPoints; ++i) {
//                double val;
//                if (!(in >> val)) {
//                    std::cerr << "Error reading scalar value for point " << i << std::endl;
//                    return false;
//                }
//                values.push_back(val);
//            }
//            m_meshData.pointData[dataName] = values;
//            return true;
//        }
//    }
//    return false; // Or handle other data types like VECTORS etc.
//}
//
//
//bool VtkReader::readCellData(std::ifstream& in, int numCells) {
//     Similar simplified reader for cell data
//    std::string line;
//    std::getline(in, line);
//
//    while (std::getline(in, line)) {
//        if (line.empty()) continue;
//        std::stringstream ss(line);
//        std::string keyword, dataName, dataType;
//        ss >> keyword >> dataName >> dataType;
//        if (toUpper(keyword) == "SCALARS") {
//            std::getline(in, line);
//            std::vector<double> values;
//            values.reserve(numCells);
//            for (int i = 0; i < numCells; ++i) {
//                double val;
//                if (!(in >> val)) {
//                    std::cerr << "Error reading scalar value for cell " << i << std::endl;
//                    return false;
//                }
//                values.push_back(val);
//            }
//            m_meshData.cellData[dataName] = values;
//            return true;
//        }
//    }
//    return false;
//}
//
//
// ===============================================
// 3. 辅助函数: 数据转换
// ===============================================
//vtkSmartPointer<vtkUnstructuredGrid> convertToVtkGrid(const InternalMesh& mesh) {
//    auto vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
//
//     1. 填充节点坐标
//    auto points = vtkSmartPointer<vtkPoints>::New();
//    points->SetDataTypeToDouble();
//    points->SetNumberOfPoints(mesh.numberOfPoints);
//    if (mesh.numberOfPoints > 0) {
//        std::memcpy(points->GetVoidPointer(0), mesh.points.data(), mesh.numberOfPoints * 3 * sizeof(double));
//    }
//    vtkGrid->SetPoints(points);
//
//     2. 填充单元拓扑
//    if (mesh.numberOfCells > 0) {
//        vtkGrid->Allocate(mesh.numberOfCells);
//        for (const auto& cell : mesh.cells) {
//            vtkIdType numIds = cell.point_ids.size();
//            vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
//            for (long long id : cell.point_ids) {
//                idList->InsertNextId(static_cast<vtkIdType>(id));
//            }
//            vtkGrid->InsertNextCell(cell.vtk_type, idList);
//        }
//    }
//
//     3. 填充节点场数据
//    for (const auto& pair : mesh.pointData) {
//        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
//        dataArray->SetName(pair.first.c_str());
//        dataArray->SetNumberOfComponents(1);
//        dataArray->SetNumberOfTuples(pair.second.size());
//        if (!pair.second.empty()) {
//            std::memcpy(dataArray->GetVoidPointer(0), pair.second.data(), pair.second.size() * sizeof(double));
//        }
//        vtkGrid->GetPointData()->AddArray(dataArray);
//    }
//
//     4. 填充单元场数据
//    for (const auto& pair : mesh.cellData) {
//        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
//        dataArray->SetName(pair.first.c_str());
//        dataArray->SetNumberOfComponents(1);
//        dataArray->SetNumberOfTuples(pair.second.size());
//        if (!pair.second.empty()) {
//            std::memcpy(dataArray->GetVoidPointer(0), pair.second.data(), pair.second.size() * sizeof(double));
//        }
//        vtkGrid->GetCellData()->AddArray(dataArray);
//    }
//
//    return vtkGrid;
//}
//
// ===============================================
// 4. 动画与可视化相关
// ===============================================
//
//class AnimationTimerCallback : public vtkCommand {
//public:
//    static AnimationTimerCallback* New() {
//        return new AnimationTimerCallback;
//    }
//
//    void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
//        this->CurrentFrameIndex = (this->CurrentFrameIndex + 1) % this->Grids->size();
//        vtkDataSet* currentGrid = (*this->Grids)[this->CurrentFrameIndex];
//        this->Mapper->SetInputData(currentGrid);
//        this->TextActor->SetInput((*this->FileNames)[this->CurrentFrameIndex].c_str());
//        this->RenderWindow->Render();
//    }
//
//    void SetData(
//        std::vector<vtkSmartPointer<vtkDataSet>>* grids,
//        const std::vector<std::string>* fileNames,
//        vtkDataSetMapper* mapper,
//        vtkTextActor* textActor,
//        vtkRenderWindow* renderWindow) {
//        this->Grids = grids;
//        this->FileNames = fileNames;
//        this->Mapper = mapper;
//        this->TextActor = textActor;
//        this->RenderWindow = renderWindow;
//    }
//
//private:
//    int CurrentFrameIndex = 0;
//    std::vector<vtkSmartPointer<vtkDataSet>>* Grids = nullptr;
//    const std::vector<std::string>* FileNames = nullptr;
//    vtkDataSetMapper* Mapper = nullptr;
//    vtkTextActor* TextActor = nullptr;
//    vtkRenderWindow* RenderWindow = nullptr;
//};
//
//void animateGrids(
//    std::vector<vtkSmartPointer<vtkDataSet>>& grids,
//    const std::vector<std::string>& fileNames,
//    double globalMin, double globalMax,
//    const char* scalarBarTitle) {
//    if (grids.empty()) {
//        std::cout << "No grids to animate." << std::endl;
//        return;
//    }
//
//    auto renderer = vtkSmartPointer<vtkRenderer>::New();
//    renderer->SetBackground(0.1, 0.2, 0.4);
//    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//    renderWindow->SetSize(1000, 800);
//    renderWindow->SetWindowName("VTK Animation Player");
//    auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    interactor->SetRenderWindow(renderWindow);
//
//    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
//    mapper->SetInputData(grids[0]);
//    mapper->SetScalarRange(globalMin, globalMax);
//
//     =================== START OF MODIFICATION ===================
//     创建一个颜色传输函数 (Color Transfer Function)
//    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
//
//     设置颜色节点: AddRGBPoint(value, R, G, B)
//     R, G, B 的范围是 0.0 到 1.0
//
//     将最小值(globalMin)设置为蓝色 (0, 0, 1)
//    ctf->AddRGBPoint(globalMin, 0.2, 0.2, 1.0);
//
//     将 0.0 设置为白色 (1, 1, 1)，这是发散式色图的关键
//    ctf->AddRGBPoint(0.0, 1.0, 1.0, 1.0);
//
//     将最大值(globalMax)设置为红色 (1, 0, 0)
//    ctf->AddRGBPoint(globalMax, 1.0, 0.2, 0.2);
//
//     将这个颜色传输函数应用到Mapper
//    mapper->SetLookupTable(ctf);
//     =================== END OF MODIFICATION =====================
//
//    if (grids[0]->GetPointData()->GetScalars()) {
//        mapper->SetScalarModeToUsePointData();
//        mapper->SelectColorArray(grids[0]->GetPointData()->GetScalars()->GetName());
//    }
//    else if (grids[0]->GetCellData()->GetScalars()) {
//        mapper->SetScalarModeToUseCellData();
//        mapper->SelectColorArray(grids[0]->GetCellData()->GetScalars()->GetName());
//    }
//
//    auto actor = vtkSmartPointer<vtkActor>::New();
//    actor->SetMapper(mapper);
//    renderer->AddActor(actor);
//
//    auto scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
//    scalarBar->SetLookupTable(mapper->GetLookupTable());
//    scalarBar->SetTitle(scalarBarTitle);
//    scalarBar->SetNumberOfLabels(5);
//    renderer->AddActor2D(scalarBar);
//
//    auto textActor = vtkSmartPointer<vtkTextActor>::New();
//    textActor->SetInput(fileNames[0].c_str());
//    textActor->GetTextProperty()->SetFontSize(16);
//    textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
//    textActor->SetPosition(20, 20);
//    renderer->AddActor2D(textActor);
//
//    auto callback = vtkSmartPointer<AnimationTimerCallback>::New();
//    callback->SetData(&grids, &fileNames, mapper, textActor, renderWindow);
//
//    interactor->AddObserver(vtkCommand::TimerEvent, callback);
//    interactor->Initialize();
//    interactor->CreateRepeatingTimer(100);
//
//    std::cout << "\nStarting animation... Close the window to stop." << std::endl;
//    renderWindow->Render();
//    interactor->Start();
//}
//
//std::pair<double, double> calculateGlobalScalarRange(const std::vector<vtkSmartPointer<vtkDataSet>>& grids, std::string& determinedArrayName) {
//    double globalMin = std::numeric_limits<double>::max();
//    double globalMax = std::numeric_limits<double>::lowest();
//    if (grids.empty()) return { 0.0, 0.0 };
//
//    const char* arrayName = nullptr;
//    if (grids[0]->GetPointData()->GetNumberOfArrays() > 0) {
//        arrayName = grids[0]->GetPointData()->GetArray(0)->GetName();
//        grids[0]->GetPointData()->SetActiveScalars(arrayName);
//    }
//    else if (grids[0]->GetCellData()->GetNumberOfArrays() > 0) {
//        arrayName = grids[0]->GetCellData()->GetArray(0)->GetName();
//        grids[0]->GetCellData()->SetActiveScalars(arrayName);
//    }
//    else {
//        return { 0.0, 0.0 };
//    }
//    determinedArrayName = arrayName;
//    std::cout << "\n--- Calculating Global Scalar Range from loaded data ('" << determinedArrayName << "') ---" << std::endl;
//
//    for (size_t i = 0; i < grids.size(); ++i) {
//        auto& grid = grids[i];
//        vtkDataArray* dataArray = grid->GetPointData()->GetArray(arrayName);
//        if (!dataArray) dataArray = grid->GetCellData()->GetArray(arrayName);
//
//        if (dataArray) {
//             Set active scalar for each grid to ensure consistency
//            if (grid->GetPointData()->GetArray(arrayName)) grid->GetPointData()->SetActiveScalars(arrayName);
//            else if (grid->GetCellData()->GetArray(arrayName)) grid->GetCellData()->SetActiveScalars(arrayName);
//
//            double range[2];
//            dataArray->GetRange(range);
//            globalMin = std::min(globalMin, range[0]);
//            globalMax = std::max(globalMax, range[1]);
//        }
//    }
//    std::cout << "--------------------------------------------------------" << std::endl;
//    return { globalMin, globalMax };
//}
//
// ===============================================
// 5. Main Function
// ===============================================
//int main() {
//     1. 设置目录并查找文件 (这部分不变)
//    std::string directoryPath = "D:/data/";
//    std::vector<std::string> vtkFilePaths;
//    try {
//        for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
//            if (entry.is_regular_file() && entry.path().extension() == ".vtk") {
//                vtkFilePaths.push_back(entry.path().string());
//            }
//        }
//    }
//    catch (const std::filesystem::filesystem_error& e) {
//        std::cerr << "Error accessing directory " << directoryPath << ": " << e.what() << std::endl;
//        return 1;
//    }
//
//    if (vtkFilePaths.empty()) {
//        std::cout << "No .vtk files found in the specified directory: " << directoryPath << std::endl;
//        return 0;
//    }
//
//    std::sort(vtkFilePaths.begin(), vtkFilePaths.end());
//    std::cout << "Found " << vtkFilePaths.size() << " VTK files." << std::endl;
//
//
//     =================== START OF MODIFICATION ===================
//     2. 将所有文件读取并存储到 vector 中 (现在使用你自己的 VtkReader)
//    std::vector<vtkSmartPointer<vtkDataSet>> allGrids;
//    std::vector<std::string> shortFileNames;
//    std::cout << "\n--- Loading all files into memory using YOUR custom VtkReader ---" << std::endl;
//
//    for (const auto& filePath : vtkFilePaths) {
//        std::cout << "Reading: " << std::filesystem::path(filePath).filename().string() << "..." << std::endl;
//
//         1. 创建你自己的 VtkReader 实例
//        VtkReader myReader;
//
//         2. 调用 readFile 方法解析文件
//        if (myReader.readFile(filePath)) {
//             3. 读取成功后，获取解析出的 InternalMesh 数据
//            const InternalMesh& meshData = myReader.getMeshData();
//
//             4. 调用转换函数，将 InternalMesh 转换为 vtkUnstructuredGrid
//            vtkSmartPointer<vtkUnstructuredGrid> grid = convertToVtkGrid(meshData);
//
//             5. 确保转换后的 grid 有效，然后添加到列表中
//            if (grid && grid->GetNumberOfPoints() > 0) {
//                allGrids.push_back(grid);
//                shortFileNames.push_back(std::filesystem::path(filePath).filename().string());
//            }
//            else {
//                std::cerr << "Failed to convert custom mesh data to VTK grid, or grid is empty: " << filePath << std::endl;
//            }
//        }
//        else {
//            std::cerr << "Failed to read file using custom VtkReader: " << filePath << std::endl;
//        }
//    }
//     ===================  END OF MODIFICATION  ===================
//
//    std::cout << "-----------------------------------" << std::endl;
//    std::cout << "Successfully loaded " << allGrids.size() << " datasets." << std::endl;
//
//
//     3. 从已加载的数据中计算全局范围 (这部分代码无需改变)
//    std::string scalarBarTitle;
//    std::pair<double, double> globalRange = calculateGlobalScalarRange(allGrids, scalarBarTitle);
//    if (globalRange.first > globalRange.second) {
//        std::cout << "\nCould not determine a valid global scalar range. Animation may not color correctly." << std::endl;
//        globalRange.first = 0.0;
//        globalRange.second = 1.0;
//        scalarBarTitle = "N/A";
//    }
//    else {
//        std::cout << "\nGlobal Scalar Range calculated: ["
//            << globalRange.first << ", " << globalRange.second << "]" << std::endl;
//    }
//
//     4. 调用动画函数 (这部分代码无需改变)
//    animateGrids(allGrids, shortFileNames, globalRange.first, globalRange.second, scalarBarTitle.c_str());
//
//    std::cout << "\nAnimation window closed. Exiting program." << std::endl;
//    return 0;
//}

//int main() {
//    // 1. 设置目录并查找文件
//    std::string directoryPath = "D:/data/";
//    std::vector<std::string> vtkFilePaths;
//    try {
//        for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
//            if (entry.is_regular_file() && entry.path().extension() == ".vtk") {
//                vtkFilePaths.push_back(entry.path().string());
//            }
//        }
//    }
//    catch (const std::filesystem::filesystem_error& e) {
//        std::cerr << "Error accessing directory " << directoryPath << ": " << e.what() << std::endl;
//        return 1;
//    }
//
//    if (vtkFilePaths.empty()) {
//        std::cout << "No .vtk files found in the specified directory: " << directoryPath << std::endl;
//        return 0;
//    }
//
//    std::sort(vtkFilePaths.begin(), vtkFilePaths.end());
//    std::cout << "Found " << vtkFilePaths.size() << " VTK files." << std::endl;
//
//    // 2. 将所有文件读取并存储到 vector 中
//    std::vector<vtkSmartPointer<vtkDataSet>> allGrids;
//    std::vector<std::string> shortFileNames;
//    std::cout << "\n--- Loading all files into memory using vtkUnstructuredGridReader ---" << std::endl;
//
//    for (const auto& filePath : vtkFilePaths) {
//        std::cout << "Reading: " << std::filesystem::path(filePath).filename().string() << "..." << std::endl;
//
//        // 使用VTK官方Reader
//        auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
//        reader->SetFileName(filePath.c_str());
//        reader->Update(); // 执行读取
//
//        if (reader->GetOutput() && reader->GetOutput()->GetNumberOfPoints() > 0) {
//            // 注意：这里需要创建一个副本，因为reader的生命周期在循环内
//            auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
//            grid->DeepCopy(reader->GetOutput());
//
//            allGrids.push_back(grid);
//            shortFileNames.push_back(std::filesystem::path(filePath).filename().string());
//        }
//        else {
//            std::cerr << "Failed to read file or file is empty: " << filePath << std::endl;
//        }
//    }
//    std::cout << "-----------------------------------" << std::endl;
//    std::cout << "Successfully loaded " << allGrids.size() << " datasets." << std::endl;
//
//    // 3. 从已加载的数据中计算全局范围 (这部分代码无需改变)
//    std::string scalarBarTitle;
//    std::pair<double, double> globalRange = calculateGlobalScalarRange(allGrids, scalarBarTitle);
//    if (globalRange.first > globalRange.second) {
//        std::cout << "\nCould not determine a valid global scalar range. Animation may not color correctly." << std::endl;
//        globalRange.first = 0.0;
//        globalRange.second = 1.0;
//        scalarBarTitle = "N/A";
//    }
//    else {
//        std::cout << "\nGlobal Scalar Range calculated: ["
//            << globalRange.first << ", " << globalRange.second << "]" << std::endl;
//    }
//
//    // 4. 调用动画函数 (这部分代码无需改变)
//    animateGrids(allGrids, shortFileNames, globalRange.first, globalRange.second, scalarBarTitle.c_str());
//
//    std::cout << "\nAnimation window closed. Exiting program." << std::endl;
//    return 0;
//}