#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
#include <limits>
// #include <filesystem> // C++17: 已被移除
#include <utility>      // 用于 std::pair
#include <cstring>      // 用于 std::memcpy
#include <cctype>       // 用于 std::toupper

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
#include <vtkCommand.h>      // 用于动画回调
#include <vtkTextActor.h>      // 用于显示文本
#include <vtkUnstructuredGridReader.h>
// ...
using namespace std;

// ===============================================
// 1. 自定义数据结构 (无变化)
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
};

// ===============================================
// 2. VtkReader 类声明 (无变化)
// ===============================================
class VtkReader {
public:
    bool readFile(const std::string& filePath);
    const InternalMesh& getMeshData() const { return m_meshData; }

private:
    enum class FileFormat { ASCII, BINARY };
    InternalMesh m_meshData;
    FileFormat m_format = FileFormat::ASCII;

    std::string toUpper(std::string s);
    std::vector<std::string> Stringsplit(const std::string& str);
    bool parseLine(const std::string& line, std::ifstream& file);

    bool readPoints(std::ifstream& in, int numPoints);
    bool readCells(std::ifstream& in, int numCells);
    bool readCellTypes(std::ifstream& in, int numCellTypes);
    bool readPointData(std::ifstream& in, int numPoints);
    bool readCellData(std::ifstream& in, int numCells);
};

// ===============================================
// 3. 辅助函数 (打印与转换)
// ===============================================

void printMeshInfo(const InternalMesh& mesh) {
    cout << "\n--- Custom Reader Summary (InternalMesh) ---" << endl;
    cout << "1. Points Information" << endl;
    cout << "   - Total Points: " << mesh.numberOfPoints << endl;
    if (mesh.numberOfPoints > 0) {
        cout << "   - Points Coordinates (showing first 5):" << endl;
        for (long long i = 0; i < mesh.numberOfPoints && i < 5; ++i) {
            cout << "     - Point " << i << ": ("
                << mesh.points[i * 3 + 0] << ", "
                << mesh.points[i * 3 + 1] << ", "
                << mesh.points[i * 3 + 2] << ")"
                << endl;
        }
    }

    cout << "\n2. Cells Information" << endl;
    cout << "   - Total Cells: " << mesh.numberOfCells << endl;
    if (mesh.numberOfCells > 0) {
        cout << "   - Cell Connectivity (showing first 5):" << endl;
        for (long long i = 0; i < mesh.numberOfCells && i < 5; ++i) {
            const Cell& cell = mesh.cells[i]; // 使用引用避免拷贝
            cout << "     - Cell " << i << " (Type ID: " << cell.vtk_type << "): ";
            cout << "Connects " << cell.point_ids.size() << " points -> [ ";
            // *** 修改: 使用迭代器循环代替 range-based for ***
            for (std::vector<long long>::const_iterator it = cell.point_ids.begin(); it != cell.point_ids.end(); ++it) {
                cout << *it << " ";
            }
            cout << "]" << endl;
        }
    }

    cout << "\n3. Field Data Information" << endl;
    if (mesh.pointData.empty()) {
        cout << "   - Point Data: No fields read." << endl;
    }
    else {
        cout << "   - Point Data Fields Found: " << mesh.pointData.size() << endl;
        // *** 修改: 使用迭代器循环代替 range-based for ***
        for (std::map<std::string, std::vector<double>>::const_iterator it = mesh.pointData.begin(); it != mesh.pointData.end(); ++it) {
            cout << "     - Field '" << it->first << "' with " << it->second.size() << " values." << endl;
        }
    }
    if (mesh.cellData.empty()) {
        cout << "   - Cell Data: No fields read." << endl;
    }
    else {
        cout << "   - Cell Data Fields Found: " << mesh.cellData.size() << endl;
        // *** 修改: 使用迭代器循环代替 range-based for ***
        for (std::map<std::string, std::vector<double>>::const_iterator it = mesh.cellData.begin(); it != mesh.cellData.end(); ++it) {
            cout << "     - Field '" << it->first << "' with " << it->second.size() << " values." << endl;
        }
    }
    cout << "------------------------------------------" << endl;
}

vtkSmartPointer<vtkUnstructuredGrid> convertToVtkGrid(const InternalMesh& mesh) {
    // *** 修改: 将 auto 替换为明确的类型 ***
    vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // 1. 填充节点坐标
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(mesh.numberOfPoints);
    memcpy(points->GetVoidPointer(0), mesh.points.data(), mesh.numberOfPoints * 3 * sizeof(double));
    vtkGrid->SetPoints(points);

    // 2. 填充单元拓扑
    vtkGrid->Allocate(mesh.numberOfCells);
    // *** 修改: 使用迭代器循环代替 range-based for ***
    for (std::vector<Cell>::const_iterator it = mesh.cells.begin(); it != mesh.cells.end(); ++it) {
        const Cell& cell = *it;
        std::vector<vtkIdType> temp_ids;
        temp_ids.reserve(cell.point_ids.size());
        for (size_t i = 0; i < cell.point_ids.size(); ++i) {
            temp_ids.push_back(static_cast<vtkIdType>(cell.point_ids[i]));
        }
        vtkGrid->InsertNextCell(cell.vtk_type, temp_ids.size(), temp_ids.data());
    }

    // 3. 填充节点场数据
    // *** 修改: 使用迭代器循环代替 range-based for ***
    for (std::map<std::string, std::vector<double>>::const_iterator it = mesh.pointData.begin(); it != mesh.pointData.end(); ++it) {
        const std::string& name = it->first;
        const std::vector<double>& data = it->second;
        vtkSmartPointer<vtkDoubleArray> dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(name.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(data.size());
        memcpy(dataArray->GetVoidPointer(0), data.data(), data.size() * sizeof(double));
        vtkGrid->GetPointData()->AddArray(dataArray);
    }

    // 4. 填充单元场数据
    // *** 修改: 使用迭代器循环代替 range-based for ***
    for (std::map<std::string, std::vector<double>>::const_iterator it = mesh.cellData.begin(); it != mesh.cellData.end(); ++it) {
        const std::string& name = it->first;
        const std::vector<double>& data = it->second;
        vtkSmartPointer<vtkDoubleArray> dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(name.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(data.size());
        memcpy(dataArray->GetVoidPointer(0), data.data(), data.size() * sizeof(double));
        vtkGrid->GetCellData()->AddArray(dataArray);
    }

    return vtkGrid;
}

// ===============================================
// 4. VtkReader 类实现
// ===============================================

bool VtkReader::readFile(const std::string& filePath) {
    ifstream file(filePath.c_str()); // 使用 .c_str() 提高旧编译器兼容性
    if (!file.is_open()) { cerr << "Error opening file: " << filePath << endl; return false; }
    string line;
    while (getline(file, line)) {
        if (!parseLine(line, file)) { cerr << "Failed to parse line: " << line << endl; return false; }
    }
    file.close();
    return true;
}

std::string VtkReader::toUpper(std::string s) {
    // *** 修改: 使用简单的 for 循环代替 std::transform 和 lambda ***
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
        cerr << "Error during parsing keyword line: " << line << ". Details: " << e.what() << endl;
        return false;
    }
    return true;
}

// readPoints, readCells, readCellTypes, readPointData, readCellData 函数内部没有使用高级特性，无需修改
bool VtkReader::readPoints(std::ifstream& in, int numPoints) {
    m_meshData.numberOfPoints = numPoints;
    m_meshData.points.reserve(numPoints * 3);
    for (int i = 0; i < numPoints; ++i) {
        string line;
        if (!getline(in, line)) { cerr << "Error: File ended unexpectedly while reading points." << endl; return false; }
        stringstream ss(line);
        double x, y, z;
        ss >> x >> y >> z;
        if (ss.fail()) { cerr << "Error: Failed to parse point coordinates from line: " << line << endl; return false; }
        m_meshData.points.push_back(x);
        m_meshData.points.push_back(y);
        m_meshData.points.push_back(z);
    }
    return true;
}

bool VtkReader::readCells(std::ifstream& in, int numCells) {
    m_meshData.numberOfCells = numCells;
    m_meshData.cells.reserve(numCells);
    for (int i = 0; i < numCells; ++i) {
        string line;
        if (!getline(in, line)) { cerr << "Error: File ended unexpectedly while reading cells." << endl; return false; }
        stringstream ss(line);
        int num_points_in_cell;
        ss >> num_points_in_cell;
        if (ss.fail()) { cerr << "Error: Failed to parse cell structure from line: " << line << endl; return false; }
        Cell current_cell;
        current_cell.point_ids.reserve(num_points_in_cell);
        for (int j = 0; j < num_points_in_cell; ++j) {
            long long point_id;
            ss >> point_id;
            current_cell.point_ids.push_back(point_id);
        }
        m_meshData.cells.push_back(current_cell);
    }
    return true;
}

bool VtkReader::readCellTypes(std::ifstream& in, int numCellTypes) {
    if (m_meshData.cells.size() != (size_t)numCellTypes) { cerr << "Error: Mismatch between cell count and cell type count." << endl; return false; }
    for (int i = 0; i < numCellTypes; ++i) {
        string line;
        if (!getline(in, line)) { cerr << "Error: File ended unexpectedly while reading cell types." << endl; return false; }
        m_meshData.cells[i].vtk_type = std::stoi(line);
    }
    return true;
}

bool VtkReader::readPointData(std::ifstream& in, int numPoints) {
    string field_line;
    // 跳过可能的空行
    while (getline(in, field_line) && field_line.empty());

    vector<string> field_parts = Stringsplit(field_line);
    if (field_parts.size() < 3 || toUpper(field_parts[0]) != "FIELD") { cerr << "Error: Expected 'FIELD Data n' line, but found: " << field_line << endl; return false; }
    int num_fields = std::stoi(field_parts[2]);

    for (int i = 0; i < num_fields; ++i) {
        string var_def_line;
        if (!getline(in, var_def_line)) { cerr << "Error: File ended while expecting a variable definition." << endl; return false; }
        vector<string> var_parts = Stringsplit(var_def_line);
        if (var_parts.size() < 4) { cerr << "Error: Malformed variable definition line: " << var_def_line << endl; return false; }
        string var_name = var_parts[0];
        int num_components = std::stoi(var_parts[1]);
        int num_tuples = std::stoi(var_parts[2]);
        if (num_tuples != numPoints) { cerr << "Warning: Tuple count mismatch for point field '" << var_name << "'." << endl; }
        long long total_values = (long long)num_tuples * num_components;
        vector<double> values;
        values.reserve(total_values);
        for (long long j = 0; j < total_values; ++j) {
            double value;
            in >> value;
            if (in.fail()) { cerr << "Error reading value for point field '" << var_name << "'" << endl; return false; }
            values.push_back(value);
        }
        m_meshData.pointData[var_name] = values;
        cout << "--> Read point field '" << var_name << "' with " << values.size() << " values." << endl;
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return true;
}

bool VtkReader::readCellData(std::ifstream& in, int numCells) {
    string field_line;
    // 跳过可能的空行
    while (getline(in, field_line) && field_line.empty());

    vector<string> field_parts = Stringsplit(field_line);
    if (field_parts.size() < 3 || toUpper(field_parts[0]) != "FIELD") { cerr << "Error: Expected 'FIELD Data n' line, but found: " << field_line << endl; return false; }
    int num_fields = std::stoi(field_parts[2]);

    for (int i = 0; i < num_fields; ++i) {
        string var_def_line;
        if (!getline(in, var_def_line)) { cerr << "Error: File ended while expecting a variable definition." << endl; return false; }
        vector<string> var_parts = Stringsplit(var_def_line);
        if (var_parts.size() < 4) { cerr << "Error: Malformed variable definition line: " << var_def_line << endl; return false; }
        string var_name = var_parts[0];
        int num_components = std::stoi(var_parts[1]);
        int num_tuples = std::stoi(var_parts[2]);
        if (num_tuples != numCells) { cerr << "Warning: Tuple count mismatch for cell field '" << var_name << "'." << endl; }
        long long total_values = (long long)num_tuples * num_components;
        vector<double> values;
        values.reserve(total_values);
        for (long long j = 0; j < total_values; ++j) {
            double value;
            in >> value;
            if (in.fail()) { cerr << "Error reading value for cell field '" << var_name << "'" << endl; return false; }
            values.push_back(value);
        }
        m_meshData.cellData[var_name] = values;
        cout << "--> Read cell field '" << var_name << "' with " << values.size() << " values." << endl;
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return true;
}
// ===============================================
// 打印 vtkUnstructuredGrid 信息的函数 (无变化)
// ===============================================
void printGridInfo(vtkSmartPointer<vtkUnstructuredGrid> grid) {
    // ... (此函数内部未使用高级特性, 保持原样)
}
// ===============================================
// 动画回调与核心函数 (将 auto 替换为显式类型)
// ===============================================
class AnimationTimerCallback : public vtkCommand {
public:
    static AnimationTimerCallback* New() {
        return new AnimationTimerCallback;
    }

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

    // *** 修改: 将 auto 替换为明确的类型 ***
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0.1, 0.2, 0.4);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(1000, 800);
    renderWindow->SetWindowName("VTK Animation Player");
    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(grids[0]);
    mapper->SetScalarRange(globalMin, globalMax);

    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->AddRGBPoint(globalMin, 0.2, 0.2, 1.0); // 蓝色
    ctf->AddRGBPoint(0.0, 1.0, 1.0, 1.0);       // 白色
    ctf->AddRGBPoint(globalMax, 1.0, 0.2, 0.2); // 红色
    mapper->SetLookupTable(ctf);

    if (grids[0]->GetPointData()->GetScalars()) {
        mapper->SetScalarModeToUsePointData();
        mapper->SelectColorArray(grids[0]->GetPointData()->GetScalars()->GetName());
    }
    else if (grids[0]->GetCellData()->GetScalars()) {
        mapper->SetScalarModeToUseCellData();
        mapper->SelectColorArray(grids[0]->GetCellData()->GetScalars()->GetName());
    }

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    renderer->AddActor(actor);

    vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle(scalarBarTitle);
    scalarBar->SetNumberOfLabels(5);
    renderer->AddActor2D(scalarBar);

    vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
    textActor->SetInput(fileNames[0].c_str());
    textActor->GetTextProperty()->SetFontSize(16);
    textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    textActor->SetPosition(20, 20);
    renderer->AddActor2D(textActor);

    vtkSmartPointer<AnimationTimerCallback> callback = vtkSmartPointer<AnimationTimerCallback>::New();
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
    if (grids.empty()) return std::make_pair(0.0, 0.0);

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
        return std::make_pair(0.0, 0.0);
    }
    determinedArrayName = arrayName;
    std::cout << "\n--- Calculating Global Scalar Range from loaded data ('" << determinedArrayName << "') ---" << std::endl;

    for (size_t i = 0; i < grids.size(); ++i) {
        vtkDataSet* grid = grids[i]; // 使用指针，因为 vector 存储的是智能指针
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
    return std::make_pair(globalMin, globalMax);
}

// 提取文件名的辅助函数 (可选)
std::string getShortFileName(const std::string& filePath) {
    size_t last_slash_pos = filePath.find_last_of("/\\");
    if (std::string::npos != last_slash_pos) {
        return filePath.substr(last_slash_pos + 1);
    }
    return filePath;
}

int main() {
    // =================== START OF MODIFICATION ===================
    // 1. 设置文件列表 (手动管理，不再使用 <filesystem>)
    // *** 您需要在这里手动添加所有要加载的 .vtk 文件的完整路径 ***
    std::vector<std::string> vtkFilePaths;
    vtkFilePaths.push_back("D:/data/vtkout000000.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000152.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000304.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000454.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000604.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000754.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000904.vtk");
    vtkFilePaths.push_back("D:/data/vtkout001054.vtk");
    // ... 在这里添加更多文件路径

    if (vtkFilePaths.empty()) {
        std::cout << "No .vtk files specified in the vtkFilePaths vector." << std::endl;
        return 0;
    }

    std::sort(vtkFilePaths.begin(), vtkFilePaths.end());
    std::cout << "Processing " << vtkFilePaths.size() << " VTK files." << std::endl;

    // 2. 将所有文件读取并存储到 vector 中 (使用你自己的 VtkReader)
    std::vector<vtkSmartPointer<vtkDataSet>> allGrids;
    std::vector<std::string> shortFileNames;
    std::cout << "\n--- Loading all files into memory using YOUR custom VtkReader ---" << std::endl;

    // *** 修改: 使用传统的索引 for 循环代替 range-based for ***
    for (size_t i = 0; i < vtkFilePaths.size(); ++i) {
        const std::string& filePath = vtkFilePaths[i];
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
    // ===================  END OF MODIFICATION  ===================

    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Successfully loaded " << allGrids.size() << " datasets." << std::endl;

    // 3. 从已加载的数据中计算全局范围 (这部分代码无需改变)
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

    // 4. 调用动画函数 (这部分代码无需改变)
    animateGrids(allGrids, shortFileNames, globalRange.first, globalRange.second, scalarBarTitle.c_str());

    std::cout << "\nAnimation window closed. Exiting program." << std::endl;
    return 0;
}
