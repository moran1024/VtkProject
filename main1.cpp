#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
#include <limits>
#include <utility>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <mutex>
#include <functional>

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
#include <vtkCommand.h>
#include <vtkTextActor.h>

using namespace std;

class OptimizedVtkReader {
public:
    bool readFile(const std::string& filePath, vtkSmartPointer<vtkUnstructuredGrid> grid);

private:
    std::string toUpper(std::string s);
    std::vector<std::string> Stringsplit(const std::string& str);
    bool parseLine(const std::string& line, std::ifstream& file, vtkSmartPointer<vtkUnstructuredGrid> grid);
    bool readPoints(std::ifstream& in, int numPoints, vtkSmartPointer<vtkPoints> points);
    bool readCellsAndTypes(std::ifstream& file, int numCells, vtkSmartPointer<vtkUnstructuredGrid> grid);
    bool readPointData(std::ifstream& in, int numPoints, vtkSmartPointer<vtkUnstructuredGrid> grid);
    bool readCellData(std::ifstream& in, int numCells, vtkSmartPointer<vtkUnstructuredGrid> grid);
};

bool OptimizedVtkReader::readFile(const std::string& filePath, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    ifstream file(filePath);
    if (!file.is_open()) { cerr << "Error opening file: " << filePath << endl; return false; }
    string line;
    // Skip file header (assuming standard VTK legacy format)
    for (int i = 0; i < 5; ++i) {
        if (!getline(file, line)) return false;
    }
    while (getline(file, line)) {
        if (!parseLine(line, file, grid)) {
            return false;
        }
    }
    file.close();
    return true;
}

std::string OptimizedVtkReader::toUpper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
    return s;
}

std::vector<std::string> OptimizedVtkReader::Stringsplit(const std::string& str) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string token;
    while (ss >> token) { result.push_back(token); }
    return result;
}

bool OptimizedVtkReader::parseLine(const std::string& line, std::ifstream& file, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    if (line.empty() || line[0] == '#') return true;
    std::vector<std::string> parts = Stringsplit(line);
    if (parts.empty()) return true;

    std::string keyword = toUpper(parts[0]);
    try {
        if (keyword == "POINTS") {
            int numPoints = std::stoi(parts[1]);
            auto points = vtkSmartPointer<vtkPoints>::New();
            points->SetDataTypeToDouble();
            points->SetNumberOfPoints(numPoints);
            grid->SetPoints(points);
            return readPoints(file, numPoints, points);
        }
        else if (keyword == "CELLS") {
            int numCells = std::stoi(parts[1]);
            return readCellsAndTypes(file, numCells, grid);
        }
        else if (keyword == "POINT_DATA") {
            return readPointData(file, std::stoi(parts[1]), grid);
        }
        else if (keyword == "CELL_DATA") {
            return readCellData(file, std::stoi(parts[1]), grid);
        }
    }
    catch (const std::exception& e) {
        cerr << "Error during parsing keyword line: " << line << ". Details: " << e.what() << endl;
        return false;
    }
    return true;
}

bool OptimizedVtkReader::readPoints(std::ifstream& in, int numPoints, vtkSmartPointer<vtkPoints> points) {
    double* dataPtr = static_cast<double*>(points->GetVoidPointer(0));
    string line;
    for (int i = 0; i < numPoints; ++i) {
        if (!getline(in, line)) { cerr << "Error: File ended unexpectedly while reading points." << endl; return false; }

        const char* p = line.c_str();
        char* end;

        dataPtr[0] = strtod(p, &end);
        if (p == end) { cerr << "Error parsing point x on line: " << line << endl; return false; }
        p = end;

        dataPtr[1] = strtod(p, &end);
        if (p == end) { cerr << "Error parsing point y on line: " << line << endl; return false; }
        p = end;

        dataPtr[2] = strtod(p, &end);
        if (p == end) { cerr << "Error parsing point z on line: " << line << endl; return false; }

        dataPtr += 3;
    }
    return true;
}

bool OptimizedVtkReader::readCellsAndTypes(std::ifstream& file, int numCells, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    grid->Allocate(numCells);
    string line;

    std::vector<vtkIdType> cell_connectivity;
    // 预估一个合理的容量，比如平均每个单元有8个点
    cell_connectivity.reserve(numCells * 8);

    // 1. 读取 CELLS 数据块
    for (int i = 0; i < numCells; ++i) {
        if (!getline(file, line)) { cerr << "Error: File ended unexpectedly while reading cells." << endl; return false; }
        const char* p = line.c_str();
        char* end;

        long long num_points_in_cell = strtoll(p, &end, 10);
        if (p == end) { cerr << "Error parsing number of points in cell: " << line << endl; return false; }
        p = end;

        for (int j = 0; j < num_points_in_cell; ++j) {
            long long point_id = strtoll(p, &end, 10);
            if (p == end) { cerr << "Error parsing point ID in cell: " << line << endl; return false; }
            p = end;
            cell_connectivity.push_back(static_cast<vtkIdType>(point_id));
        }
    }

    // 2. 读取 CELL_TYPES 数据块
    // 首先读取并丢弃 "CELL_TYPES <num>" 这一行
    if (!getline(file, line)) { cerr << "Error: Missing CELL_TYPES line" << endl; return false; }

    std::vector<int> cell_types;
    cell_types.reserve(numCells);
    for (int i = 0; i < numCells; ++i) {
        int type;
        file >> type;
        if (file.fail()) { cerr << "Error reading cell type for cell " << i << endl; return false; }
        cell_types.push_back(type);
    }
    // 消费掉最后一行可能剩余的换行符
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // 3. 将单元插入到 VTK 网格中
    vtkIdType const* p_conn = cell_connectivity.data();
    for (int i = 0; i < numCells; ++i) {
        int cellType = cell_types[i];

        // 使用 switch 语句替换 vtkCellTypes::GetNumberOfPoints
        vtkIdType n_pts = 0;
        switch (cellType)
        {
        case VTK_TETRA: n_pts = 4; break;
        case VTK_HEXAHEDRON: n_pts = 8; break;
        case VTK_WEDGE: n_pts = 6; break;
        case VTK_PYRAMID: n_pts = 5; break;
        case VTK_QUAD: n_pts = 4; break;
        case VTK_TRIANGLE: n_pts = 3; break;
        case VTK_LINE: n_pts = 2; break;
        case VTK_VERTEX: n_pts = 1; break;
            // 根据需要添加其他单元类型
        default:
            cerr << "Warning: Unsupported or invalid cell type encountered: " << cellType << endl;
            break;
        }

        if (n_pts > 0) {
            grid->InsertNextCell(cellType, n_pts, p_conn);
            // 将连接性指针向前移动 n_pts 个位置
            p_conn += n_pts;
        }
    }
    return true;
}


bool OptimizedVtkReader::readPointData(std::ifstream& in, int numPoints, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    string field_line;
    while (getline(in, field_line) && field_line.empty());
    vector<string> field_parts = Stringsplit(field_line);
    if (field_parts.size() < 3 || toUpper(field_parts[0]) != "FIELD") { return false; } // Simplified error handling
    int num_fields = std::stoi(field_parts[2]);
    for (int i = 0; i < num_fields; ++i) {
        string var_def_line;
        if (!getline(in, var_def_line)) return false;
        vector<string> var_parts = Stringsplit(var_def_line);
        if (var_parts.size() < 4) return false;
        string var_name = var_parts[0];
        long long num_tuples = std::stoll(var_parts[2]);

        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(var_name.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(num_tuples);

        double* dataPtr = static_cast<double*>(dataArray->GetVoidPointer(0));
        for (long long j = 0; j < num_tuples; ++j) {
            in >> dataPtr[j];
            if (in.fail()) return false;
        }
        grid->GetPointData()->AddArray(dataArray);
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return true;
}

bool OptimizedVtkReader::readCellData(std::ifstream& in, int numCells, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    string field_line;
    while (getline(in, field_line) && field_line.empty());
    vector<string> field_parts = Stringsplit(field_line);
    if (field_parts.size() < 3 || toUpper(field_parts[0]) != "FIELD") { return false; }
    int num_fields = std::stoi(field_parts[2]);
    for (int i = 0; i < num_fields; ++i) {
        string var_def_line;
        if (!getline(in, var_def_line)) return false;
        vector<string> var_parts = Stringsplit(var_def_line);
        if (var_parts.size() < 4) return false;
        string var_name = var_parts[0];
        long long num_tuples = std::stoll(var_parts[2]);

        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(var_name.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(num_tuples);

        double* dataPtr = static_cast<double*>(dataArray->GetVoidPointer(0));
        for (long long j = 0; j < num_tuples; ++j) {
            in >> dataPtr[j];
            if (in.fail()) return false;
        }
        grid->GetCellData()->AddArray(dataArray);
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return true;
}


// ====================== Animation and Main Logic (No changes needed below) ======================

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

std::string getShortFileName(const std::string& filePath) {
    size_t last_slash_pos = filePath.find_last_of("/\\");
    if (std::string::npos != last_slash_pos) { return filePath.substr(last_slash_pos + 1); }
    return filePath;
}

std::mutex g_cout_mutex;
void print_safe(const std::string& message) {
    std::lock_guard<std::mutex> lock(g_cout_mutex);
    cout << message << endl;
}

struct FileReadResult {
    std::string filePath;
    vtkSmartPointer<vtkUnstructuredGrid> grid;
    bool operator<(const FileReadResult& other) const { return filePath < other.filePath; }
};

void load_files_worker(const std::vector<std::string>& files_to_process, std::vector<FileReadResult>& results) {
    for (const auto& filePath : files_to_process) {
        print_safe("Thread [" + std::to_string(std::hash<std::thread::id>()(std::this_thread::get_id())) + "] processing: " + getShortFileName(filePath));
        OptimizedVtkReader myReader;
        auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        if (myReader.readFile(filePath, grid)) {
            if (grid && grid->GetNumberOfPoints() > 0) {
                results.push_back({ filePath, grid });
            }
            else {
                print_safe("Warning: Grid is empty for " + filePath);
            }
        }
        else {
            print_safe("Warning: Failed to read file with custom reader: " + filePath);
        }
    }
}

int main() {
    std::vector<std::string> vtkFilePaths;
    vtkFilePaths.push_back("D:/data/vtkout000000.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000152.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000304.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000454.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000604.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000754.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000904.vtk");
    vtkFilePaths.push_back("D:/data/vtkout001054.vtk");
    // ... add other file paths ...

    if (vtkFilePaths.empty()) { cout << "No .vtk files specified." << endl; return 0; }
    std::sort(vtkFilePaths.begin(), vtkFilePaths.end());
    cout << "Found " << vtkFilePaths.size() << " VTK files to process." << endl << endl;

    // --- Single-Threaded Benchmark ---
    cout << "--- Starting Single-Threaded Benchmark ---" << endl;
    auto start_time_st = std::chrono::high_resolution_clock::now();
    std::vector<FileReadResult> single_thread_results;
    for (const auto& filePath : vtkFilePaths) {
        cout << "Reading: " << getShortFileName(filePath) << "..." << endl;
        OptimizedVtkReader myReader;
        auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        if (myReader.readFile(filePath, grid)) {
            if (grid && grid->GetNumberOfPoints() > 0) {
                single_thread_results.push_back({ filePath, grid });
            }
        }
    }
    auto end_time_st = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_st = end_time_st - start_time_st;
    cout << "--- Single-Threaded Benchmark Finished ---" << endl;
    cout << "Successfully loaded " << single_thread_results.size() << " datasets." << endl;
    cout << "Time taken: " << duration_st.count() << " seconds." << endl << endl;

    // --- Multi-Threaded Performance Test ---
    cout << "--- Starting Multi-Threaded Performance Test ---" << endl;
    auto start_time_mt = std::chrono::high_resolution_clock::now();
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 2;
    num_threads = std::min(num_threads, (unsigned int)vtkFilePaths.size());
    print_safe("Using " + std::to_string(num_threads) + " threads.");
    std::vector<std::thread> threads;
    std::vector<std::vector<FileReadResult>> results_per_thread(num_threads);
    size_t start_index = 0;
    for (unsigned int i = 0; i < num_threads; ++i) {
        size_t files_for_this_thread = vtkFilePaths.size() / num_threads + (i < vtkFilePaths.size() % num_threads ? 1 : 0);
        size_t end_index = start_index + files_for_this_thread;
        std::vector<std::string> files_subset(vtkFilePaths.begin() + start_index, vtkFilePaths.begin() + end_index);
        threads.push_back(std::thread(load_files_worker, files_subset, std::ref(results_per_thread[i])));
        start_index = end_index;
    }
    for (auto& t : threads) { t.join(); }
    std::vector<FileReadResult> multi_thread_results;
    for (const auto& thread_results : results_per_thread) {
        multi_thread_results.insert(multi_thread_results.end(), thread_results.begin(), thread_results.end());
    }
    auto end_time_mt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_mt = end_time_mt - start_time_mt;
    cout << "--- Multi-Threaded Test Finished ---" << endl;
    cout << "Successfully loaded " << multi_thread_results.size() << " datasets." << endl;
    cout << "Time taken: " << duration_mt.count() << " seconds." << endl << endl;

    // --- Summary and Visualization ---
    cout << "================= Performance Summary =================" << endl;
    cout << "Single-Threaded Time: " << duration_st.count() << " s" << endl;
    cout << "Multi-Threaded Time:  " << duration_mt.count() << " s" << endl;
    if (duration_mt.count() > 1e-9) {
        double speedup = duration_st.count() / duration_mt.count();
        cout.precision(2);
        cout << "Speedup: " << std::fixed << speedup << "x" << endl;
    }
    cout << "=====================================================" << endl << endl;

    std::sort(multi_thread_results.begin(), multi_thread_results.end());
    std::vector<vtkSmartPointer<vtkDataSet>> allGrids;
    std::vector<std::string> shortFileNames;
    for (const auto& res : multi_thread_results) {
        allGrids.push_back(res.grid);
        shortFileNames.push_back(getShortFileName(res.filePath));
    }

    if (!allGrids.empty()) {
        std::string scalarBarTitle;
        std::pair<double, double> globalRange = calculateGlobalScalarRange(allGrids, scalarBarTitle);
        animateGrids(allGrids, shortFileNames, globalRange.first, globalRange.second, scalarBarTitle.c_str());
    }
    else {
        cout << "No datasets were loaded, cannot start visualization." << endl;
    }

    cout << "\nAnimation window closed. Exiting program." << endl;
    return 0;
}
