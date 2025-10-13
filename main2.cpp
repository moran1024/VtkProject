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
#include <vtkIdTypeArray.h>
#include <vtkUnsignedCharArray.h>

using namespace std;

// =================================================================================
// =================== 重构后的 OptimizedVtkReader 类 =========================
// =================================================================================

class OptimizedVtkReader {
public:
    bool readFile(const std::string& filePath, vtkSmartPointer<vtkUnstructuredGrid> grid);

private:
    std::string toUpper(std::string s);
    std::vector<std::string> Stringsplit(const std::string& str);

    bool readDataset(std::ifstream& file, vtkSmartPointer<vtkUnstructuredGrid> grid);

    bool readPoints(std::ifstream& in, int numPoints, vtkSmartPointer<vtkPoints> points);

    // MODIFIED: 专门为 VTK 8.2.0 修改了此函数
    bool readCellsAndTypes(std::ifstream& file, int numCells, int totalCellSize, vtkSmartPointer<vtkUnstructuredGrid> grid);

    bool readPointData(std::ifstream& in, int numPoints, vtkSmartPointer<vtkUnstructuredGrid> grid);
    bool readCellData(std::ifstream& in, int numCells, vtkSmartPointer<vtkUnstructuredGrid> grid);
};

// 更稳健的头文件解析
bool OptimizedVtkReader::readFile(const std::string& filePath, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    ifstream file(filePath);
    if (!file.is_open()) { cerr << "Error opening file: " << filePath << endl; return false; }

    string line;
    // 跳过文件头和注释，直到找到 DATASET 关键字行
    while (getline(file, line)) {
        std::string upper_line = toUpper(line);
        if (upper_line.rfind("DATASET", 0) == 0) {
            // 找到了数据集的开头，现在读取文件的剩余部分
            return readDataset(file, grid);
        }
    }

    cerr << "Error: Could not find 'DATASET' keyword in file: " << filePath << endl;
    return false;
}

// 处理文件头之后的主要解析逻辑
bool OptimizedVtkReader::readDataset(std::ifstream& file, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    string line;
    try {
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::vector<std::string> parts = Stringsplit(line);
            if (parts.empty()) continue;

            std::string keyword = toUpper(parts[0]);

            if (keyword == "POINTS") {
                int numPoints = std::stoi(parts[1]);
                auto points = vtkSmartPointer<vtkPoints>::New();
                points->SetDataTypeToDouble();
                points->SetNumberOfPoints(numPoints);
                if (!readPoints(file, numPoints, points)) return false;
                grid->SetPoints(points);
            }
            else if (keyword == "CELLS") {
                int numCells = std::stoi(parts[1]);
                int totalCellSize = std::stoi(parts[2]); // 单元格列表中的整数总数
                if (!readCellsAndTypes(file, numCells, totalCellSize, grid)) return false;
            }
            else if (keyword == "POINT_DATA") {
                if (!readPointData(file, std::stoi(parts[1]), grid)) return false;
            }
            else if (keyword == "CELL_DATA") {
                if (!readCellData(file, std::stoi(parts[1]), grid)) return false;
            }
        }
    }
    catch (const std::exception& e) {
        cerr << "Error during file parsing. Details: " << e.what() << endl;
        return false;
    }
    return true;
}

// 此函数已经经过了很好的优化，无需修改
bool OptimizedVtkReader::readPoints(std::ifstream& in, int numPoints, vtkSmartPointer<vtkPoints> points) {
    double* dataPtr = static_cast<double*>(points->GetVoidPointer(0));
    string line;
    for (int i = 0; i < numPoints; ++i) {
        if (!getline(in, line)) { cerr << "Error: File ended unexpectedly while reading points." << endl; return false; }
        const char* p = line.c_str();
        char* end;
        dataPtr[0] = strtod(p, &end); if (p == end) { cerr << "Error parsing point x on line: " << line << endl; return false; } p = end;
        dataPtr[1] = strtod(p, &end); if (p == end) { cerr << "Error parsing point y on line: " << line << endl; return false; } p = end;
        dataPtr[2] = strtod(p, &end); if (p == end) { cerr << "Error parsing point z on line: " << line << endl; return false; }
        dataPtr += 3;
    }
    return true;
}

// ############# 已再次修改此函数以兼容 VTK 8.2.0 #############
bool OptimizedVtkReader::readCellsAndTypes(std::ifstream& file, int numCells, int totalCellSize, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    // 1. 将所有单元格连接性数据读入一个 vtkIdTypeArray 缓冲区
    auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
    connectivity->SetNumberOfValues(totalCellSize);
    vtkIdType* connPtr = connectivity->GetPointer(0);

    string line;
    for (int i = 0; i < numCells; ++i) {
        if (!getline(file, line)) { cerr << "Error reading cell connectivity for cell " << i << endl; return false; }

        const char* p = line.c_str();
        char* end;
        long long num_points_in_cell = strtoll(p, &end, 10);
        if (p == end) { cerr << "Error parsing number of points in cell: " << line << endl; return false; }
        p = end;

        *connPtr++ = static_cast<vtkIdType>(num_points_in_cell);

        for (int j = 0; j < num_points_in_cell; ++j) {
            long long point_id = strtoll(p, &end, 10);
            if (p == end) { cerr << "Error parsing point ID in cell: " << line << endl; return false; }
            p = end;
            *connPtr++ = static_cast<vtkIdType>(point_id);
        }
    }

    // 2. 将所有单元格类型读入一个 std::vector<int> 缓冲区
    // 跳过 "CELL_TYPES <num>" 这一行
    if (!getline(file, line)) { cerr << "Error: Missing CELL_TYPES line" << endl; return false; }

    std::vector<int> cell_types;
    cell_types.reserve(numCells);
    for (int i = 0; i < numCells; ++i) {
        int type;
        file >> type;
        if (file.fail()) { cerr << "Error reading cell type for cell " << i << endl; return false; }
        cell_types.push_back(type);
    }
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // 3. 创建一个 vtkCellArray 并使用 SetCells 进行批量导入
    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetCells(numCells, connectivity);

    // 4. 在非结构化网格上设置单元格和类型 (使用 int* 指针)
    grid->SetCells(cell_types.data(), cellArray);

    return true;
}

// 添加了对常见 "SCALARS" 关键字的支持
bool OptimizedVtkReader::readPointData(std::ifstream& in, int numPoints, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    string data_header_line;
    while (getline(in, data_header_line) && data_header_line.empty());

    vector<string> parts = Stringsplit(data_header_line);
    if (parts.empty()) return false;

    string keyword = toUpper(parts[0]);

    if (keyword == "SCALARS") {
        if (parts.size() < 2) return false;
        string var_name = parts[1];

        // 跳过 LOOKUP_TABLE 这一行
        getline(in, data_header_line);

        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(var_name.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(numPoints);

        double* dataPtr = dataArray->GetPointer(0);
        for (int i = 0; i < numPoints; ++i) {
            in >> dataPtr[i];
            if (in.fail()) return false;
        }
        grid->GetPointData()->AddArray(dataArray);
        grid->GetPointData()->SetActiveScalars(var_name.c_str());
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    }
    else if (keyword == "FIELD") {
        // 原来的 FIELD 逻辑保持不变
        int num_fields = std::stoi(parts[2]);
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
    }
    else {
        cerr << "Unsupported POINT_DATA type: " << keyword << endl;
        return false;
    }
    return true;
}

// 添加了对常见 "SCALARS" 关键字的支持
bool OptimizedVtkReader::readCellData(std::ifstream& in, int numCells, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    string data_header_line;
    while (getline(in, data_header_line) && data_header_line.empty());

    vector<string> parts = Stringsplit(data_header_line);
    if (parts.empty()) return false;

    string keyword = toUpper(parts[0]);

    if (keyword == "SCALARS") {
        if (parts.size() < 2) return false;
        string var_name = parts[1];

        getline(in, data_header_line); // 跳过 LOOKUP_TABLE

        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(var_name.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(numCells);

        double* dataPtr = dataArray->GetPointer(0);
        for (int i = 0; i < numCells; ++i) {
            in >> dataPtr[i];
            if (in.fail()) return false;
        }
        grid->GetCellData()->AddArray(dataArray);
        grid->GetCellData()->SetActiveScalars(var_name.c_str());
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    }
    else if (keyword == "FIELD") {
        // 原有逻辑
        int num_fields = std::stoi(parts[2]);
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
    }
    else {
        cerr << "Unsupported CELL_DATA type: " << keyword << endl;
        return false;
    }
    return true;
}

// 未修改的辅助函数
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


// ====================== 动画和主程序逻辑 (以下无需更改) ======================

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
    ctf->AddRGBPoint((globalMin + globalMax) / 2.0, 1.0, 1.0, 1.0); // 白色
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
        vtkDataSet* grid = grids[i];
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
    // 请确保将路径修改为您的文件实际存放的位置
    vtkFilePaths.push_back("D:/data/vtkout000000.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000152.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000304.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000454.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000604.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000754.vtk");
    vtkFilePaths.push_back("D:/data/vtkout000904.vtk");
    vtkFilePaths.push_back("D:/data/vtkout001054.vtk");
    // ... 在此添加其他文件路径 ...

    if (vtkFilePaths.empty()) { cout << "No .vtk files specified." << endl; return 0; }
    std::sort(vtkFilePaths.begin(), vtkFilePaths.end());
    cout << "Found " << vtkFilePaths.size() << " VTK files to process." << endl << endl;

    // --- 单线程基准测试 ---
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

    // --- 多线程性能测试 ---
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

    // --- 总结与可视化 ---
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

