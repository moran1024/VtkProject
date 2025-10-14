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
#include <future>

using namespace std;

// =================================================================================
// =================== 高性能重构: OptimizedVtkReader 类 ======================
// =================================================================================
// 策略: 将整个文件一次性读入内存缓冲区，然后使用指针进行极速解析。
// 这将 I/O 开销降至最低，并利用了 C-style 函数 (strtod, strtoll, strstr) 的高性能。

class OptimizedVtkReader {
public:
    // 公共 API 保持不变
    bool readFile(const std::string& filePath, vtkSmartPointer<vtkUnstructuredGrid> grid);

private:
    // 辅助函数: 在内存缓冲区中导航
    void skipLine(const char*& p);
    void skipWhitespace(const char*& p);
    std::string getNextLine(const char*& p);
    std::string toUpper(std::string s); // 复用

    // 核心解析逻辑，现在操作内存指针而不是 ifstream
    bool parseBuffer(const char* buffer, vtkSmartPointer<vtkUnstructuredGrid> grid);

    // 解析函数已修改为接受 const char*&
    bool readPoints(const char*& cursor, int numPoints, vtkSmartPointer<vtkPoints> points);
    bool readCellsAndTypes(const char*& cursor, int numCells, int totalCellSize, vtkSmartPointer<vtkUnstructuredGrid> grid);
    bool readPointData(const char*& cursor, int numPoints, vtkSmartPointer<vtkUnstructuredGrid> grid);
    bool readCellData(const char*& cursor, int numCells, vtkSmartPointer<vtkUnstructuredGrid> grid);
};

// --- 实现 ---

// 主入口函数: 将整个文件读入内存
bool OptimizedVtkReader::readFile(const std::string& filePath, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    ifstream file(filePath, std::ios::binary);
    if (!file) {
        cerr << "Error opening file: " << filePath << endl;
        return false;
    }

    // 获取文件大小并将整个文件读入字符串缓冲区
    file.seekg(0, std::ios::end);
    std::string buffer(file.tellg(), '\0');
    file.seekg(0);
    file.read(&buffer[0], buffer.size());

    if (buffer.empty()) {
        cerr << "Error: File is empty or could not be read: " << filePath << endl;
        return false;
    }

    // 从内存缓冲区开始解析
    return parseBuffer(buffer.c_str(), grid);
}

// 核心调度程序: 在内存缓冲区中查找关键字并调用相应的解析器
bool OptimizedVtkReader::parseBuffer(const char* buffer, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    const char* cursor = buffer;

    // 1. 查找 DATASET 关键字
    cursor = strstr(cursor, "DATASET");
    if (!cursor) {
        cerr << "Error: Could not find 'DATASET' keyword in file." << endl;
        return false;
    }
    skipLine(cursor); // 移动到下一行

    try {
        // 循环直到缓冲区末尾
        while (cursor && *cursor != '\0') {
            skipWhitespace(cursor);
            if (*cursor == '\0') break;

            const char* line_start = cursor;
            skipLine(cursor);
            std::string line(line_start, cursor - line_start);

            std::stringstream ss(line);
            std::string keyword;
            ss >> keyword;
            keyword = toUpper(keyword);

            if (keyword == "POINTS") {
                int numPoints;
                ss >> numPoints;
                auto points = vtkSmartPointer<vtkPoints>::New();
                points->SetDataTypeToDouble();
                points->SetNumberOfPoints(numPoints);
                if (!readPoints(cursor, numPoints, points)) return false;
                grid->SetPoints(points);
            }
            else if (keyword == "CELLS") {
                int numCells, totalCellSize;
                ss >> numCells >> totalCellSize;
                if (!readCellsAndTypes(cursor, numCells, totalCellSize, grid)) return false;
            }
            else if (keyword == "POINT_DATA") {
                int numPoints;
                ss >> numPoints;
                if (!readPointData(cursor, numPoints, grid)) return false;
            }
            else if (keyword == "CELL_DATA") {
                int numCells;
                ss >> numCells;
                if (!readCellData(cursor, numCells, grid)) return false;
            }
        }
    }
    catch (const std::exception& e) {
        cerr << "Error during buffer parsing. Details: " << e.what() << endl;
        return false;
    }
    return true;
}

// 从内存缓冲区极速读取点
bool OptimizedVtkReader::readPoints(const char*& cursor, int numPoints, vtkSmartPointer<vtkPoints> points) {
    double* dataPtr = static_cast<double*>(points->GetVoidPointer(0));
    char* end;

    for (int i = 0; i < numPoints; ++i) {
        skipWhitespace(cursor);
        if (*cursor == '\0') { cerr << "Error: Buffer ended unexpectedly while reading points." << endl; return false; }

        dataPtr[0] = strtod(cursor, &end); if (cursor == end) return false; cursor = end;
        dataPtr[1] = strtod(cursor, &end); if (cursor == end) return false; cursor = end;
        dataPtr[2] = strtod(cursor, &end); if (cursor == end) return false; cursor = end;
        dataPtr += 3;
    }
    return true;
}

// 从内存缓冲区极速读取单元和类型
bool OptimizedVtkReader::readCellsAndTypes(const char*& cursor, int numCells, int totalCellSize, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
    connectivity->SetNumberOfValues(totalCellSize);
    vtkIdType* connPtr = connectivity->GetPointer(0);
    char* end;

    // 1. 读取单元连接性
    for (int i = 0; i < numCells; ++i) {
        skipWhitespace(cursor);
        if (*cursor == '\0') { cerr << "Error reading cell connectivity for cell " << i << endl; return false; }

        long long num_points_in_cell = strtoll(cursor, &end, 10);
        if (cursor == end) return false; cursor = end;
        *connPtr++ = static_cast<vtkIdType>(num_points_in_cell);

        for (int j = 0; j < num_points_in_cell; ++j) {
            long long point_id = strtoll(cursor, &end, 10);
            if (cursor == end) return false; cursor = end;
            *connPtr++ = static_cast<vtkIdType>(point_id);
        }
    }

    // 2. 查找 CELL_TYPES 关键字并读取类型
    cursor = strstr(cursor, "CELL_TYPES");
    if (!cursor) { cerr << "Error: Missing CELL_TYPES line" << endl; return false; }
    skipLine(cursor); // 跳过 "CELL_TYPES <num>" 这一行

    std::vector<int> cell_types;
    cell_types.reserve(numCells);
    for (int i = 0; i < numCells; ++i) {
        skipWhitespace(cursor);
        if (*cursor == '\0') { cerr << "Error reading cell type for cell " << i << endl; return false; }
        long long type = strtoll(cursor, &end, 10);
        if (cursor == end) return false; cursor = end;
        cell_types.push_back(static_cast<int>(type));
    }

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetCells(numCells, connectivity);
    grid->SetCells(cell_types.data(), cellArray);

    return true;
}

// 从内存缓冲区极速读取点数据
bool OptimizedVtkReader::readPointData(const char*& cursor, int numPoints, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    skipWhitespace(cursor);
    const char* line_start = cursor;
    skipLine(cursor);
    std::string line(line_start, cursor - line_start);
    std::stringstream ss(line);

    string keyword;
    ss >> keyword;
    keyword = toUpper(keyword);

    if (keyword == "SCALARS") {
        string var_name;
        ss >> var_name;
        skipLine(cursor); // 跳过 LOOKUP_TABLE 行

        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(var_name.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(numPoints);

        double* dataPtr = dataArray->GetPointer(0);
        char* end;
        for (int i = 0; i < numPoints; ++i) {
            skipWhitespace(cursor);
            if (*cursor == '\0') return false;
            dataPtr[i] = strtod(cursor, &end);
            if (cursor == end) return false; cursor = end;
        }
        grid->GetPointData()->AddArray(dataArray);
        grid->GetPointData()->SetActiveScalars(var_name.c_str());

    }
    else if (keyword == "FIELD") {
        string fieldName;
        int num_fields;
        ss >> fieldName >> num_fields;

        for (int i = 0; i < num_fields; ++i) {
            std::string var_def_line = getNextLine(cursor);
            std::stringstream var_ss(var_def_line);
            string var_name;
            int num_comp, num_tuples;
            string data_type;
            var_ss >> var_name >> num_comp >> num_tuples >> data_type;

            auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
            dataArray->SetName(var_name.c_str());
            dataArray->SetNumberOfComponents(num_comp);
            dataArray->SetNumberOfTuples(num_tuples);

            double* dataPtr = dataArray->GetPointer(0);
            char* end;
            for (int j = 0; j < num_tuples * num_comp; ++j) {
                skipWhitespace(cursor);
                if (*cursor == '\0') return false;
                dataPtr[j] = strtod(cursor, &end);
                if (cursor == end) return false; cursor = end;
            }
            grid->GetPointData()->AddArray(dataArray);
        }
    }
    else {
        cerr << "Unsupported POINT_DATA type: " << keyword << endl;
        return false;
    }
    return true;
}

// 从内存缓冲区极速读取单元数据
bool OptimizedVtkReader::readCellData(const char*& cursor, int numCells, vtkSmartPointer<vtkUnstructuredGrid> grid) {
    skipWhitespace(cursor);
    const char* line_start = cursor;
    skipLine(cursor);
    std::string line(line_start, cursor - line_start);
    std::stringstream ss(line);

    string keyword;
    ss >> keyword;
    keyword = toUpper(keyword);

    if (keyword == "SCALARS") {
        string var_name;
        ss >> var_name;
        skipLine(cursor); // 跳过 LOOKUP_TABLE 行

        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(var_name.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(numCells);

        double* dataPtr = dataArray->GetPointer(0);
        char* end;
        for (int i = 0; i < numCells; ++i) {
            skipWhitespace(cursor);
            if (*cursor == '\0') return false;
            dataPtr[i] = strtod(cursor, &end);
            if (cursor == end) return false; cursor = end;
        }
        grid->GetCellData()->AddArray(dataArray);
        grid->GetCellData()->SetActiveScalars(var_name.c_str());
    }
    // ================== [FIXED] 使用正确的内存指针解析逻辑 ==================
    else if (keyword == "FIELD") {
        string fieldName; // VTK 文件格式的一部分，但此处未使用
        int num_fields;
        ss >> fieldName >> num_fields;

        for (int i = 0; i < num_fields; ++i) {
            std::string var_def_line = getNextLine(cursor);
            std::stringstream var_ss(var_def_line);
            string var_name;
            int num_comp, num_tuples;
            string data_type;
            var_ss >> var_name >> num_comp >> num_tuples >> data_type;

            auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
            dataArray->SetName(var_name.c_str());
            dataArray->SetNumberOfComponents(num_comp);
            dataArray->SetNumberOfTuples(num_tuples);

            double* dataPtr = dataArray->GetPointer(0);
            char* end;
            for (int j = 0; j < num_tuples * num_comp; ++j) {
                skipWhitespace(cursor);
                if (*cursor == '\0') {
                    cerr << "Error: Buffer ended unexpectedly while reading CELL_DATA FIELD " << var_name << endl;
                    return false;
                }
                dataPtr[j] = strtod(cursor, &end);
                if (cursor == end) {
                    cerr << "Error: Failed to parse number in CELL_DATA FIELD " << var_name << endl;
                    return false;
                }
                cursor = end;
            }
            // 关键区别：将数组添加到 CellData
            grid->GetCellData()->AddArray(dataArray);
        }
    }
    // ======================== [FIX END] ===========================
    else {
        cerr << "Unsupported CELL_DATA type: " << keyword << endl;
        return false;
    }
    return true;
}

// --- 辅助函数实现 ---
void OptimizedVtkReader::skipLine(const char*& p) {
    while (*p != '\n' && *p != '\0') {
        p++;
    }
    if (*p == '\n') {
        p++;
    }
}

void OptimizedVtkReader::skipWhitespace(const char*& p) {
    while (*p != '\0' && isspace(static_cast<unsigned char>(*p))) {
        p++;
    }
}

std::string OptimizedVtkReader::getNextLine(const char*& p) {
    skipWhitespace(p);
    const char* start = p;
    skipLine(p);
    const char* end = p;
    // 去除行尾的空白符
    if (end > start) {
        end--; // 指向最后一个字符
        while (end > start && isspace(static_cast<unsigned char>(*end))) {
            end--;
        }
        end++; // 指向最后一个非空白字符的后面
    }
    return std::string(start, end - start);
}

std::string OptimizedVtkReader::toUpper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
    return s;
}



// AnimationTimerCallback 类
class AnimationTimerCallback : public vtkCommand {
public:
    static AnimationTimerCallback* New() {
        return new AnimationTimerCallback;
    }

    // 【核心修改】在这里设置活动标量
    void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
        this->CurrentFrameIndex = (this->CurrentFrameIndex + 1) % this->Grids->size();
        vtkDataSet* currentGrid = (*this->Grids)[this->CurrentFrameIndex];

        // --- 新增代码 ---
        // 在将数据送入Mapper之前，确保当前Grid的活动标量是正确的
        if (currentGrid->GetPointData()->GetArray(this->ScalarArrayName.c_str())) {
            currentGrid->GetPointData()->SetActiveScalars(this->ScalarArrayName.c_str());
        }
        else if (currentGrid->GetCellData()->GetArray(this->ScalarArrayName.c_str())) {
            currentGrid->GetCellData()->SetActiveScalars(this->ScalarArrayName.c_str());
        }
        // --- 新增代码结束 ---

        this->Mapper->SetInputData(currentGrid);
        this->TextActor->SetInput((*this->FileNames)[this->CurrentFrameIndex].c_str());
        this->RenderWindow->Render();
    }

    // 【修改】更新 SetData 函数以接收数组名称
    void SetData(
        std::vector<vtkSmartPointer<vtkDataSet>>* grids,
        const std::vector<std::string>* fileNames,
        vtkDataSetMapper* mapper,
        vtkTextActor* textActor,
        vtkRenderWindow* renderWindow,
        const std::string& scalarArrayName) { // <-- 新增参数
        this->Grids = grids;
        this->FileNames = fileNames;
        this->Mapper = mapper;
        this->TextActor = textActor;
        this->RenderWindow = renderWindow;
        this->ScalarArrayName = scalarArrayName; // <-- 保存数组名称
    }

private:
    int CurrentFrameIndex = 0;
    std::vector<vtkSmartPointer<vtkDataSet>>* Grids = nullptr;
    const std::vector<std::string>* FileNames = nullptr;
    vtkDataSetMapper* Mapper = nullptr;
    vtkTextActor* TextActor = nullptr;
    vtkRenderWindow* RenderWindow = nullptr;
    std::string ScalarArrayName; // <-- 新增成员变量
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
    callback->SetData(&grids, &fileNames, mapper, textActor, renderWindow, std::string(scalarBarTitle));

    interactor->AddObserver(vtkCommand::TimerEvent, callback);
    interactor->Initialize();
    interactor->CreateRepeatingTimer(100);

    std::cout << "\nStarting animation... Close the window to stop." << std::endl;
    renderWindow->Render();
    interactor->Start();
}

std::pair<double, double> calculateGlobalScalarRange(
    const std::vector<vtkSmartPointer<vtkDataSet>>& grids,
    std::string& determinedArrayName)
{
    if (grids.empty()) {
        return std::make_pair(0.0, 0.0);
    }

    // --- 1. 首先，在主线程中确定要使用的数组名称 ---
    const char* arrayName = nullptr;
    if (grids[0]->GetPointData()->GetNumberOfArrays() > 0) {
        arrayName = grids[0]->GetPointData()->GetArray(0)->GetName();
    }
    else if (grids[0]->GetCellData()->GetNumberOfArrays() > 0) {
        arrayName = grids[0]->GetCellData()->GetArray(0)->GetName();
    }
    else {
        
        return std::make_pair(0.0, 0.0);
    }
    determinedArrayName = arrayName;
    std::cout << "\n--- Calculating Global Scalar Range from loaded data ('" << determinedArrayName << "') ---" << std::endl;

    // --- 2. 定义并行任务 (Worker Lambda) ---
    //   [&] 捕获外部变量(grids, arrayName)的引用
    auto worker = [&](size_t start, size_t end) -> std::pair<double, double> {
        double localMin = std::numeric_limits<double>::max();
        double localMax = std::numeric_limits<double>::lowest();

        for (size_t i = start; i < end; ++i) {
            vtkDataSet* grid = grids[i];

            // 【关键修正】: 在循环内部为当前 grid 获取 dataArray
            vtkDataArray* dataArray = grid->GetPointData()->GetArray(arrayName);
            if (!dataArray) {
                dataArray = grid->GetCellData()->GetArray(arrayName);
            }

            // 现在可以安全地使用 dataArray
            if (dataArray) {
                double range[2];
                dataArray->GetRange(range);
                localMin = std::min(localMin, range[0]);
                localMax = std::max(localMax, range[1]);
            }
        }
        return { localMin, localMax };
        };

    // --- 3. 分发任务到线程池 ---
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (grids.size() < num_threads) {
        num_threads = grids.size();
    }
    std::vector<std::future<std::pair<double, double>>> futures;
    size_t start_index = 0;

    for (unsigned int i = 0; i < num_threads; ++i) {
        size_t files_for_this_thread = grids.size() / num_threads + (i < grids.size() % num_threads ? 1 : 0);
        if (files_for_this_thread == 0) continue;
        size_t end_index = start_index + files_for_this_thread;
        futures.push_back(std::async(std::launch::async, worker, start_index, end_index));
        start_index = end_index;
    }

    // --- 4. 收集并合并结果 (Reduction) ---
    double globalMin = std::numeric_limits<double>::max();   // 【改进】: 只在这里定义一次
    double globalMax = std::numeric_limits<double>::lowest(); // 【改进】: 只在这里定义一次

    for (auto& f : futures) {
        auto local_range = f.get(); // f.get()会等待线程完成并返回结果
        globalMin = std::min(globalMin, local_range.first);
        globalMax = std::max(globalMax, local_range.second);
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    return { globalMin, globalMax };
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
    cout << "--- Starting Single-Threaded Benchmark (with new high-performance reader) ---" << endl;
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
