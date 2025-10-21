#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <algorithm> // 为了使用 std::min
#include "H5Cpp.h"

// 为了方便使用，声明命名空间
using namespace std;
using namespace H5;

// --- 数据结构 (未改变) ---
struct Triplet { float x, y, z; };
struct AttributeDetails { string value; string typeName; vector<hsize_t> dims; string length; string padding; string charset; };
struct ElementClassInfo { AttributeDetails elementType; AttributeDetails sectionCategory; bool isValid = false; };
struct StepClassInfo { AttributeDetails Description; AttributeDetails Domain; AttributeDetails Index; bool isValid = false; };
struct FrameClassInfo { AttributeDetails Description; AttributeDetails Inc_Mode; AttributeDetails LoadCase; AttributeDetails Time_Freq; bool isValid = false; };
struct FLUIDClassInfo { AttributeDetails ComponentLabels; AttributeDetails Invariants; AttributeDetails Position; AttributeDetails Type; bool isValid = false; };

// --- 辅助函数和读取函数 (未改变) ---

AttributeDetails readAttribute(const H5Object& obj, const string& attrName) {
    AttributeDetails details;
    if (!obj.attrExists(attrName)) { return details; }
    try {
        Attribute attr = obj.openAttribute(attrName);
        DataType type = attr.getDataType();
        DataSpace space = attr.getSpace();
        H5T_class_t type_class = type.getClass();
        switch (type_class) {
        case H5T_INTEGER: details.typeName = "Integer"; break;
        case H5T_FLOAT:   details.typeName = "Float";   break;
        case H5T_STRING:  details.typeName = "String";  break;
        default:          details.typeName = "Unknown"; break;
        }
        int rank = space.getSimpleExtentNdims();
        if (rank > 0) {
            details.dims.resize(rank);
            space.getSimpleExtentDims(details.dims.data());
        }
        else {
            details.dims.assign(1, 1);
        }
        if (type_class == H5T_STRING) {
            attr.read(type, details.value);
            StrType str_type = attr.getStrType();
            details.length = str_type.isVariableStr() ? "variable" : to_string(str_type.getSize());
            details.padding = (str_type.getStrpad() == H5T_STR_NULLTERM) ? "NULLTERM" : "SPACEPAD";
            details.charset = (str_type.getCset() == H5T_CSET_ASCII) ? "ASCII" : "UTF8";
        }
        else if (type_class == H5T_INTEGER) {
            long long temp_val = 0;
            attr.read(PredType::NATIVE_LLONG, &temp_val);
            details.value = to_string(temp_val);
        }
        else if (type_class == H5T_FLOAT) {
            double temp_val = 0;
            attr.read(PredType::NATIVE_DOUBLE, &temp_val);
            stringstream ss;
            ss << scientific << temp_val;
            details.value = ss.str();
        }
    }
    catch (const H5::Exception& error) {
        cerr << "Error reading attribute: " << attrName << endl;
        error.printErrorStack();
        return {};
    }
    return details;
}

template <typename T>
vector<T> ReadFlattenedDataset(const H5File& file, const string& dataset_path, const PredType& predType) {
    vector<T> data;
    if (!H5Lexists(file.getId(), dataset_path.c_str(), H5P_DEFAULT)) {
        cerr << "Warning: Dataset or Link does not exist at path: " << dataset_path << endl;
        return data;
    }
    try {
        DataSet dataset = file.openDataSet(dataset_path);
        DataSpace dataspace = dataset.getSpace();
        int rank = dataspace.getSimpleExtentNdims();
        if (rank == 0) return {};
        vector<hsize_t> dims(rank);
        dataspace.getSimpleExtentDims(dims.data());
        hsize_t total_size = 1;
        for (hsize_t dim : dims) {
            total_size *= dim;
        }
        if (total_size > 0) {
            data.resize(total_size);
            dataset.read(data.data(), predType);
        }
    }
    catch (const H5::Exception& error) {
        cerr << "HDF5 error while reading dataset: " << dataset_path << endl;
        error.printErrorStack();
        return {};
    }
    return data;
}

ElementClassInfo ReadElementClassAttributes(const H5File& file, const string& group_path) {
    ElementClassInfo info;
    Group group = file.openGroup(group_path);
    info.elementType = readAttribute(group, "ElementType");
    info.sectionCategory = readAttribute(group, "SectionCatergory");
    info.isValid = true;
    return info;
}

StepClassInfo ReadStepAttributes(const H5File& file, const string& group_path) {
    StepClassInfo info;
    Group group = file.openGroup(group_path);
    info.Description = readAttribute(group, "Description");
    info.Domain = readAttribute(group, "Domain");
    info.Index = readAttribute(group, "Index");
    info.isValid = true;
    return info;
}

FrameClassInfo ReadFrameAttributes(const H5File& file, const string& group_path) {
    FrameClassInfo info;
    Group group = file.openGroup(group_path);
    info.Description = readAttribute(group, "Description");
    info.Inc_Mode = readAttribute(group, "Inc/Mode");
    info.LoadCase = readAttribute(group, "LoadCase");
    info.Time_Freq = readAttribute(group, "Time/Freq");
    info.isValid = true;
    return info;
}

FLUIDClassInfo ReadFLUIDAttributes(const H5File& file, const string& group_path) {
    FLUIDClassInfo info;
    Group group = file.openGroup(group_path);
    info.ComponentLabels = readAttribute(group, "ComponentLabels");
    info.Invariants = readAttribute(group, "Invariants");
    info.Position = readAttribute(group, "Position");
    info.Type = readAttribute(group, "Type");
    info.isValid = true;
    return info;
}

void printAttribute(const string& attrName, const AttributeDetails& details) {
    if (details.typeName.empty()) return;
    cout << "--- " << attrName << " Details ---" << endl;
    cout << "  Value: " << details.value << endl;
    cout << "  Type: " << details.typeName << endl;
    if (details.typeName == "String") {
        cout << "    - Length: " << details.length << ", Padding: " << details.padding << ", Charset: " << details.charset << endl;
    }
    cout << "  Dims: ";
    for (hsize_t d : details.dims) { cout << d << " "; }
    cout << endl;
}

void printFluidInfo(const string& name, const FLUIDClassInfo& info) {
    cout << "\n=========================================" << endl;
    cout << "  Attributes for Dataset Group: " << name << endl;
    cout << "=========================================" << endl;
    if (!info.isValid) { cout << "Invalid or unreadable." << endl; return; }
    printAttribute("ComponentLabels", info.ComponentLabels);
    printAttribute("Invariants", info.Invariants);
    printAttribute("Position", info.Position);
    printAttribute("Type", info.Type);
}

template<typename T>
void printFirstN(const vector<T>& data, const string& label, int n = 5, int precision = 1) {
    cout << "First " << n << " " << label << ": ";
    int count = min((int)data.size(), n);
    if (count == 0) {
        cout << "[empty]" << endl;
        return;
    }
    cout << fixed << setprecision(precision);
    for (int i = 0; i < count; ++i) {
        cout << data[i] << " ";
    }
    cout.unsetf(ios_base::floatfield);
    cout << "..." << endl;
}


// --- 主函数 ---

int main() {
    // 关闭 HDF5 库的自动错误打印，由我们自己的 try-catch 处理
    H5::Exception::dontPrint();

    const string file_name = "D:/data/fluid1.h5";
    const string part_name = "FLUID";
    const string step_name = "STEP-1";
    const int frame_index = 1;

    try {
        // 在程序开始时只打开一次文件
        H5File file(file_name, H5F_ACC_RDONLY);

        // --- 读取节点信息 ---
        cout << "--- Reading Node Info ---" << endl;
        auto node_labels = ReadFlattenedDataset<int>(file, "/Parts/" + part_name + "/Nodes/Labels", PredType::NATIVE_INT);
        auto flat_coords = ReadFlattenedDataset<float>(file, "/Parts/" + part_name + "/Nodes/Coordinates", PredType::NATIVE_FLOAT);
        printFirstN(node_labels, "node labels");
        cout << "Total node labels: " << node_labels.size() << endl << endl;

        cout << "First 5 node coordinates: ";
        // **【修改部分】**
        // 1. 临时设置精度为2位小数
        cout << fixed << setprecision(2);
        int num_coords_to_print = min((int)flat_coords.size() / 3, 5);
        for (int i = 0; i < num_coords_to_print; ++i) {
            cout << "(" << flat_coords[i * 3] << ", " << flat_coords[i * 3 + 1] << ", " << flat_coords[i * 3 + 2] << ") ";
        }
        // 2. 恢复默认的浮点数格式
        cout.unsetf(ios_base::floatfield);
        cout << endl;
        // **【修改结束】**

        cout << "Total coordinate sets: " << flat_coords.size() / 3 << endl << endl;

        // --- 读取单元信息 ---
        cout << "--- Reading Element Info (ElementClass:0) ---" << endl;
        string element_class_path = "/Parts/" + part_name + "/Elements/ElementClass:0";
        ElementClassInfo elem_info = ReadElementClassAttributes(file, element_class_path);
        printAttribute("ElementType", elem_info.elementType);
        printAttribute("SectionCategory", elem_info.sectionCategory);
        auto element_labels = ReadFlattenedDataset<int>(file, element_class_path + "/Labels", PredType::NATIVE_INT);
        auto connectivities = ReadFlattenedDataset<int>(file, element_class_path + "/Connectivities", PredType::NATIVE_INT);
        printFirstN(element_labels, "element labels");
        cout << "Total element labels: " << element_labels.size() << endl << endl;
        int nodes_per_element = 0;
        if (!element_labels.empty()) {
            nodes_per_element = connectivities.size() / element_labels.size();
        }
        cout << "First 5 connectivity tuples:" << endl;
        int num_elements_to_print = min((int)element_labels.size(), 5);
        for (int i = 0; i < num_elements_to_print; ++i) {
            cout << "  Element " << element_labels[i] << ": [ ";
            for (int j = 0; j < nodes_per_element; ++j) {
                cout << connectivities[i * nodes_per_element + j] << " ";
            }
            cout << "]" << endl;
        }
        cout << "\nTotal elements: " << element_labels.size() << endl;
        if (nodes_per_element > 0) {
            cout << "Nodes per element: " << nodes_per_element << endl;
        }
        cout << endl;

        // --- 读取 Step 和 Frame 的属性 ---
        cout << "--- Reading Step-1 Attributes ---" << endl;
        string step_path = "/Steps/" + step_name;
        StepClassInfo step_info = ReadStepAttributes(file, step_path);
        printAttribute("Description", step_info.Description);
        printAttribute("Domain", step_info.Domain);
        printAttribute("Index", step_info.Index);
        cout << endl;
        cout << "--- Reading Frame:" << frame_index << " Attributes ---" << endl;
        string frame_path = step_path + "/Frames/Frame:" + to_string(frame_index);
        FrameClassInfo frame_info = ReadFrameAttributes(file, frame_path);
        printAttribute("Description", frame_info.Description);
        printAttribute("Inc/Mode", frame_info.Inc_Mode);
        printAttribute("LoadCase", frame_info.LoadCase);
        printAttribute("Time/Freq", frame_info.Time_Freq);

        // --- 读取流体数据集的属性和真实数据 ---
        const vector<string> dataset_names = { "BURNF", "DENSITY", "EVF", "S", "VF" };
        for (const auto& name : dataset_names) {
            FLUIDClassInfo fluid_info = ReadFLUIDAttributes(file, frame_path + "/" + name);
            printFluidInfo(name, fluid_info);

            string data_path = frame_path + "/" + name + "/" + part_name + "-1/ElementClass:0/LocationIndex:1/Real";
            auto real_data = ReadFlattenedDataset<float>(file, data_path, PredType::NATIVE_FLOAT);
            cout << "  -> Total " << name << " dataset size: " << real_data.size() << endl;

            int num_components = 1;
            if (!fluid_info.ComponentLabels.dims.empty() && fluid_info.ComponentLabels.dims[0] > 0) {
                num_components = fluid_info.ComponentLabels.dims[0];
            }

            if (num_components > 1) {
                cout << "  > First 5 " << name << " data points:" << endl;
                int num_points_to_print = min((int)real_data.size() / num_components, 5);

                cout << fixed << setprecision(3);
                for (int i = 0; i < num_points_to_print; ++i) {
                    cout << "    [ ";
                    for (int j = 0; j < num_components; ++j) {
                        cout << real_data[i * num_components + j] << " ";
                    }
                    cout << "]" << endl;
                }
                cout.unsetf(ios_base::floatfield);
            }
            else {
                cout << "  > ";
                if (name == "DENSITY") {
                    printFirstN(real_data, name + " data", 5, 3);
                }
                else {
                    printFirstN(real_data, name + " data", 5, 1);
                }
            }
        }

    }
    catch (const H5::Exception& error) {
        cerr << "A top-level HDF5 error occurred:" << endl;
        error.printErrorStack();
        return 1;
    }

    return 0;
}
