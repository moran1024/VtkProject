#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include "H5Cpp.h"

// Using declarations for convenience
using namespace std;
using namespace H5;

// Struct to hold 3D coordinate data
struct Triplet {
    float x, y, z;
};

// Struct to hold Element Class attribute information
struct ElementClassInfo {
    string elementType;
    string sectionCategory;
    bool isValid = false; // Flag to check if read was successful
};

// --- Reader Functions ---

vector<int> Read_Int_Data_From_File(const string& file_name, const string& dataset_path) {
    vector<int> data;
    try {
        H5File file(file_name, H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet(dataset_path);
        DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims);
        data.resize(dims[0]);
        dataset.read(data.data(), PredType::NATIVE_INT);
    }
    catch (const H5::Exception& error) {
        error.printErrorStack();
        return {};
    }
    return data;
}

vector<Triplet> Read_Coordinates_From_File(const string& file_name, const string& dataset_path) {
    vector<Triplet> coordinates;
    try {
        H5File file(file_name, H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet(dataset_path);
        DataSpace dataspace = dataset.getSpace();
        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims);
        coordinates.resize(dims[0]);
        dataset.read(coordinates.data(), PredType::NATIVE_FLOAT);
    }
    catch (const H5::Exception& error) {
        error.printErrorStack();
        return {};
    }
    return coordinates;
}

vector<vector<int>> Read_Connectivities_From_File(const string& file_name, const string& dataset_path) {
    vector<vector<int>> connectivities_2d;
    try {
        H5File file(file_name, H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet(dataset_path);
        DataSpace dataspace = dataset.getSpace();
        if (dataspace.getSimpleExtentNdims() != 2) {
            cerr << "Error: Connectivities dataset is not 2D." << endl;
            return {};
        }
        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims);
        vector<int> flat_connectivities(dims[0] * dims[1]);
        dataset.read(flat_connectivities.data(), PredType::NATIVE_INT);
        connectivities_2d.resize(dims[0]);
        for (hsize_t i = 0; i < dims[0]; ++i) {
            auto start_it = flat_connectivities.begin() + i * dims[1];
            auto end_it = start_it + dims[1];
            connectivities_2d[i] = vector<int>(start_it, end_it);
        }
    }
    catch (const H5::Exception& error) {
        error.printErrorStack();
        return {};
    }
    return connectivities_2d;
}

ElementClassInfo ReadElementClassAttributes(const string& file_name, const string& part_name, int class_index) {
    ElementClassInfo info;
    // CORRECTED: Path construction logic is fixed.
    string group_path = "/Parts/" + part_name + "/Elements/ElementClass:" + to_string(class_index);

    try {
        H5File file(file_name, H5F_ACC_RDONLY);
        Group group = file.openGroup(group_path);

        Attribute elementTypeAttr = group.openAttribute("ElementType");
        DataType type = elementTypeAttr.getDataType();
        elementTypeAttr.read(type, info.elementType);

        Attribute sectionCategoryAttr = group.openAttribute("SectionCategory");
        DataType type2 = sectionCategoryAttr.getDataType();
        sectionCategoryAttr.read(type2, info.sectionCategory);

        info.isValid = true;
    }
    catch (const H5::Exception& error) {
        cerr << "HDF5 Error reading attributes from " << group_path << ":" << endl;
        error.printErrorStack();
        return {};
    }
    return info;
}

// --- Main Execution ---

int main() {
    const string file_name = "D:/data/fluid0.h5";
    const string part_name = "FLUID";

    // --- Read Node Information ---
    cout << "--- Nodes ---" << endl;
    vector<int> node_labels = Read_Int_Data_From_File(file_name, "/Parts/" + part_name + "/Nodes/Labels");
    cout << "First 5 node labels: ";
    for (int i = 0; i < 5 && i < node_labels.size(); i++) {
        cout << node_labels[i] << " ";
    }
    cout << "\nTotal node labels: " << node_labels.size() << endl << endl;

    vector<Triplet> coordinates = Read_Coordinates_From_File(file_name, "/Parts/" + part_name + "/Nodes/Coordinates");
    cout << "First 5 node coordinates: ";
    for (int i = 0; i < 5 && i < coordinates.size(); i++) {
        cout << "(" << coordinates[i].x << ", " << coordinates[i].y << ", " << coordinates[i].z << ") ";
    }
    cout << "\nTotal coordinate sets: " << coordinates.size() << endl << endl;

    // --- Read Element Information ---
    int class_index = 0;
    cout << "--- Elements (in ElementClass:" << class_index << ") ---" << endl;

    // ADDED: Call the function to read attributes
    ElementClassInfo attributes = ReadElementClassAttributes(file_name, part_name, class_index);
    if (attributes.isValid) {
        cout << "Attributes:" << endl;
        cout << "  - ElementType: " << attributes.elementType << endl;
        cout << "  - SectionCategory: " << attributes.sectionCategory << endl << endl;
    }

    string element_class_path = "/Parts/" + part_name + "/Elements/ElementClass:" + to_string(class_index);
    vector<int> element_labels = Read_Int_Data_From_File(file_name, element_class_path + "/Labels");
    cout << "First 5 element labels: ";
    for (int i = 0; i < 5 && i < element_labels.size(); i++) {
        cout << element_labels[i] << " ";
    }
    cout << "\nTotal element labels: " << element_labels.size() << endl << endl;

    vector<vector<int>> connectivities = Read_Connectivities_From_File(file_name, element_class_path + "/Connectivities");
    cout << "First 5 connectivity tuples:" << endl;
    for (int i = 0; i < 5 && i < connectivities.size(); ++i) {
        cout << "  Element " << element_labels[i] << ": [ ";
        for (int node_id : connectivities[i]) {
            cout << node_id << " ";
        }
        cout << "]" << endl;
    }
    cout << "\nTotal elements: " << connectivities.size() << endl;
    if (!connectivities.empty()) {
        cout << "Nodes per element: " << connectivities[0].size() << endl;
    }

    return 0;
}
