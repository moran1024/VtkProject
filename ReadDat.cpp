#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
#include <limits>

using namespace std;


struct InternalMesh {
    long long numberOfPoints = 0;
    long long numberOfElements = 0;
    std::vector<std::string> variableNames;
    std::vector<double> points; // 存储 (x1,y1,z1, x2,y2,z2, ...)
    std::map<std::string, std::vector<double>> pointData;
    std::vector<long long> connectivity; // 存储单元连接信息
    std::string title;
    std::string zoneName;
    int strandID = 0;
    double solutionTime = 0.0;
    std::string zoneType;
    std::string dataPacking;
};

class Dat_Reader {
public:
    bool readFile(const std::string& filePath);
    const InternalMesh& getMeshData() const { return m_meshData; }

private:
    InternalMesh m_meshData;
    void trim(std::string& s);
    std::vector<std::string> extractQuotedStrings(std::string& line);
    void parseZoneHeader(const std::string& buffer);
    bool readBlockData(std::ifstream& in);
};


void Dat_Reader::trim(std::string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r"));
    s.erase(s.find_last_not_of(" \t\n\r") + 1);
}

std::vector<std::string> Dat_Reader::extractQuotedStrings(std::string& line) {
    std::vector<std::string> result;
    std::stringstream ss(line);
    std::string segment;
    while (std::getline(ss, segment, '"')) {
        if (std::getline(ss, segment, '"')) {
            result.push_back(segment);
        }
    }
    return result;
}


void Dat_Reader::parseZoneHeader(const std::string& buffer) {
    std::string processed_buffer = buffer;

    std::string upper_buffer = buffer;
    std::transform(upper_buffer.begin(), upper_buffer.end(), upper_buffer.begin(), ::toupper);
    size_t zone_pos = upper_buffer.find("ZONE");

    if (zone_pos != std::string::npos) {
        processed_buffer.erase(zone_pos, 4);
    }


    std::stringstream ss(processed_buffer); 
    std::string segment;
    while (std::getline(ss, segment, ',')) {
        size_t eq_pos = segment.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = segment.substr(0, eq_pos);
            std::string value = segment.substr(eq_pos + 1);
            trim(key);
            trim(value);
            std::transform(key.begin(), key.end(), key.begin(), ::toupper);
            value.erase(remove(value.begin(), value.end(), '"'), value.end());

            if (key == "T") m_meshData.zoneName = value;
            else if (key == "NODES" || key == "N") m_meshData.numberOfPoints = std::stoll(value); // 兼容简写 N
            else if (key == "ELEMENTS" || key == "E") m_meshData.numberOfElements = std::stoll(value); // 兼容简写 E
            else if (key == "STRANDID") m_meshData.strandID = std::stoi(value);
            else if (key == "SOLUTIONTIME") m_meshData.solutionTime = std::stod(value);
            else if (key == "ZONETYPE") m_meshData.zoneType = value;
            else if (key == "DATAPACKING") m_meshData.dataPacking = value;
        }
    }
}


bool Dat_Reader::readFile(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return false;
    }

    std::string line;
    std::string zoneHeaderBuffer;
    bool inZoneHeader = false;

    while (true) { 
        std::streampos pos_before_getline = file.tellg();

        if (!std::getline(file, line)) {
            break; 
        }

        trim(line);
        if (line.empty()) continue;

        std::string upperLine = line;
        std::transform(upperLine.begin(), upperLine.end(), upperLine.begin(), ::toupper);

        if (upperLine.rfind("TITLE", 0) == 0) {
            line.erase(0, line.find('=') + 1);
            trim(line);
            line.erase(remove(line.begin(), line.end(), '"'), line.end());
            m_meshData.title = line;
        }
        else if (upperLine.rfind("VARIABLES", 0) == 0) {
            line.erase(0, line.find('=') + 1);
            auto vars = extractQuotedStrings(line);
            m_meshData.variableNames.insert(m_meshData.variableNames.end(), vars.begin(), vars.end());

            while (file.peek() == '"' || (file.peek() == ' ' && file.seekg(1, ios_base::cur) && file.peek() == '"' && file.seekg(-1, ios_base::cur))) {
                std::getline(file, line);
                trim(line);
                vars = extractQuotedStrings(line);
                m_meshData.variableNames.insert(m_meshData.variableNames.end(), vars.begin(), vars.end());
            }
        }
        else if (upperLine.rfind("ZONE", 0) == 0) {
            inZoneHeader = true;
            zoneHeaderBuffer = line;
        }
        else if (inZoneHeader) {
            if (line.find('=') != std::string::npos) {
                zoneHeaderBuffer += "," + line;
            }
            else {
                inZoneHeader = false;
                parseZoneHeader(zoneHeaderBuffer);

                file.seekg(pos_before_getline);

                if (!readBlockData(file)) return false;
                break; 
            }
        }
    }
    file.close();
    return true;
}

bool Dat_Reader::readBlockData(std::ifstream& in) {
    long long num_nodes = m_meshData.numberOfPoints;
    if (num_nodes == 0 || m_meshData.variableNames.empty()) {
        std::cerr << "Error: Node count or variables not defined before data block." << std::endl;
        return false;
    }
    if (m_meshData.dataPacking != "BLOCK") {
        std::cerr << "Error: Only DATAPACKING=BLOCK is supported." << std::endl;
        return false;
    }

    std::cout << "Reading BLOCK data for " << num_nodes << " nodes." << std::endl;
    std::vector<std::vector<double>> data_columns(m_meshData.variableNames.size());

    for (size_t i = 0; i < m_meshData.variableNames.size(); ++i) {
        data_columns[i].reserve(num_nodes);
        for (long long j = 0; j < num_nodes; ++j) {
            float value;
            in >> value;
            if (in.fail()) {
                std::cerr << "Error reading data for variable '" << m_meshData.variableNames[i]
                    << "', node " << j + 1 << std::endl;
                return false;
            }
            data_columns[i].push_back(static_cast<double>(value));
        }
    }
    std::cout << "--> Finished reading all raw data columns." << std::endl;

    m_meshData.points.reserve(num_nodes * 3);
    const auto& x_coords = data_columns[0];
    const auto& y_coords = data_columns[1];
    const auto& z_coords = data_columns[2];
    for (long long i = 0; i < num_nodes; ++i) {
        m_meshData.points.push_back(x_coords[i]);
        m_meshData.points.push_back(y_coords[i]);
        m_meshData.points.push_back(z_coords[i]);
    }

    for (size_t i = 3; i < m_meshData.variableNames.size(); ++i) {
        m_meshData.pointData[m_meshData.variableNames[i]] = data_columns[i];
    }

    if (m_meshData.numberOfElements > 0) {
        std::cout << "Reading connectivity for " << m_meshData.numberOfElements << " elements." << std::endl;
        int nodesPerElement = 0;
        if (m_meshData.zoneType == "FEBrick") nodesPerElement = 8;
        else if (m_meshData.zoneType == "FETetrahedron") nodesPerElement = 4;
        else if (m_meshData.zoneType == "FEQuad") nodesPerElement = 4;
        else if (m_meshData.zoneType == "FETriangle") nodesPerElement = 3;
        else std::cerr << "Warning: Unsupported ZONETYPE for connectivity reading: " << m_meshData.zoneType << std::endl;

        if (nodesPerElement > 0) {
            long long totalConnectivityEntries = m_meshData.numberOfElements * nodesPerElement;
            m_meshData.connectivity.reserve(totalConnectivityEntries);
            for (long long i = 0; i < totalConnectivityEntries; ++i) {
                long long nodeIndex;
                in >> nodeIndex;
                if (in.fail()) {
                    std::cerr << "Error reading connectivity data at entry " << i << std::endl;
                    return false;
                }
                m_meshData.connectivity.push_back(nodeIndex - 1);
            }
            std::cout << "--> Finished reading connectivity." << std::endl;
        }
    }
    return true;
}

void printDatInfo(const InternalMesh& mesh) {
    cout << "\n--- Tecplot DAT File Summary ---" << endl;
    cout << "Title: " << mesh.title << endl;
    cout << "Zone Name: " << mesh.zoneName << endl;
    cout << "Solution Time: " << mesh.solutionTime << endl;
    cout << "Strand ID: " << mesh.strandID << endl;
    cout << "Zone Type: " << mesh.zoneType << endl;
    cout << "Data Packing: " << mesh.dataPacking << endl;
    cout << "Variables (" << mesh.variableNames.size() << "): ";
    for (const auto& name : mesh.variableNames) cout << name << " | ";
    cout << endl;
    cout << "Number of Points: " << mesh.numberOfPoints << endl;
    cout << "Number of Elements: " << mesh.numberOfElements << endl;

    if (!mesh.connectivity.empty() && mesh.numberOfElements > 0) {
        long long nodesPerElement_ll = mesh.connectivity.size() / mesh.numberOfElements;
        int nodesPerElement = static_cast<int>(nodesPerElement_ll);
        cout << "Connectivity for first element (" << nodesPerElement << " nodes): ";
        for (int i = 0; i < nodesPerElement && i < mesh.connectivity.size(); ++i) {
            cout << mesh.connectivity[i] << " ";
        }
        cout << endl;
    }

    cout << "\n--- First 10 Coordinate Values ---" << endl;
    long long points_to_print = std::min(10LL, mesh.numberOfPoints);

    if (points_to_print > 0) {
        cout << "X: ";
        for (long long i = 0; i < points_to_print; ++i) {
            cout << mesh.points[i * 3] << " ";
        }
        cout << endl;

        cout << "Y: ";
        for (long long i = 0; i < points_to_print; ++i) {
            cout << mesh.points[i * 3 + 1] << " ";
        }
        cout << endl;

        cout << "Z: ";
        for (long long i = 0; i < points_to_print; ++i) {
            cout << mesh.points[i * 3 + 2] << " ";
        }
        cout << endl;
    }
    else {
        cout << "No point data available to print." << endl;
    }

    if (!mesh.pointData.empty()) {
        auto const& first_field = mesh.pointData.begin();
        cout << "\nFirst data field '" << first_field->first << "' has " << first_field->second.size() << " values." << endl;
    }
    cout << "---------------------------------" << endl;
}


int main() {
    Dat_Reader reader;
    
    std::string filePath = "D:/桌面/数据读取课题/数据文件/tecplotdat文件/flowcontour/flowcontour      0.dat"; // 例如: "D:/data/flow.dat"
    if (reader.readFile(filePath)) {
        std::cout << "\nFile reading process completed successfully." << std::endl;
        const auto& meshData = reader.getMeshData();
        printDatInfo(meshData);
    }
    else {
        std::cout << "\nFile reading process failed." << std::endl;
    }
    return 0;
}


