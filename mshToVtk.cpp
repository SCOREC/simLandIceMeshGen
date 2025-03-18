#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

/* 
initial broken (bad input file parsing) version written by the free version of chatgpt in March 2025 given the following prompt:
> write C++ code that converts a jigsaw .msh file to a paraview legacy vtk file

the response was:
> Here's a C++ program that reads a .msh file (from JIGSAW) and converts it into a ParaView legacy .vtk file. It assumes that:
> 
>     The first block contains vertex coordinates.
>     The second block contains edges, indexed into the list of vertices.
>     Indices in the .msh file are 1-based and need to be converted to 0-based for VTK.
> 
> The program reads the input .msh file, processes the vertices and edges, and writes the output in the VTK legacy format.
*/

struct Vertex {
    double x, y, z;
};

struct Edge {
    int v1, v2;
};

void convertMshToVtk(const std::string& inputFile, const std::string& outputFile) {
    std::ifstream mshFile(inputFile);
    std::ofstream vtkFile(outputFile);
    
    if (!mshFile.is_open() || !vtkFile.is_open()) {
        std::cerr << "Error opening file!\n";
        return;
    }
    
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::string line;
    bool readingVertices = true;

    while (std::getline(mshFile, line)) {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        
        while (std::getline(iss, token, ';')) {
          tokens.push_back(token);
        }
        
        if (tokens[0][0] == '#') continue;
        if (tokens.size() == 2 && ! tokens[0].compare(0,5,"MSHID") ) continue;
        if (tokens.size() == 1 && ! tokens[0].compare(0,4,"NDIMS") ) continue;
        if (tokens.size() == 1 && ! tokens[0].compare(0,5,"POINT") ) continue;
        
        if (tokens.size() == 1 && ! tokens[0].compare(0,5,"EDGE2")) { // Transition to edges
            readingVertices = false;
            continue;
        }
        
        if (readingVertices) {
            if (tokens.size() >= 2) { // 2D (x, y) or 3D (x, y, z)
                Vertex v;
                v.x = std::stod(tokens[0]);
                v.y = std::stod(tokens[1]);
                v.z = (tokens.size() > 2) ? std::stod(tokens[2]) : 0.0;
                vertices.push_back(v);
            }
        } else {
            if (tokens.size() >= 2) { // Edge indices
                Edge e;
                e.v1 = std::stoi(tokens[0]);
                e.v2 = std::stoi(tokens[1]);
                edges.push_back(e);
            }
        }
    }

    mshFile.close();
    
    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Converted from .msh format\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";
    
    // Write points
    vtkFile << "POINTS " << vertices.size() << " float\n";
    for (const auto& v : vertices) {
        vtkFile << v.x << " " << v.y << " " << v.z << "\n";
    }
    
    // Write lines
    vtkFile << "LINES " << edges.size() << " " << edges.size() * 3 << "\n";
    for (const auto& e : edges) {
        vtkFile << "2 " << e.v1 << " " << e.v2 << "\n";
    }
    
    vtkFile.close();
    std::cout << "Conversion complete: " << outputFile << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.msh output.vtk\n";
        return 1;
    }
    
    convertMshToVtk(argv[1], argv[2]);
    return 0;
}

