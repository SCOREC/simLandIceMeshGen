#include <netcdfWriter.h>
#include <netcdf>
#include <iostream>
#include <vector>
#include <map>
#include <set>

#include "MeshSim.h"
#include "SimParasolidKrnl.h"
#include "SimMessages.h"
#include "SimInfo.h"
#include "SimInfoCodes.h"
#include "SimMeshingInfoCodes.h"

using namespace netCDF;
using namespace netCDF::exceptions;

int writeMeshSimToNetCDF(pMesh mesh, pGModel model, std::string outputFileName) {
  try {
    // Create the NetCDF file
    NcFile dataFile(outputFileName, NcFile::replace);

    // Get mesh dimensions
    int numVertices = M_numVertices(mesh);
    int numCells = M_numFaces(mesh);  // In 2D, cells are faces
    int vertexDegree = 3;  // Maximum number of cells adjacent to a vertex

    // Define dimensions
    NcDim nVerticesDim = dataFile.addDim("nVertices", numVertices);
    NcDim nCellsDim = dataFile.addDim("nCells", numCells);
    NcDim vertexDegreeDim = dataFile.addDim("vertexDegree", vertexDegree);

    // Define variables
    NcVar xVertexVar = dataFile.addVar("xVertex", ncDouble, nVerticesDim);
    NcVar yVertexVar = dataFile.addVar("yVertex", ncDouble, nVerticesDim);
    NcVar zVertexVar = dataFile.addVar("zVertex", ncDouble, nVerticesDim);

    NcVar xCellVar = dataFile.addVar("xCell", ncDouble, nCellsDim);
    NcVar yCellVar = dataFile.addVar("yCell", ncDouble, nCellsDim);
    NcVar zCellVar = dataFile.addVar("zCell", ncDouble, nCellsDim);

    std::vector<NcDim> cellsOnVertexDims;
    cellsOnVertexDims.push_back(nVerticesDim);
    cellsOnVertexDims.push_back(vertexDegreeDim);
    NcVar cellsOnVertexVar = dataFile.addVar("cellsOnVertex", ncInt, cellsOnVertexDims);

    NcVar meshDensityVar = dataFile.addVar("meshDensity", ncDouble, nCellsDim);

    NcVar geomModelIdCellVar = dataFile.addVar("geomModelIdCell", ncInt, nCellsDim);
    NcVar geomModelDimCellVar = dataFile.addVar("geomModelDimCell", ncInt, nCellsDim);
    NcVar geomModelIdVertexVar = dataFile.addVar("geomModelIdVertex", ncInt, nVerticesDim);
    NcVar geomModelDimVertexVar = dataFile.addVar("geomModelDimVertex", ncInt, nVerticesDim);

    // Add global attributes
    dataFile.putAtt("on_a_sphere", "NO");
    dataFile.putAtt("sphere_radius", ncDouble, 0.0);

    // Allocate arrays for vertex data
    std::vector<double> xVertex(numVertices);
    std::vector<double> yVertex(numVertices);
    std::vector<double> zVertex(numVertices, 0.0);  // Set to 0 for planar mesh
    std::vector<int> geomModelIdVertex(numVertices);
    std::vector<int> geomModelDimVertex(numVertices);

    // Map to store vertex index by vertex pointer
    std::map<pVertex, int> vertexIndexMap;

    // Iterate over vertices and renumber them
    VIter vertices = M_vertexIter(mesh);
    pVertex vertex;
    int vertexIdx = 0;
    while ((vertex = VIter_next(vertices))) {
      double xyz[3];
      V_coord(vertex, xyz);

      xVertex[vertexIdx] = xyz[0];
      yVertex[vertexIdx] = xyz[1];
      zVertex[vertexIdx] = 0.0;  // 2D mesh, z = 0

      // Get geometric classification
      pGEntity gent = EN_whatIn((pEntity)vertex);
      geomModelIdVertex[vertexIdx] = GEN_tag(gent);
      geomModelDimVertex[vertexIdx] = GEN_type(gent);

      // Store mapping and renumber
      vertexIndexMap[vertex] = vertexIdx;
      EN_setID((pEntity)vertex, vertexIdx);

      vertexIdx++;
    }
    VIter_delete(vertices);

    // Write vertex coordinate data
    xVertexVar.putVar(xVertex.data());
    yVertexVar.putVar(yVertex.data());
    zVertexVar.putVar(zVertex.data());
    geomModelIdVertexVar.putVar(geomModelIdVertex.data());
    geomModelDimVertexVar.putVar(geomModelDimVertex.data());

    // Allocate arrays for cell data
    std::vector<double> xCell(numCells);
    std::vector<double> yCell(numCells);
    std::vector<double> zCell(numCells, 0.0);  // Set to 0 for planar mesh
    std::vector<double> meshDensity(numCells);
    std::vector<int> geomModelIdCell(numCells);
    std::vector<int> geomModelDimCell(numCells);

    // Map to store cell index by face pointer
    std::map<pFace, int> faceIndexMap;

    // Iterate over all faces (cells in 2D)
    FIter faces = M_faceIter(mesh);
    pFace face;
    int cellIdx = 0;
    while ((face = FIter_next(faces))) {
      // Calculate cell center
      double center[3] = {0.0, 0.0, 0.0};
      pPList faceVerts = F_vertices(face, 1);
      int numFaceVerts = PList_size(faceVerts);

      void* iter = 0;
      pVertex fv;
      while ((fv = (pVertex)PList_next(faceVerts, &iter))) {
        double xyz[3];
        V_coord(fv, xyz);
        center[0] += xyz[0];
        center[1] += xyz[1];
      }
      center[0] /= numFaceVerts;
      center[1] /= numFaceVerts;

      xCell[cellIdx] = center[0];
      yCell[cellIdx] = center[1];
      zCell[cellIdx] = 0.0;  // 2D mesh, z = 0

      // Calculate mesh density (average edge length for this face)
      pPList faceEdges = F_edges(face, 1, 0);
      int numEdges = PList_size(faceEdges);
      double totalEdgeLength = 0.0;
      void* edgeIter = 0;
      pEdge edge;
      while ((edge = (pEdge)PList_next(faceEdges, &edgeIter))) {
        totalEdgeLength += E_length(edge);
      }
      meshDensity[cellIdx] = (numEdges > 0) ? (totalEdgeLength / numEdges) : 0.0;
      PList_delete(faceEdges);

      // Get geometric classification
      pGEntity gent = EN_whatIn((pEntity)face);
      geomModelIdCell[cellIdx] = GEN_tag(gent);
      geomModelDimCell[cellIdx] = GEN_type(gent);

      // Store mapping
      faceIndexMap[face] = cellIdx;
      PList_delete(faceVerts);
      cellIdx++;
    }
    FIter_delete(faces);

    // Write cell data
    xCellVar.putVar(xCell.data());
    yCellVar.putVar(yCell.data());
    zCellVar.putVar(zCell.data());
    meshDensityVar.putVar(meshDensity.data());
    geomModelIdCellVar.putVar(geomModelIdCell.data());
    geomModelDimCellVar.putVar(geomModelDimCell.data());

    // Build cellsOnVertex connectivity (1-based indexing)
    // Initialize with 0 (indicating no cell)
    std::vector<int> cellsOnVertex(numVertices * vertexDegree, 0);

    // For each vertex, find adjacent cells
    vertices = M_vertexIter(mesh);
    while ((vertex = VIter_next(vertices))) {
      int vIdx = vertexIndexMap[vertex];

      // Get faces (cells) adjacent to this vertex
      pPList adjFaces = V_faces(vertex);
      int numAdjFaces = PList_size(adjFaces);

      void* iter = 0;
      pFace adjFace;
      int localCellIdx = 0;
      while ((adjFace = (pFace)PList_next(adjFaces, &iter)) && localCellIdx < vertexDegree) {
        int cellId = faceIndexMap[adjFace];
        // Store as 1-based index
        cellsOnVertex[vIdx * vertexDegree + localCellIdx] = cellId + 1;
        localCellIdx++;
      }

      PList_delete(adjFaces);
    }
    VIter_delete(vertices);

    // Write cellsOnVertex connectivity
    cellsOnVertexVar.putVar(cellsOnVertex.data());

    std::cout << "Successfully wrote NetCDF file: " << outputFileName << std::endl;
    std::cout << "  Number of vertices: " << numVertices << std::endl;
    std::cout << "  Number of cells: " << numCells << std::endl;

    // File is automatically closed when dataFile goes out of scope
    return 0;

  } catch (NcException& e) {
    std::cerr << "NetCDF Error: " << e.what() << std::endl;
    return 1;
  }
}
