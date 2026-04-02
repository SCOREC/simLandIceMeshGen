#include <netcdfWriter.h>
#include <netcdf.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cstring> //strlen

#include "MeshSim.h"
#include "SimParasolidKrnl.h"
#include "SimMessages.h"
#include "SimInfo.h"
#include "SimInfoCodes.h"
#include "SimMeshingInfoCodes.h"

// Error handling macro for NetCDF C API
#define NC_CHECK(e) { \
  int nc_stat = (e); \
  if (nc_stat != NC_NOERR) { \
    std::cerr << "NetCDF Error: " << nc_strerror(nc_stat) << std::endl; \
    return 1; \
  } \
}

int writeMeshSimToNetCDF(pMesh mesh, pGModel model, std::string outputFileName) {
  // Get mesh dimensions
  int numVertices = M_numVertices(mesh);
  int numCells = M_numFaces(mesh);  // In 2D, cells are faces
  int vertexDegree = 3;  // Maximum number of cells adjacent to a vertex

  // Create the NetCDF file
  int ncid;
  NC_CHECK(nc_create(outputFileName.c_str(), NC_CLOBBER, &ncid));

  // Define dimensions
  int nVerticesDimID, nCellsDimID, vertexDegreeDimID;
  NC_CHECK(nc_def_dim(ncid, "nVertices", numVertices, &nVerticesDimID));
  NC_CHECK(nc_def_dim(ncid, "nCells", numCells, &nCellsDimID));
  NC_CHECK(nc_def_dim(ncid, "vertexDegree", vertexDegree, &vertexDegreeDimID));

  // Define variables
  int xVertexVarID, yVertexVarID, zVertexVarID;
  int xCellVarID, yCellVarID, zCellVarID;
  int cellsOnVertexVarID;
  int meshDensityVarID;
  int geomModelIdCellVarID, geomModelDimCellVarID;
  int geomModelIdVertexVarID, geomModelDimVertexVarID;

  NC_CHECK(nc_def_var(ncid, "xVertex", NC_DOUBLE, 1, &nVerticesDimID, &xVertexVarID));
  NC_CHECK(nc_def_var(ncid, "yVertex", NC_DOUBLE, 1, &nVerticesDimID, &yVertexVarID));
  NC_CHECK(nc_def_var(ncid, "zVertex", NC_DOUBLE, 1, &nVerticesDimID, &zVertexVarID));

  NC_CHECK(nc_def_var(ncid, "xCell", NC_DOUBLE, 1, &nCellsDimID, &xCellVarID));
  NC_CHECK(nc_def_var(ncid, "yCell", NC_DOUBLE, 1, &nCellsDimID, &yCellVarID));
  NC_CHECK(nc_def_var(ncid, "zCell", NC_DOUBLE, 1, &nCellsDimID, &zCellVarID));

  int cellsOnVertexDimIDs[2] = {nVerticesDimID, vertexDegreeDimID};
  NC_CHECK(nc_def_var(ncid, "cellsOnVertex", NC_INT, 2, cellsOnVertexDimIDs, &cellsOnVertexVarID));

  NC_CHECK(nc_def_var(ncid, "meshDensity", NC_DOUBLE, 1, &nCellsDimID, &meshDensityVarID));

  NC_CHECK(nc_def_var(ncid, "geomModelIdCell", NC_INT, 1, &nCellsDimID, &geomModelIdCellVarID));
  NC_CHECK(nc_def_var(ncid, "geomModelDimCell", NC_INT, 1, &nCellsDimID, &geomModelDimCellVarID));
  NC_CHECK(nc_def_var(ncid, "geomModelIdVertex", NC_INT, 1, &nVerticesDimID, &geomModelIdVertexVarID));
  NC_CHECK(nc_def_var(ncid, "geomModelDimVertex", NC_INT, 1, &nVerticesDimID, &geomModelDimVertexVarID));

  // Add global attributes
  const char* on_a_sphere_val = "NO";
  NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "on_a_sphere", strlen(on_a_sphere_val), on_a_sphere_val));
  double sphere_radius_val = 0.0;
  NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &sphere_radius_val));

  // End define mode
  NC_CHECK(nc_enddef(ncid));

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
  NC_CHECK(nc_put_var_double(ncid, xVertexVarID, xVertex.data()));
  NC_CHECK(nc_put_var_double(ncid, yVertexVarID, yVertex.data()));
  NC_CHECK(nc_put_var_double(ncid, zVertexVarID, zVertex.data()));
  NC_CHECK(nc_put_var_int(ncid, geomModelIdVertexVarID, geomModelIdVertex.data()));
  NC_CHECK(nc_put_var_int(ncid, geomModelDimVertexVarID, geomModelDimVertex.data()));

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
  NC_CHECK(nc_put_var_double(ncid, xCellVarID, xCell.data()));
  NC_CHECK(nc_put_var_double(ncid, yCellVarID, yCell.data()));
  NC_CHECK(nc_put_var_double(ncid, zCellVarID, zCell.data()));
  NC_CHECK(nc_put_var_double(ncid, meshDensityVarID, meshDensity.data()));
  NC_CHECK(nc_put_var_int(ncid, geomModelIdCellVarID, geomModelIdCell.data()));
  NC_CHECK(nc_put_var_int(ncid, geomModelDimCellVarID, geomModelDimCell.data()));

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
  NC_CHECK(nc_put_var_int(ncid, cellsOnVertexVarID, cellsOnVertex.data()));

  // Close the NetCDF file
  NC_CHECK(nc_close(ncid));

  std::cout << "Successfully wrote NetCDF file: " << outputFileName << std::endl;
  std::cout << "  Number of vertices: " << numVertices << std::endl;
  std::cout << "  Number of cells: " << numCells << std::endl;

  return 0;
}
