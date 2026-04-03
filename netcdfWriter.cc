#include <netcdfWriter.h>
#include <netcdf.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cstring> //strlen
#include <cassert> //assert

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
  const int numDualVertices = M_numFaces(mesh);
  const int numDualCells = M_numVertices(mesh);
  const int dualVertexDegree = 3;

  // Create the NetCDF file
  int ncid;
  NC_CHECK(nc_create(outputFileName.c_str(), NC_CLOBBER, &ncid));

  // Define dimensions
  int nDualVerticesDimID, nDualCellsDimID, dualVertexDegreeDimID;
  NC_CHECK(nc_def_dim(ncid, "nVertices", numDualVertices, &nDualVerticesDimID));
  NC_CHECK(nc_def_dim(ncid, "nCells", numDualCells, &nDualCellsDimID));
  NC_CHECK(nc_def_dim(ncid, "vertexDegree", dualVertexDegree, &dualVertexDegreeDimID));

  // Define variables
  int xDualVertexVarID, yDualVertexVarID, zDualVertexVarID;
  int xDualCellVarID, yDualCellVarID, zDualCellVarID;
  int dualCellsOnDualVertexVarID; //primal elements to primal vertices
  int meshDensityVarID;
  int geomModelIdDualCellVarID, geomModelDimDualCellVarID;
  int geomModelIdDualVertexVarID, geomModelDimDualVertexVarID;

  NC_CHECK(nc_def_var(ncid, "xVertex", NC_DOUBLE, 1, &nDualVerticesDimID, &xDualVertexVarID));
  NC_CHECK(nc_def_var(ncid, "yVertex", NC_DOUBLE, 1, &nDualVerticesDimID, &yDualVertexVarID));
  NC_CHECK(nc_def_var(ncid, "zVertex", NC_DOUBLE, 1, &nDualVerticesDimID, &zDualVertexVarID));

  NC_CHECK(nc_def_var(ncid, "xCell", NC_DOUBLE, 1, &nDualCellsDimID, &xDualCellVarID));
  NC_CHECK(nc_def_var(ncid, "yCell", NC_DOUBLE, 1, &nDualCellsDimID, &yDualCellVarID));
  NC_CHECK(nc_def_var(ncid, "zCell", NC_DOUBLE, 1, &nDualCellsDimID, &zDualCellVarID));

  int dualCellsOnDualVertexDimIDs[2] = {nDualVerticesDimID, dualVertexDegreeDimID};
  NC_CHECK(nc_def_var(ncid, "cellsOnVertex", NC_INT, 2, dualCellsOnDualVertexDimIDs, &dualCellsOnDualVertexVarID));

  NC_CHECK(nc_def_var(ncid, "meshDensity", NC_DOUBLE, 1, &nDualCellsDimID, &meshDensityVarID));

  NC_CHECK(nc_def_var(ncid, "geomModelIdCell", NC_INT, 1, &nDualCellsDimID, &geomModelIdDualCellVarID));
  NC_CHECK(nc_def_var(ncid, "geomModelDimCell", NC_INT, 1, &nDualCellsDimID, &geomModelDimDualCellVarID));
  NC_CHECK(nc_def_var(ncid, "geomModelIdVertex", NC_INT, 1, &nDualVerticesDimID, &geomModelIdDualVertexVarID));
  NC_CHECK(nc_def_var(ncid, "geomModelDimVertex", NC_INT, 1, &nDualVerticesDimID, &geomModelDimDualVertexVarID));

  // Add global attributes
  const char* on_a_sphere_val = "NO"; //only considering planar meshes for landice at the moment
  NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "on_a_sphere", strlen(on_a_sphere_val), on_a_sphere_val));
  double sphere_radius_val = 0.0;
  NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &sphere_radius_val));

  // End define mode
  NC_CHECK(nc_enddef(ncid));

  // Allocate arrays for vertex data
  std::vector<double> xDualCell(numDualCells);
  std::vector<double> yDualCell(numDualCells);
  std::vector<double> zDualCell(numDualCells, 0.0);  // Set to 0 for planar mesh
  std::vector<int> geomModelIdDualCell(numDualCells);
  std::vector<int> geomModelDimDualCell(numDualCells);

  // Iterate over vertices and renumber them
  VIter vertices = M_vertexIter(mesh);
  pVertex vertex;
  int vertexIdx = 0;
  while ((vertex = VIter_next(vertices))) {
    double xyz[3];
    V_coord(vertex, xyz);

    xDualCell[vertexIdx] = xyz[0];
    yDualCell[vertexIdx] = xyz[1];
    zDualCell[vertexIdx] = 0.0;  // 2D mesh, z = 0

    // Get geometric classification
    pGEntity gent = EN_whatIn((pEntity)vertex);
    geomModelIdDualCell[vertexIdx] = GEN_tag(gent);
    geomModelDimDualCell[vertexIdx] = GEN_type(gent);

    vertexIdx++;
  }
  VIter_delete(vertices);

  // write dual cell data
  NC_CHECK(nc_put_var_double(ncid, xDualCellVarID, xDualCell.data()));
  NC_CHECK(nc_put_var_double(ncid, yDualCellVarID, yDualCell.data()));
  NC_CHECK(nc_put_var_double(ncid, zDualCellVarID, zDualCell.data()));
  NC_CHECK(nc_put_var_int(ncid, geomModelIdDualCellVarID, geomModelIdDualCell.data()));
  NC_CHECK(nc_put_var_int(ncid, geomModelDimDualCellVarID, geomModelDimDualCell.data()));

  // jigsaw_to_netcdf sets density to all ones, so will we
  std::vector<double> meshDensity(numDualCells, 1);
  NC_CHECK(nc_put_var_double(ncid, meshDensityVarID, meshDensity.data()));

  // Allocate arrays for cell data
  std::vector<double> xDualVertex(numDualVertices);
  std::vector<double> yDualVertex(numDualVertices);
  std::vector<double> zDualVertex(numDualVertices, 0.0);  // Set to 0 for planar mesh
  std::vector<int> geomModelIdDualVertex(numDualVertices);
  std::vector<int> geomModelDimDualVertex(numDualVertices);
  std::vector<int> dualCellsOnDualVertex(numDualVertices * dualVertexDegree, 0); // 1-based indexing required by mpas

  FIter faces = M_faceIter(mesh);
  pFace face;
  int cellIdx = 0;
  while ((face = FIter_next(faces))) {
    // Calculate cell center
    double center[3] = {0.0, 0.0, 0.0};
    pPList faceVerts = F_vertices(face, 1);
    int numFaceVerts = PList_size(faceVerts);
    assert(numFaceVerts == dualVertexDegree);

    void* iter = 0;
    pVertex fv;
    int downVtxIndex = 0;
    while ((fv = (pVertex)PList_next(faceVerts, &iter))) {
      double xyz[3];
      V_coord(fv, xyz);
      center[0] += xyz[0];
      center[1] += xyz[1];
      const int vtxId = EN_id((pEntity)fv);
      dualCellsOnDualVertex[cellIdx * dualVertexDegree + downVtxIndex++] = vtxId + 1; // 1-based indexing
    }
    center[0] /= numFaceVerts;
    center[1] /= numFaceVerts;

    xDualVertex[cellIdx] = center[0];
    yDualVertex[cellIdx] = center[1];
    zDualVertex[cellIdx] = 0.0;  // 2D mesh, z = 0

    // Get geometric classification
    pGEntity gent = EN_whatIn((pEntity)face);
    geomModelIdDualVertex[cellIdx] = GEN_tag(gent);
    geomModelDimDualVertex[cellIdx] = GEN_type(gent);

    // Store mapping
    PList_delete(faceVerts);
    cellIdx++;
  }
  FIter_delete(faces);

  // Write dual vertex data
  NC_CHECK(nc_put_var_double(ncid, xDualVertexVarID, xDualVertex.data()));
  NC_CHECK(nc_put_var_double(ncid, yDualVertexVarID, yDualVertex.data()));
  NC_CHECK(nc_put_var_double(ncid, zDualVertexVarID, zDualVertex.data()));
  NC_CHECK(nc_put_var_int(ncid, geomModelIdDualVertexVarID, geomModelIdDualVertex.data()));
  NC_CHECK(nc_put_var_int(ncid, geomModelDimDualVertexVarID, geomModelDimDualVertex.data()));
  NC_CHECK(nc_put_var_int(ncid, dualCellsOnDualVertexVarID, dualCellsOnDualVertex.data()));

  // Close the NetCDF file
  NC_CHECK(nc_close(ncid));

  std::cout << "Successfully wrote NetCDF file: " << outputFileName << std::endl;
  std::cout << "  Number of vertices: " << numDualVertices << std::endl;
  std::cout << "  Number of cells: " << numDualCells << std::endl;

  return 0;
}
