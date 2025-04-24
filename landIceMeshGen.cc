#include "MeshSim.h"
#include "SimAdvModel.h"
#include "SimDisplay.h"
#include "SimInfo.h"
#include "SimModel.h"
#include "SimUtil.h"
#include <iostream>
#include <math.h>

#include <algorithm> //[min|max]element
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits> //std::numeric_limits
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "splineInterpolation.h"

void messageHandler(int type, const char *msg);

struct PlaneBounds {
  double minX;
  double maxX;
  double minY;
  double maxY;
};

struct GeomInfo {
  int numVtx;
  int numEdges;
  std::vector<double> vtx_x;
  std::vector<double> vtx_y;
  std::vector<std::array<int, 2>> edges;
};

std::tuple<std::string, int> readKeyValue(std::ifstream &in,
                                          bool debug = true) {
  std::string key, value;
  std::getline(in, key, '=');
  std::getline(in, value);
  if (debug)
    std::cout << "line: " << key << " " << value << std::endl;
  const int val = std::stoi(value);
  return {key, val};
}

void skipLine(std::ifstream &in, bool debug = true) {
  std::string line;
  std::getline(in, line);
  if (debug)
    std::cout << "skip line: " << line << std::endl;
}

PlaneBounds getBoundingPlane(GeomInfo &geom) {
  auto minX = std::min_element(geom.vtx_x.begin(), geom.vtx_x.end());
  auto maxX = std::max_element(geom.vtx_x.begin(), geom.vtx_x.end());
  auto minY = std::min_element(geom.vtx_y.begin(), geom.vtx_y.end());
  auto maxY = std::max_element(geom.vtx_y.begin(), geom.vtx_y.end());
  return {*minX, *maxX, *minY, *maxY};
}

std::array<double, 3> readPoint(std::ifstream &in, bool debug = true) {
  std::array<double, 3> pt;
  std::string value;
  std::getline(in, value, ';');
  pt[0] = std::stoi(value);
  std::getline(in, value, ';');
  pt[1] = std::stoi(value);
  std::getline(in, value);
  pt[2] = std::stoi(value);
  return pt;
}

std::array<int, 2> readEdge(std::ifstream &in, bool debug = true) {
  std::array<int, 2> edge;
  std::string value;
  std::getline(in, value, ';');
  edge[0] = std::stoi(value);
  std::getline(in, value, ';');
  edge[1] = std::stoi(value);
  // not using the id of the edge
  std::getline(in, value);
  return edge;
}

/*
Removed getline in the following functions to use for VTK files because of
better whitespace handling
*/
std::array<double, 3> readPointVtk(std::ifstream &in, bool debug = true) {
  std::array<double, 3> pt;
  in >> pt[0] >> pt[1] >> pt[2];
  return pt;
}

std::array<int, 2> readEdgeVtk(std::ifstream &in, bool debug = true) {
  int numPoints;
  in >> numPoints;
  assert(numPoints == 2);

  std::array<int, 2> edge;
  in >> edge[0] >> edge[1];

  return edge;
}
GeomInfo readVtkGeom(std::string fname, bool debug = false) {
  std::ifstream vtkFile(fname);
  if (!vtkFile.is_open()) {
    fprintf(stderr, "failed to open VTK geom file %s\n", fname.c_str());
    exit(EXIT_FAILURE);
  }

  GeomInfo geom;

  // Version
  skipLine(vtkFile, debug);
  // Title
  skipLine(vtkFile, debug);
  // Format of VTK
  std::string format;
  vtkFile >> format;
  // Check for ASCII for now
  assert(format == "ASCII");
  // Skip to next line
  skipLine(vtkFile, debug);

  // Dataset Type
  std::string keyword, datasetType;
  vtkFile >> keyword >> datasetType;
  assert(keyword == "DATASET");
  // Check for Polydata for now
  assert(datasetType == "POLYDATA");

  // Read points
  int numPoints;
  std::string dataType;
  vtkFile >> keyword >> numPoints >> dataType;
  assert(keyword == "POINTS");
  geom.numVtx = numPoints;

  // DID not change
  geom.vtx_x.reserve(geom.numVtx);
  geom.vtx_y.reserve(geom.numVtx);

  // point coordinates
  for (int i = 0; i < geom.numVtx; i++) {
    auto pt = readPointVtk(vtkFile, debug);
    geom.vtx_x.push_back(pt[0]);
    geom.vtx_y.push_back(pt[1]);
    if (debug)
      std::cout << "pt " << geom.vtx_x[i] << ", " << geom.vtx_y[i] << std::endl;
  }

  // Read lines
  vtkFile >> keyword >> geom.numEdges;
  int totalIndexCount;
  vtkFile >> totalIndexCount;
  assert(keyword == "LINES");

  geom.edges.reserve(geom.numEdges);

  // edge indices
  for (int i = 0; i < geom.numEdges; i++) {
    geom.edges.push_back(readEdgeVtk(vtkFile, debug));
    if (debug)
      std::cout << "edge " << geom.edges[i][0] << ", " << geom.edges[i][1]
                << std::endl;
  }

  return geom;
}

GeomInfo readJigGeom(std::string fname, bool debug = false) {
  std::ifstream mshFile(fname);
  if (!mshFile.is_open()) {
    fprintf(stderr, "failed to open jigsaw geom file %s\n", fname.c_str());
    exit(EXIT_FAILURE);
  }

  GeomInfo geom;

  // header - skip
  skipLine(mshFile, debug);
  // MSHID - skip
  skipLine(mshFile, debug);
  // NDIMS
  {
    auto [key, value] = readKeyValue(mshFile, debug);
    if (debug)
      std::cout << "key: " << key << " val: " << value << std::endl;
    assert(value == 2);
  }
  // POINT
  {
    auto [key, value] = readKeyValue(mshFile, debug);
    if (debug)
      std::cout << "key: " << key << " val: " << value << std::endl;
    geom.numVtx = value;
  }
  geom.vtx_x.reserve(geom.numVtx);
  geom.vtx_y.reserve(geom.numVtx);
  // point coordinates
  for (int i = 0; i < geom.numVtx; i++) {
    auto pt = readPoint(mshFile, debug);
    geom.vtx_x.push_back(pt[0]);
    geom.vtx_y.push_back(pt[1]);
    if (debug)
      std::cout << "pt " << geom.vtx_x[i] << ", " << geom.vtx_y[i] << std::endl;
  }
  // EDGE
  {
    auto [key, value] = readKeyValue(mshFile, debug);
    if (debug)
      std::cout << "key: " << key << " val: " << value << std::endl;
    geom.numEdges = value;
  }
  geom.edges.reserve(geom.numEdges);
  // edge indices
  for (int i = 0; i < geom.numEdges; i++) {
    geom.edges.push_back(readEdge(mshFile, debug));
    if (debug)
      std::cout << "edge " << geom.edges[i][0] << ", " << geom.edges[i][1]
                << std::endl;
  }

  return geom;
}

bool isPtCoincident(double ax, double ay, double bx, double by,
                    double tolSquared = 1) {
  double xDelta = std::abs(ax - bx);
  double yDelta = std::abs(ay - by);
  double length = xDelta * xDelta + yDelta * yDelta;
  return (length < tolSquared);
}

bool checkVertexUse(GeomInfo &geom, bool debug = false) {
  std::map<int, int> vtxCounter;
  for (int i = 0; i < geom.numVtx; i++)
    vtxCounter[i] = 0;
  for (auto e : geom.edges) {
    vtxCounter[e[0]]++;
    vtxCounter[e[1]]++;
  }
  bool isOk = true;
  for (auto p : vtxCounter) {
    if (p.second != 2) {
      std::cout << "vtx " << p.first << " uses " << p.second << "\n";
      isOk = false;
    }
  }
  return isOk;
}

GeomInfo cleanJigGeom(GeomInfo &dirty, double coincidentVtxToleranceSquared,
                      bool debug = false) {
  assert(checkVertexUse(dirty));
  // trying to check the the dirty geom has a chain of edges
  assert(dirty.numEdges == dirty.numVtx);
  // the first four edges form a loop
  assert(dirty.edges[0][0] == dirty.edges[3][1]);
  // the remaining edges form a loop
  assert(dirty.edges[4][0] == dirty.edges[dirty.numEdges - 1][1]);

  GeomInfo clean;
  clean.vtx_x.reserve(dirty.numVtx);
  clean.vtx_y.reserve(dirty.numVtx);
  // Look for vertices that are nearly coincident
  clean.vtx_x.push_back(dirty.vtx_x[0]);
  clean.vtx_y.push_back(dirty.vtx_y[0]);
  for (int i = 1; i < dirty.numVtx; i++) {
    auto close =
        isPtCoincident(dirty.vtx_x[i - 1], dirty.vtx_y[i - 1], dirty.vtx_x[i],
                       dirty.vtx_y[i], coincidentVtxToleranceSquared);
    if (!close) {
      clean.vtx_x.push_back(dirty.vtx_x[i]);
      clean.vtx_y.push_back(dirty.vtx_y[i]);
    } else {
      if (debug) {
        std::cout << "coincident pt " << i - 1 << " (" << dirty.vtx_x[i - 1]
                  << ", " << dirty.vtx_y[i - 1] << ") " << i << " ("
                  << dirty.vtx_x[i] << ", " << dirty.vtx_y[i] << ")\n";
      }
    }
  }
  clean.numVtx = clean.vtx_x.size();

  // loops have an equal number of verts and edges
  assert(dirty.numVtx >= 4); // there must be a bounding box
  clean.edges.reserve(dirty.numVtx);
  // the first loop is the rectangular boundary
  for (int i = 0; i < 3; i++)
    clean.edges.push_back({i, i + 1});
  clean.edges.push_back({3, 0}); // close the loop
  if (dirty.numVtx > 4) {        // there is another loop
    // the second loop is the grounding line
    for (int i = 4; i < clean.numVtx - 1; i++)
      clean.edges.push_back({i, i + 1});
    clean.edges.push_back({clean.numVtx - 1, 4}); // close the loop
  }

  clean.numEdges = clean.edges.size();
  assert(clean.numEdges == clean.numVtx);
  return clean;
}

std::array<double, 3> subtractPts(double a[3], double b[3]) {
  return {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
}

std::array<double, 3> getNormal(pGEdge first, pGEdge second) {
  // the tail of edge first is the head of edge second
  assert(GE_vertex(first, 1) == GE_vertex(second, 0));
  pGVertex src = GE_vertex(first, 1);
  pGVertex uDest = GE_vertex(first, 0);
  pGVertex vDest = GE_vertex(second, 1);
  double srcPt[3];
  GV_point(src, srcPt);
  double uDestPt[3];
  GV_point(uDest, uDestPt);
  double vDestPt[3];
  GV_point(vDest, vDestPt);
  auto u = subtractPts(uDestPt, srcPt);
  auto v = subtractPts(vDestPt, srcPt);
  return {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
          u[0] * v[1] - u[1] * v[0]};
}

std::string getFileExtension(const std::string &filename) {
  size_t dotPos = filename.rfind('.');
  if (dotPos != std::string::npos) {
    return filename.substr(dotPos);
  }
  return "";
}

pGEdge fitCurveToContour(pGRegion region, pGVertex first, pGVertex last, const int numPts, std::vector<double>& pts) {
  auto curve = SCurve_createPiecewiseLinear(numPts, pts.data());
  pGEdge edge = GR_createEdge(region, first, last, curve, 1);
  return edge;
}

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: <jigsaw .msh or .vtk file> <output prefix> "
                 "<coincidentVtxToleranceSquared> <stride>\n";
    std::cerr << "coincidentVtxToleranceSquared is the mininum allowed "
                 "distance (squared) between adjacent vertices in the "
                 "input.  Vertices within the specified distance will "
                 "be merged.\n";
    std::cerr << "stride defines the maximum number of points to be "
                 "used to define a single geometric model edge.\n"
                 "  A stride less than 1 is invalid.\n"
                 "  A stride of 1 will create one edge per point "
                 "where edge i will be defined by points i-1 "
                 "and i.\n"
                 "  A stride of N will create one edge for sequences of "
                 "N points.  It the number of points is not evenly "
                 "divisible by N then the last edge will contain the "
                 "remaining points.\n";
    return 1;
  }
  assert(argc == 5);

  GeomInfo dirty;
  pGVertex *vertices; // array to store the returned model vertices
  pGFace *faces;      // array to store the returned model faces
  pGRegion region;    // pointer to returned model region
  pGIPart part;
  pGModel model; // pointer to the complete model

  std::string filename = argv[1];
  std::string ext = getFileExtension(filename);

  if (ext == ".vtk") {
    dirty = readVtkGeom(argv[1]);
  } else if (ext == ".msh") {
    dirty = readJigGeom(argv[1]);
  } else {
    std::cerr << "Unsupported file extension: " << ext << "\n";
    return 1;
  }
  auto geom = cleanJigGeom(dirty, std::stof(argv[3]), true);
  const auto prefix = std::string(argv[2]);
  std::string modelFileName = prefix + ".smd";
  std::string meshFileName = prefix + ".sms";

  const auto debug = false;

  // You will want to place a try/catch around all SimModSuite calls,
  // as errors are thrown.
  try {
    Sim_logOn("simMeshGen.log");
    SimModel_start(); // Call before Sim_readLicenseFile
    // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
    // pass in the location of a file containing your keys.  For a release
    // product, use Sim_registerKey()
    Sim_readLicenseFile(0);
    // Tessellation of GeomSim geometry requires Meshing to have started
    MS_init();

    Sim_setMessageHandler(messageHandler);
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    model = GM_new(1);
    part = GM_rootPart(model);
    region = GIP_outerRegion(part);

    vertices = new pGVertex[4];

    // TODO generalize face creation
    if (geom.numEdges > 4) {
      faces = new pGFace[2];
    } else {
      faces = new pGFace[1];
    }

    // First we'll add the vertices
    int i;
    for (i = 0; i < 4; i++) {
      double vtx[3] = {geom.vtx_x[i], geom.vtx_y[i], 0};
      vertices[i] = GR_createVertex(region, vtx);
      if (debug)
        std::cout << "vtx " << i << " (" << vtx[0] << " , " << vtx[1] << ")\n";
    }

    std::vector<pGEdge> edges;
    // Now we'll add the edges
    double point0[3], point1[3]; // xyz locations of the two vertices
    pCurve linearCurve;
    for (i = 0; i < 4; i++) {
      const auto startVertIdx = geom.edges[i][0];
      const auto endVertIdx = geom.edges[i][1];
      auto startVert = vertices[startVertIdx];
      auto endVert = vertices[endVertIdx];
      GV_point(startVert, point0);
      GV_point(endVert, point1);
      linearCurve = SCurve_createLine(point0, point1);
      auto edge = GR_createEdge(region, startVert, endVert, linearCurve, 1);
      edges.push_back(edge);
      if (debug) {
        std::cout << "edge " << i << " (" << point0[0] << " , " << point0[1]
                  << ")"
                  << ",(" << point1[0] << " , " << point1[1] << ")\n";
      }
    }

    const int stride = std::stoi(argv[4]);
    assert(stride > 0);
    const int firstPt = 4;
    std::cout << "numPts " << geom.numVtx-4 << " lastPt " << geom.numVtx << "\n";
    double pt[3] = {geom.vtx_x[firstPt], geom.vtx_y[firstPt], 0};
    std::cout << "creatingVtx " << pt[0] << " " << pt[1] << "\n";
    pGVertex firstVtx = GR_createVertex(region, pt);
    pGVertex prevVtx = firstVtx;
    int prevVtxIdx = firstPt;
    for(i=5; i<=geom.numVtx; i++) {
      if(i%stride == 0 || i == geom.numVtx) {
        const int isLastPt = (i == geom.numVtx ? 1 : 0);
        const int numPts = i - prevVtxIdx + 1;
        std::vector<double> pts(numPts*3);
        int idx = 0;
        for(int j=prevVtxIdx; j<=i && j<geom.numVtx; j++) {
          pts[idx++] = geom.vtx_x[j];
          pts[idx++] = geom.vtx_y[j];
          pts[idx++] = 0;
        }
        pGVertex vtx;
        if(isLastPt) {
          pts[idx++] = geom.vtx_x[firstPt];
          pts[idx++] = geom.vtx_y[firstPt];
          pts[idx++] = 0;
          vtx = firstVtx;
        } else {
          double pt[3] = {geom.vtx_x[i], geom.vtx_y[i], 0};
          std::cout << "creatingVtx " << pt[0] << " " << pt[1] << "\n";
          vtx = GR_createVertex(region, pt);
        }
        std::cout << "range " << prevVtxIdx << " " << i << " numPts " << numPts << " isLastPt " << isLastPt << "\n";
        { double first[3];
          GV_point(prevVtx, first);
          double last[3];
          GV_point(vtx, last);
          std::cout << "start " << first[0] << " " << first[1] << "\n";
          std::cout << "end " << last[0] << " " << last[1] << "\n";
          std::cout << "pts: ";
          for(int j=0; j<pts.size(); j++) {
            std::cout << pts.at(j) << " ";
          } std::cout << "\n";
        }
        auto edge = fitCurveToContour(region, prevVtx, vtx, numPts, pts);
        edges.push_back(edge);
        prevVtx = vtx;
        prevVtxIdx = i;
      }
    }

    auto planeBounds = getBoundingPlane(geom);

    // Now add the faces
    double corner[3], xPt[3],
        yPt[3];        // the points defining the surface of the face
    pGEdge *faceEdges; // the array of edges connected to the face
    int *faceDirs;     // the direction of the edge with respect to the face

    // When defining the loop, will always start with the first edge in the
    // faceEdges array
    pSurface planarSurface;

    // **************
    // Create the face between the bounding rectangle and the grounding line
    // (water)
    // **************
    // Define the surface
    corner[0] = planeBounds.minX;
    corner[1] = planeBounds.minY;
    corner[2] = 0;
    xPt[0] = planeBounds.maxX;
    xPt[1] = planeBounds.minY;
    xPt[2] = 0;
    yPt[0] = planeBounds.minX;
    yPt[1] = planeBounds.maxY;
    yPt[2] = 0;

    const int faceDirectionFwd = 1;
    const int faceDirectionRev = 0;
    const int sameNormal = 1;
    const int oppositeNormal = 0;

    // Create the face
    faceEdges = new pGEdge[edges.size()];
    faceDirs = new int[edges.size()];
    // the first four edges define the outer bounding rectangle
    for (i = 0; i < 4; i++) {
      faceDirs[i] = faceDirectionFwd; // clockwise
      faceEdges[i] = edges.at(i);
    }
    if (edges.size() > 4) {
      // the remaining edges define the grounding line
      // TODO generalize loop creation
      int j = edges.size() - 1;
      for (i = 4; i < edges.size(); i++) {
        faceDirs[i] = faceDirectionRev; // counter clockwise
        // all edges are input in counter clockwise order,
        // reverse the order so the face is on the left (simmetrix requirement)
        faceEdges[i] = edges.at(j--);
      }

      int numLoopsOuterFace = 2;
      int loopFirstEdgeIdx[2] = {0, 4};
      planarSurface = SSurface_createPlane(corner, xPt, yPt);
      faces[0] = GR_createFace(region, edges.size(), faceEdges, faceDirs,
                               numLoopsOuterFace, loopFirstEdgeIdx,
                               planarSurface, sameNormal);
      std::cout << "faces[0] area: " << GF_area(faces[0], 0.2) << "\n";
      assert(GF_area(faces[0], 0.2) > 0);
    } else {
      int numLoopsOuterFace = 1;
      int loopFirstEdgeIdx[1] = {0};
      planarSurface = SSurface_createPlane(corner, xPt, yPt);
      faces[0] = GR_createFace(region, edges.size(), faceEdges, faceDirs,
                               numLoopsOuterFace, loopFirstEdgeIdx,
                               planarSurface, sameNormal);
      std::cout << "faces[0] area: " << GF_area(faces[0], 0.2) << "\n";
      assert(GF_area(faces[0], 0.2) > 0);
    }

    if (edges.size() > 4) {
      // **************
      // Create the 'ice' face bounded by the grounding line
      // **************
      planarSurface = SSurface_createPlane(corner, xPt, yPt);
      const int numEdgesInnerFace = edges.size() - 4;
      const int numLoopsInnerFace = 1;
      int loopFirstEdgeIdx[1] = {0};
      int j = 4;
      for (i = 0; i < numEdgesInnerFace; i++) {
        faceDirs[i] = faceDirectionFwd; // clockwise
        faceEdges[i] = edges.at(j++);
      }
      faces[1] = GR_createFace(region, numEdgesInnerFace, faceEdges, faceDirs,
                               numLoopsInnerFace, loopFirstEdgeIdx,
                               planarSurface, sameNormal);
      std::cout << "faces[1] area: " << GF_area(faces[1], 0.2) << "\n";
      assert(GF_area(faces[1], 0.2) > 0);
    }

    auto isValid = GM_isValid(model, 2, NULL);
    if (!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    } else {
      std::cout << "Model is valid.\n";
    }

    std::cout << "Number of vertices in model: " << GM_numVertices(model)
              << std::endl;
    std::cout << "Number of edges in model: " << GM_numEdges(model)
              << std::endl;
    std::cout << "Number of faces in model: " << GM_numFaces(model)
              << std::endl;
    std::cout << "Number of regions in model: " << GM_numRegions(model)
              << std::endl;
    GM_write(model, modelFileName.c_str(), 0, 0);

    // This next section creates a surface mesh from the model.  You can comment
    // out this section if you don't want to mesh
    pMesh mesh = M_new(0, model);
    pACase meshCase = MS_newMeshCase(model);

    pModelItem domain = GM_domain(model);
    // find the smallest size of the geometric model edges
    auto minGEdgeLen = std::numeric_limits<double>::max();
    for (i = 0; i < edges.size(); i++) {
      auto len = GE_length(faceEdges[i]);
      if (len < minGEdgeLen)
        minGEdgeLen = len;
    }
    std::cout << "Min geometric model edge length: " << minGEdgeLen
              << std::endl;
    const auto contourMeshSize = minGEdgeLen * 32;
    const auto globMeshSize = contourMeshSize * 64;
    std::cout << "Contour absolute mesh size target: " << contourMeshSize
              << std::endl;
    std::cout << "Global absolute mesh size target: " << globMeshSize
              << std::endl;
    MS_setMeshSize(meshCase, domain, 1, globMeshSize, NULL);
    for (i = 4; i < edges.size(); i++)
      MS_setMeshSize(meshCase, faceEdges[i], 1, contourMeshSize, NULL);

    {
      GFIter fIter = GM_faceIter(model);
      pGFace modelFace;
      while (modelFace = GFIter_next(fIter)) {
        const double area = GF_area(modelFace, 0.2);
        std::cout << "face area: " << area << "\n";
        assert(area > 0);
      }
      GFIter_delete(fIter);
    }

    pSurfaceMesher surfMesh = SurfaceMesher_new(meshCase, mesh);
    SurfaceMesher_execute(surfMesh, progress);
    SurfaceMesher_delete(surfMesh);
    std::cout << "Number of mesh faces in surface: " << M_numFaces(mesh)
              << std::endl;

    M_write(mesh, meshFileName.c_str(), 0, progress);
    std::cout << "Number of mesh regions in volume: " << M_numRegions(mesh)
              << std::endl;
    MS_deleteMeshCase(meshCase);
    M_release(mesh);
    // end of meshing section
    
    delete[] faceEdges;
    delete[] faceDirs;
    delete[] vertices;
    delete[] faces;
    // cleanup
    GM_release(model);
    Progress_delete(progress);
    MS_exit();
    Sim_unregisterAllKeys();
    SimModel_stop();
    Sim_logOff();

  } catch (pSimInfo err) {
    std::cerr << "SimModSuite error caught:" << std::endl;
    std::cerr << "  Error code: " << SimInfo_code(err) << std::endl;
    std::cerr << "  Error string: " << SimInfo_toString(err) << std::endl;
    SimInfo_delete(err);
    return 1;
  } catch (...) {
    std::cerr << "Unhandled exception caught" << std::endl;
    return 1;
  }
  return 0;
}

void messageHandler(int type, const char *msg) {
  switch (type) {
  case Sim_InfoMsg:
    std::cout << "Info: " << msg << std::endl;
    break;
  case Sim_DebugMsg:
    std::cout << "Debug: " << msg << std::endl;
    break;
  case Sim_WarningMsg:
    std::cout << "Warning: " << msg << std::endl;
    break;
  case Sim_ErrorMsg:
    std::cout << "Error: " << msg << std::endl;
    break;
  }
  return;
}
