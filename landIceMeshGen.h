#ifndef LANDICEMESHGEN_H
#define LANDICEMESHGEN_H

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
#include <functional> //std::function

#include "splineInterpolation.h"

struct GeomInfo {
  int numVtx;
  int numEdges;
  std::vector<double> vtx_x;
  std::vector<double> vtx_y;
  std::vector<std::array<int, 2>> edges;
};

struct ModelTopo {
  std::vector<pGVertex> vertices;
  std::vector<pGEdge> edges;
  std::vector<pGFace> faces;
  std::vector<pGEdge> faceEdges; // the array of edges connected to the face
  std::vector<int> faceDirs; // the direction of the edge with respect to the face
  pGRegion region;
  pGIPart part;
  pGModel model;

};


GeomInfo readVtkGeom(std::string fname, bool debug = false);
GeomInfo readJigGeom(std::string fname, bool debug = false);

void convertMetersToKm(GeomInfo &geom);
GeomInfo cleanGeom(GeomInfo &dirty, double coincidentVtxToleranceSquared,
                      bool debug = false);

void createBoundingBoxGeom(ModelTopo& mdlTopo, GeomInfo& geom, bool debug=false);


std::tuple<std::vector<int>,std::vector<int>>
discoverTopology(GeomInfo& geom, double angleTol, double onCurveAngleTol, bool debug = false);
void createEdges(ModelTopo& mdlTopo, GeomInfo& geom, std::vector<int>& isPtOnCurve, std::vector<int>& isMdlVtx, bool debug=false);
void createFaces(ModelTopo& mdlTopo, GeomInfo& geom);
void printModelInfo(pGModel model);
void createMesh(ModelTopo mdlTopo, std::string& meshFileName, pProgress progress);
#endif
