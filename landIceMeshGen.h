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

//FIXME - make this a class
struct GeomInfo {
  int numVtx;
  int numEdges;
  std::vector<double> vtx_x;
  std::vector<double> vtx_y;
  std::vector<std::array<int, 2>> edges;
  int firstContourPt;
  int getPrevPtIdx(int i) {
    if(i == firstContourPt) {
      return vtx_x.size()-1;
    } else {
      return i-1;
    }
  }
  int getNextPtIdx(int i) {
    if(i == vtx_x.size()-1) {
      return firstContourPt;
    } else {
      return i+1;
    }
  }
  void reverseContourPoints() {
    std::reverse(vtx_x.begin()+firstContourPt, vtx_x.end());
    std::reverse(vtx_y.begin()+firstContourPt, vtx_y.end());
    edges.clear(); //indices become invalid
  }
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
discoverTopology(GeomInfo& geom, double coincidentPtTolSquared, double angleTol, double onCurveAngleTol, bool debug = false);
void createEdges(ModelTopo& mdlTopo, GeomInfo& geom, std::vector<int>& isPtOnCurve, std::vector<int>& isMdlVtx, const bool debug=false, const int firstContourPt=4);
void createFaces(ModelTopo& mdlTopo, GeomInfo& geom);
void printModelInfo(pGModel model);
void createMesh(ModelTopo mdlTopo, std::string& meshFileName, pProgress progress);

class OnCurve {
  public:
  OnCurve(double onCurveAngleTol);
  //similar to scorec/tomms @ 2f97d13 (simapis-mod branch)
  int operator()(double tc_m1, double tc, double tc_p1);
  double getLowerTolTC() const { return tc_angle_lower; }
  double getUpperTolTC() const { return tc_angle_upper; }
  private:
  const double deg_angle_lower;
  const double deg_angle_upper;
  const double tc_angle_lower;
  const double tc_angle_upper;
};

#endif
