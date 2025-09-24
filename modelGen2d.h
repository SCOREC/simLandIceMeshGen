#ifndef MODELGEN2D_H
#define MODELGEN2D_H

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
#include <math.h>

#include "splineInterpolation.h"

struct PointClassification {
  std::vector<int> dim;
  std::vector<int> id;
  std::vector<int> splineIdx;
  void writeToOsh(std::string filename);
  PointClassification(int n) : dim(n), id(n), splineIdx(n) {}
};

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

struct PlaneBounds {
  double minX;
  double maxX;
  double minY;
  double maxY;
};

PlaneBounds getBoundingPlane(GeomInfo &geom);


GeomInfo readVtkGeom(std::string fname, bool debug = false);
GeomInfo readJigGeom(std::string fname, bool debug = false);

double getLengthSquared(double ax, double ay, double bx, double by);
bool isPtCoincident(double ax, double ay, double bx, double by,
                    double tolSquared = 1);

int findFirstPt(std::vector<int>& prop, const int offset, const int match);

void convertMetersToKm(GeomInfo &geom);
GeomInfo cleanGeom(GeomInfo &dirty, double coincidentVtxToleranceSquared,
                      bool debug = false);

std::tuple<std::vector<int>,std::vector<int>>
discoverTopology(GeomInfo& geom, double coincidentPtTolSquared, double angleTol, double onCurveAngleTol, bool debug = false);

class OnCurve {
  public:
  OnCurve(double onCurveAngleTol, bool isDebug=false);
  //similar to scorec/tomms @ 2f97d13 (simapis-mod branch)
  int operator()(double tc_m1, double tc, double tc_p1);
  double getLowerTolTC() const { return tc_angle_lower; }
  double getUpperTolTC() const { return tc_angle_upper; }
  private:
  const double deg_angle_lower;
  const double deg_angle_upper;
  const double tc_angle_lower;
  const double tc_angle_upper;
  const bool debug;
};

#endif
