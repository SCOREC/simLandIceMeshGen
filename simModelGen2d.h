#ifndef SIMMODELGEN2D_H
#define SIMMODELGEN2D_H

#include "MeshSim.h"
#include "SimAdvModel.h"
#include "SimDisplay.h"
#include "SimInfo.h"
#include "SimModel.h"
#include "SimUtil.h"
#include <iostream>
#include <math.h>

#include <string>
#include <vector>

#include "modelGen2d.h"

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



void createBoundingBoxGeom(ModelTopo& mdlTopo, GeomInfo& geom, SplineInterp::SplineInfo& splines, bool debug=false);


void createEdges(ModelTopo& mdlTopo, GeomInfo& geom, PointClassification& ptClass, SplineInterp::SplineInfo& splines, std::vector<int>& isPtOnCurve, std::vector<int>& isMdlVtx, const bool debug=false);
void createFaces(ModelTopo& mdlTopo, GeomInfo& geom, bool debug=false);
void printModelInfo(pGModel model);
void createMesh(ModelTopo mdlTopo, std::string& meshFileName, pProgress progress, bool debug=false);


#endif
