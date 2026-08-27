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

// per-outer-boundary-point (indexed like geom.vtx_x/y) handle of the model
// vertex/edge it is classified on, captured during createEdges.
struct BoundaryClassification {
  std::vector<pGEntity> handle;
  // nextSegmentEdge[k] is the model edge spanning boundary points k and
  // k+1 (wrapping at numVtx-1 -> 0). A model edge can be a
  // single segment between two adjacent model vertices, in which case
  // handle[k] and handle[k+1] are both vertices with no edge between them
  // to infer from.
  std::vector<pGEdge> nextSegmentEdge;
  BoundaryClassification(int n) : handle(n), nextSegmentEdge(n) {}
};

void createEdges(ModelTopo& mdlTopo, GeomInfo& geom, PointClassification& ptClass, BoundaryClassification& bndClass, SplineInterp::SplineInfo& splines, std::vector<int>& isPtOnCurve, std::vector<int>& isMdlVtx, const bool debug=false);
void createFace(ModelTopo& mdlTopo, PlaneBounds& planeBounds, bool debug=false);
void createFaces(ModelTopo& mdlTopo, PlaneBounds& planeBounds, bool hasSingleContour, int numOuterEdges, bool debug=false);
void printModelInfo(pGModel model);
void specifyBoundaryTriangleMesh(pMesh mesh, GeomInfo& outerGeom, BoundaryClassification& bndClassOuter, pGFace outerFace, bool debug=false);
pMesh createMesh(ModelTopo mdlTopo, GeomInfo& outerGeom, BoundaryClassification& bndClassOuter, pGFace outerFace, std::string& meshFileName, pProgress progress, bool debug=false);

#endif
