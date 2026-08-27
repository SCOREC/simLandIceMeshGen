#include "simModelGen2d.h"
#include "Quadtree.h"
#include <map>

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



double getPt2PtEdgeLength(pGEdge edge) {
  pGVertex start = GE_vertex(edge, 1);
  pGVertex end = GE_vertex(edge, 0);
  double startPt[3];
  GV_point(start, startPt);
  double endPt[3];
  GV_point(end, endPt);
  auto lenSq = getLengthSquared(startPt[0], startPt[1], endPt[0], endPt[1]);
  return std::sqrt(lenSq);
}

pGEdge fitCurveToContourSimInterp(bool isLinearSpline, pGRegion region, pGVertex first, pGVertex last,
                         std::vector<double>& pts, bool debug=false) {
  assert(pts.size() % 3 == 0); //pts must contain coordinates x1,y1,z1, x2,y2,z2, ...
  const int numPts = pts.size()/3;
  assert(numPts > 1);
  pCurve curve;
  if( isLinearSpline || numPts == 2 || numPts == 3) {
    curve = SCurve_createPiecewiseLinear(numPts, &pts[0]); //TODO - replace withe bspline?
  } else {
    const int order = 4;
    curve = SCurve_createInterpolatedBSpline(order, numPts, &pts[0], NULL);
  }
  pGEdge edge = GR_createEdge(region, first, last, curve, 1);
  if(numPts>=4 && debug) {
    const auto p2pLength = getPt2PtEdgeLength(edge);
    const auto eLength = GE_length(edge);
    if( eLength > 1.5*p2pLength ) {
      std::cerr << "Warning: curve length " << eLength << " is more than 1.5 times longer than the end point to end point length " << p2pLength << "\n";
    }
  }
  return edge;
}

void printModelInfo(pGModel model) {
  std::cout << "Number of vertices in model: " << GM_numVertices(model)
    << std::endl;
  std::cout << "Number of edges in model: " << GM_numEdges(model)
    << std::endl;
  std::cout << "Number of faces in model: " << GM_numFaces(model)
    << std::endl;
  std::cout << "Number of regions in model: " << GM_numRegions(model)
    << std::endl;
}

void setClassification(GeomInfo& geom, PointClassification& ptClass, BoundaryClassification& bndClass, const int firstPt, const int numPts, pGVertex startingMdlVtx, pGVertex endMdlVtx, pGEdge edge, const int splineIdx) {
  const auto vtxDim = 0;
  ptClass.id.at(firstPt) = GEN_tag(startingMdlVtx);
  ptClass.dim.at(firstPt) = vtxDim;
  ptClass.splineIdx.at(firstPt) = splineIdx;
  bndClass.handle.at(firstPt) = startingMdlVtx;

  auto ptIdx = geom.getNextPtIdx(firstPt); //handle wrap around in indexing
  if(numPts > 2) {
    const auto edgeTag = GEN_tag(edge);
    const auto edgeDim = 1;
    auto ptCount = 1;
    while(ptCount < numPts-1) {
      ptClass.id.at(ptIdx) = edgeTag;
      ptClass.dim.at(ptIdx) = edgeDim;
      ptClass.splineIdx.at(ptIdx) = splineIdx;
      bndClass.handle.at(ptIdx) = edge;
      ptIdx = geom.getNextPtIdx(ptIdx);
      ptCount++;
    }
  }

  ptClass.id.at(ptIdx) = GEN_tag(endMdlVtx);
  ptClass.dim.at(ptIdx) = vtxDim;
  ptClass.splineIdx.at(ptIdx) = splineIdx;
  bndClass.handle.at(ptIdx) = endMdlVtx;

  //record which model edge spans each segment [firstPt..ptIdx]
  //this handles recording the edge formed by two adjacent points
  //that are model vertices
  auto segPtIdx = firstPt;
  for (int segCount = 0; segCount < numPts-1; segCount++) {
    bndClass.nextSegmentEdge.at(segPtIdx) = edge;
    segPtIdx = geom.getNextPtIdx(segPtIdx);
  }
}

void createEdges(ModelTopo& mdlTopo, GeomInfo& geom, PointClassification& ptClass, BoundaryClassification& bndClass, SplineInterp::SplineInfo& splines, std::vector<int>& isPtOnCurve, std::vector<int>& isMdlVtx, const bool debug) {
  if(geom.numVtx <= 0) { //no contour
    return;
  }

  enum class State {MdlVtx = 0, OnCurve = 1, NotOnCurve = 2};
  enum class Action {Init, Advance, Line, Curve, LinearSpline};
  typedef std::pair<State,Action> psa; // next state, action
  using func=std::function<psa(int pt)>;
  using funcIntBool=std::function<psa(int pt, bool)>;


  pGVertex firstMdlVtx;
  int firstPtIdx;
  int startingCurvePtIdx;
  pGVertex startingMdlVtx;
  std::vector<int> ptsOnCurve;
  funcIntBool createCurve = [&](int pt, bool isLinearSpline=false) {
    assert(ptsOnCurve.size() >= 2);
    double vtx[3] = {geom.vtx_x[pt], geom.vtx_y[pt], 0};
    pGVertex endMdlVtx;
    if(pt == firstPtIdx) { //wrap around
      endMdlVtx = firstMdlVtx;
    } else {
      endMdlVtx = GR_createVertex(mdlTopo.region, vtx);
      mdlTopo.vertices.push_back(endMdlVtx);
    }

    std::vector<double> pts(ptsOnCurve.size()*3);
    for(int i=0, j = 0; j<ptsOnCurve.size(); j++, i+=3) {
      const int ptIdx = ptsOnCurve.at(j);
      pts[i]   = geom.vtx_x[ptIdx];
      pts[i+1] = geom.vtx_y[ptIdx];
      pts[i+2] = 0;
    }

    double first[3];
    GV_point(startingMdlVtx, first);
    const double tol = 1e-12;
    if( ! isPtCoincident(pts[0], pts[1], first[0], first[1], tol)) {
      std::cerr << "first model vtx does not match first point on curve!... exiting\n";
      exit(EXIT_FAILURE);
    }

    double last[3];
    GV_point(endMdlVtx, last);
    const int lpIdx = (ptsOnCurve.size()-1)*3;
    if( ! isPtCoincident(pts[lpIdx], pts[lpIdx+1], last[0], last[1], tol) ) {
      std::cerr << "last model vtx does not match last point on curve!... exiting\n";
      exit(EXIT_FAILURE);
    }

    auto edge = fitCurveToContourSimInterp(isLinearSpline, mdlTopo.region, startingMdlVtx, endMdlVtx, pts, debug);
    setClassification(geom, ptClass, bndClass, startingCurvePtIdx, ptsOnCurve.size(), startingMdlVtx, endMdlVtx, edge, splines.size());
    if(isLinearSpline) {
      splines.addSpline(SplineInterp::attach_piecewise_linear_curve(pts));
    } else {
      splines.addSpline(SplineInterp::fitCubicSplineToPoints(pts));
    }
    mdlTopo.edges.push_back(edge);

    if (debug) {
      std::cerr << "edge " << mdlTopo.edges.size()
        << " splineIdx " << splines.size() << " "
        << " isLinearSpline " << isLinearSpline << " "
        << " range " << startingCurvePtIdx << " " << pt
        << " numPts " << ptsOnCurve.size() << "\n";
      std::cout << "start " << first[0] << " " << first[1] << "\n";
      std::cout << "end " << last[0] << " " << last[1] << "\n";
      std::cout << "x,y,z\n";
      for(int j=0; j<pts.size(); j+=3) {
        std::cerr << pts.at(j) << ", " << pts.at(j+1) << ", " << pts.at(j+2) << "\n";
      }
    }

    startingMdlVtx = endMdlVtx;
    startingCurvePtIdx = pt;
    ptsOnCurve.clear();  //FIXME - find a cheaper way
    ptsOnCurve.push_back(pt);
    return psa{State::MdlVtx,Action::Curve};
  };
  func createLinearSpline = [&](int pt) {
    return createCurve(pt,true);
  };
  func createBSpline = [&](int pt) {
    return createCurve(pt,false);
  };
  func createLine = [&](int pt) {
    assert(ptsOnCurve.size() == 1);
    ptsOnCurve.push_back(pt);
    auto ignored = createCurve(pt, true);
    return psa{State::MdlVtx,Action::Line};
  };
  func advance = [&](int pt) {
    ptsOnCurve.push_back(pt);
    return psa{State::OnCurve,Action::Advance};
  };
  func advanceLinearSpline = [&](int pt) {
    ptsOnCurve.push_back(pt);
    return psa{State::NotOnCurve,Action::Advance};
  };
  func createCurveFromPriorPt = [&](int pt) {
    if(ptsOnCurve.size() == 1 ) {
      return createLine(pt);
    } else if (ptsOnCurve.size() >= 2 ) {
      //we are not adding the current point yet, so there must be
      //at least two points in the list to form a curve
      auto ignored = createBSpline(ptsOnCurve.back());
      ptsOnCurve.push_back(pt);
      return psa{State::NotOnCurve,Action::LinearSpline};
    } else {
      std::cerr << "createCurveFromPriorPt: no points on the curve.... exiting\n";
      exit(EXIT_FAILURE);
    }
  };
  func createLinearSplineFromPriorPt = [&](int pt) {
    if(ptsOnCurve.size() == 1 ) {
      return createLine(pt);
    } else if (ptsOnCurve.size() >= 2 ) {
      //we are not adding the current point yet, so there must be
      //at least two points in the list to form a curve
      auto ignored = createLinearSpline(ptsOnCurve.back());
      ptsOnCurve.push_back(pt);
      return psa{State::OnCurve,Action::LinearSpline};
    } else {
      std::cerr << "createLinearSplineFromPriorPt: no points on the curve.... exiting\n";
      exit(EXIT_FAILURE);
    }
  };
  func createCurveFromCurrentPt = [&](int pt) {
    ptsOnCurve.push_back(pt);
    auto ret = createBSpline(pt);
    return ret;
  };
  func createLinearSplineFromCurrentPt = [&](int pt) {
    ptsOnCurve.push_back(pt);
    auto ret = createLinearSpline(pt);
    return ret;
  };
  func fail = [&](int pt) {
    std::cerr << "bad state.... exiting\n";
    exit(EXIT_FAILURE);
    return psa{State::OnCurve,Action::Advance};
  };
  typedef std::pair<State,State> pss; // current state, next state
  std::map<pss,func> machine = {
    {{State::MdlVtx,State::MdlVtx}, createLine},
    {{State::MdlVtx,State::OnCurve}, advance},
    {{State::MdlVtx,State::NotOnCurve}, advanceLinearSpline},
    {{State::OnCurve,State::MdlVtx}, createCurveFromCurrentPt},
    {{State::OnCurve,State::OnCurve}, advance},
    {{State::OnCurve,State::NotOnCurve}, createCurveFromPriorPt},
    {{State::NotOnCurve,State::MdlVtx}, createLinearSplineFromCurrentPt},
    {{State::NotOnCurve,State::OnCurve}, createLinearSplineFromPriorPt},
    {{State::NotOnCurve,State::NotOnCurve}, advanceLinearSpline}
  };

  const int isVtx=1;
  firstPtIdx = startingCurvePtIdx = findFirstPt(isMdlVtx, geom.firstContourPt, isVtx);
  if(firstPtIdx == -1) {
    std::cerr << "Error: at least one point must be marked as a model vertex... exiting\n";
    assert(false);
    exit(EXIT_FAILURE);
  }
  double vtx[3] = {geom.vtx_x[startingCurvePtIdx], geom.vtx_y[startingCurvePtIdx], 0};
  firstMdlVtx = startingMdlVtx = GR_createVertex(mdlTopo.region, vtx);
  mdlTopo.vertices.push_back(firstMdlVtx);
  ptsOnCurve.push_back(startingCurvePtIdx);

  State state = State::MdlVtx;
  int ptsVisited = 0; //don't count the first vertex until we close the loop
  int ptIdx = startingCurvePtIdx+1;
  while(ptsVisited < isMdlVtx.size()-geom.firstContourPt) {
    State nextState;
    if(isMdlVtx[ptIdx] == 1) {
      nextState = State::MdlVtx;
    } else if (isPtOnCurve[ptIdx] == 1) {
      nextState = State::OnCurve;
    } else if (isPtOnCurve[ptIdx] == 0) {
      nextState = State::NotOnCurve;
    } else {
      exit(EXIT_FAILURE);
    }
    psa res = machine[{state,nextState}](ptIdx);
    state = res.first;
    ptsVisited++;
    ptIdx = geom.getNextPtIdx(ptIdx);
  }

}

void createFace(ModelTopo& mdlTopo, PlaneBounds& planeBounds, bool debug) {
  const double corner[3] = {planeBounds.minX, planeBounds.minY, 0};
  const double xPt[3] = {planeBounds.maxX, planeBounds.minY, 0};
  const double yPt[3] = {planeBounds.minX, planeBounds.maxY, 0};
  const int faceDirectionFwd = 1;
  const int sameNormal = 1;

  pSurface planarSurface = SSurface_createPlane(corner, xPt, yPt);
  const int numEdgesInnerFace = mdlTopo.edges.size();
  const int numLoopsInnerFace = 1;
  int loopFirstEdgeIdx[1] = {0};
  for (int i = 0; i < numEdgesInnerFace; i++) {
    mdlTopo.faceDirs.push_back(faceDirectionFwd); // clockwise
    mdlTopo.faceEdges.push_back(mdlTopo.edges.at(i));
  }
  mdlTopo.faces.push_back(GR_createFace(mdlTopo.region, numEdgesInnerFace,
        mdlTopo.faceEdges.data(),
        mdlTopo.faceDirs.data(),
        numLoopsInnerFace, loopFirstEdgeIdx,
        planarSurface, sameNormal));
  if(debug) {
    std::cout << "faces[1] area: " << GF_area(mdlTopo.faces.at(0), 0.2) << "\n";
  }
  assert(GF_area(mdlTopo.faces.at(0), 0.2) > 0);
}

void createFaces(ModelTopo& mdlTopo, PlaneBounds& planeBounds, bool hasSingleContour, int numOuterEdges, bool debug) {
  // Now add the faces
  double corner[3], xPt[3], yPt[3]; // the points defining the surface of the face

  // When defining the loop, will always start with the first edge in the
  // faceEdges array
  pSurface planarSurface;

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

  if (!hasSingleContour) {
    // **************
    // Create the face between the bounding rectangle and the grounding line
    // (water)
    // **************
    // the first numOuterEdges edges define the outer contour
    for (int i = 0; i < numOuterEdges; i++) {
      mdlTopo.faceDirs.push_back(faceDirectionFwd); // clockwise
      mdlTopo.faceEdges.push_back(mdlTopo.edges.at(i));
    }
    // the remaining edges define the grounding line
    // TODO generalize loop creation
    int j = mdlTopo.edges.size() - 1;
    for (int i = numOuterEdges; i < mdlTopo.edges.size(); i++) {
      mdlTopo.faceDirs.push_back(faceDirectionRev); // counter clockwise
      // all edges are input in counter clockwise order,
      // reverse the order so the face is on the left (simmetrix requirement)
      mdlTopo.faceEdges.push_back(mdlTopo.edges.at(j--));
    }
    int numLoopsOuterFace = 2;
    int loopFirstEdgeIdx[2] = {0, numOuterEdges};
    planarSurface = SSurface_createPlane(corner, xPt, yPt);
    mdlTopo.faces.push_back(GR_createFace(mdlTopo.region, mdlTopo.edges.size(),
          mdlTopo.faceEdges.data(),
          mdlTopo.faceDirs.data(),
          numLoopsOuterFace, loopFirstEdgeIdx,
          planarSurface, sameNormal));
    if(debug) {
      std::cout << "faces[0] area: " << GF_area(mdlTopo.faces[0], 0.2) << "\n";
    }
    assert(GF_area(mdlTopo.faces[0], 0.2) > 0);

    mdlTopo.faceEdges.clear();
    mdlTopo.faceDirs.clear();

    // **************
    // Create the 'ice' face bounded by the grounding line
    // **************
    planarSurface = SSurface_createPlane(corner, xPt, yPt);
    const int numEdgesInnerFace = mdlTopo.edges.size() - numOuterEdges;
    const int numLoopsInnerFace = 1;
    int loopFirstEdgeInnerIdx[1] = {0};
    int k = numOuterEdges;
    for (int i = 0; i < numEdgesInnerFace; i++) {
      mdlTopo.faceDirs.push_back(faceDirectionFwd); // clockwise
      mdlTopo.faceEdges.push_back(mdlTopo.edges.at(k++));
    }
    mdlTopo.faces.push_back(GR_createFace(mdlTopo.region, numEdgesInnerFace,
          mdlTopo.faceEdges.data(),
          mdlTopo.faceDirs.data(),
          numLoopsInnerFace, loopFirstEdgeInnerIdx,
          planarSurface, sameNormal));
    if(debug) {
      std::cout << "faces[1] area: " << GF_area(mdlTopo.faces[1], 0.2) << "\n";
    }
    assert(GF_area(mdlTopo.faces[1], 0.2) > 0);
  } else {
    // **************
    // Single contour: the contour itself is the domain boundary
    // **************
    // the contour is CCW; faceDirectionFwd places the face on the left
    for (int i = 0; i < (int)mdlTopo.edges.size(); i++) {
      mdlTopo.faceDirs.push_back(faceDirectionFwd);
      mdlTopo.faceEdges.push_back(mdlTopo.edges.at(i));
    }
    int numLoopsOuterFace = 1;
    int loopFirstEdgeIdx[1] = {0};
    planarSurface = SSurface_createPlane(corner, xPt, yPt);
    mdlTopo.faces.push_back(GR_createFace(mdlTopo.region, mdlTopo.edges.size(),
          mdlTopo.faceEdges.data(),
          mdlTopo.faceDirs.data(),
          numLoopsOuterFace, loopFirstEdgeIdx,
          planarSurface, sameNormal));
    if(debug) {
      std::cout << "faces[0] area: " << GF_area(mdlTopo.faces[0], 0.2) << "\n";
    }
    assert(GF_area(mdlTopo.faces[0], 0.2) > 0);
  }
}

void specifyBoundaryTriangleMesh(pMesh mesh, GeomInfo& outerGeom, BoundaryClassification& bndClassOuter, pGFace outerFace, bool debug) {
  const auto numAllVtx = (int)outerGeom.all_vertices_x.size();

  //map from boundary traversal order to all_vertices
  std::vector<int> boundaryPosToAllIdx(outerGeom.numVtx, -1);
  for (int i = 0; i < numAllVtx; i++) {
    const auto pos = outerGeom.boundaryOrder.at(i);
    if (pos >= 0) {
      boundaryPosToAllIdx.at(pos) = i;
    }
  }
  for (int k = 0; k < outerGeom.numVtx; k++) {
    if (boundaryPosToAllIdx.at(k) < 0) {
      std::cerr << "ERROR: boundary traversal position " << k
                 << " has no corresponding all_vertices point (boundaryOrder "
                    "is missing an entry)... exiting\n";
      exit(EXIT_FAILURE);
    }
  }

  //specify a mesh vertex for every point (boundary and interior), tagged by
  //its all_vertices index so triangles can reference vertices by that same
  //index when specifying mesh faces below.
  for (int i = 0; i < numAllVtx; i++) {
    double xyz[3] = {outerGeom.all_vertices_x[i], outerGeom.all_vertices_y[i], 0};
    const auto pos = outerGeom.boundaryOrder.at(i);
    pGEntity ent = (pos >= 0) ? bndClassOuter.handle.at(pos) : (pGEntity)outerFace;
    if (ent == nullptr) {
      std::cerr << "ERROR: all_vertices point " << i << " (boundary position "
                 << pos << ") has no model entity to classify against... "
                    "exiting\n";
      exit(EXIT_FAILURE);
    }
    if (debug && pos >= 0 && GEN_type(ent) == Gedge) {
      double closest[3], param;
      GE_closestPoint((pGEdge)ent, xyz, closest, &param);
      const auto dx = xyz[0]-closest[0], dy = xyz[1]-closest[1];
      const auto dist = std::sqrt(dx*dx+dy*dy);
      const auto tol = GEN_tolerance(ent);
      if (dist > tol) {
        std::cerr << "WARNING: all_vertices point " << i
                   << " (boundary position " << pos << ") is " << dist
                   << " from its classified model edge, tolerance is "
                   << tol << "\n";
      }
    }
    MS_specifyVertex(mesh, xyz, NULL, ent, i);
  }

  //create a mesh edge for every boundary segment for explicit 
  //classification on the model edges
  for (int k = 0; k < outerGeom.numVtx; k++) {
    const auto kNext = (k + 1) % outerGeom.numVtx;
    pGEdge edgeEnt = bndClassOuter.nextSegmentEdge.at(k);
    if (edgeEnt == nullptr) {
      std::cerr << "ERROR: boundary segment " << k << "->" << kNext
                 << " has no model edge to classify against... exiting\n";
      exit(EXIT_FAILURE);
    }
    int vertTags[2] = {boundaryPosToAllIdx.at(k), boundaryPosToAllIdx.at(kNext)};
    MS_specifyEdge(mesh, vertTags, edgeEnt, -1);
  }

  //specify a mesh face for every input triangle
  for (auto& tri : outerGeom.triangles) {
    for (int v : tri) {
      if (v < 0 || v >= numAllVtx) {
        std::cerr << "ERROR: triangle references out-of-range vertex index "
                   << v << " (numAllVtx=" << numAllVtx << ")... exiting\n";
        exit(EXIT_FAILURE);
      }
    }
    int vertTags[3] = {tri[0], tri[1], tri[2]};
    MS_specifyFace(mesh, 3, vertTags, outerFace, -1);
  }

  if (debug) {
    std::cerr << "specified " << numAllVtx << " boundary-triangle mesh "
                 "vertices and " << outerGeom.triangles.size() << " mesh faces\n";
  }
}

pMesh createMesh(ModelTopo mdlTopo, GeomInfo& outerGeom, BoundaryClassification& bndClassOuter, pGFace outerFace, std::string& meshFileName, pProgress progress, bool debug) {
  pMesh mesh = M_new(0, mdlTopo.model);
  if (outerGeom.hasBoundaryTriangles()) {
    specifyBoundaryTriangleMesh(mesh, outerGeom, bndClassOuter, outerFace, debug);
  }
  pACase meshCase = MS_newMeshCase(mdlTopo.model);

  pModelItem domain = GM_domain(mdlTopo.model);
  const auto relativeSizeType = 2;
  MS_setMeshSize(meshCase, domain, relativeSizeType, 0.1, NULL);
  const auto calcCurvatureFromEdgesAndFaces = 3;
  MS_setMeshCurv(meshCase, domain, relativeSizeType, 0.001, calcCurvatureFromEdgesAndFaces);
  MS_setMinCurvSize(meshCase, domain, relativeSizeType, 0.0005);
  //make the transition from the fine mesh to the coarse mesh slower, default
  //rate ~= 2/3
  const auto slowGradationRate = 0.3;
  MS_setGlobalSizeGradationRate(meshCase, slowGradationRate);

  pSurfaceMesher surfMesh = SurfaceMesher_new(meshCase, mesh);
  SurfaceMesher_execute(surfMesh, progress);
  SurfaceMesher_delete(surfMesh);
  std::cout << "Number of mesh faces in surface: " << M_numFaces(mesh)
    << std::endl;

  M_write(mesh, meshFileName.c_str(), 0, progress);
  std::cout << "Number of mesh regions in volume: " << M_numRegions(mesh)
    << std::endl;
  MS_deleteMeshCase(meshCase);
  return mesh;
}

