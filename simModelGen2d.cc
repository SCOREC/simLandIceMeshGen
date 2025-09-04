#include "simModelGen2d.h"
#include "Quadtree.h"
#include <map>
#include <numeric> //accumulate

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

void createEdges(ModelTopo& mdlTopo, GeomInfo& geom, std::vector<int>& isPtOnCurve, std::vector<int>& isMdlVtx, const bool debug) {
  if(geom.numVtx <= geom.firstContourPt) { // no internal contour
    return;
  }

  enum class State {MdlVtx = 0, OnCurve = 1, NotOnCurve = 2};
  enum class Action {Init, Advance, Line, Curve, LinearSpline};
  typedef std::pair<State,Action> psa; // next state, action
  using func=std::function<psa(int pt)>;
  using funcIntBool=std::function<psa(int pt, bool)>;

  auto numMdlVerts = std::accumulate(isMdlVtx.begin()+geom.firstContourPt, isMdlVtx.end(), 0);
  std::vector<SplineInterp::BSpline2d> splines;
  splines.reserve(numMdlVerts+1);

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
    if(isLinearSpline) {
      splines.emplace_back(SplineInterp::attach_piecewise_linear_curve(pts));
    } else {
      splines.emplace_back(SplineInterp::fitCubicSplineToPoints(pts));
    }
    mdlTopo.edges.push_back(edge);

    if (debug) {
      std::cerr << "edge " << mdlTopo.edges.size()
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
    auto ignored = createBSpline(pt);
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

void createBoundingBoxGeom(ModelTopo& mdlTopo, GeomInfo& geom, bool debug) {
  // TODO generalize face creation
  if (geom.numEdges > 4) {
    mdlTopo.faces.reserve(2);
  } else {
    mdlTopo.faces.reserve(1);
  }

  // First we'll add the vertices
  int i;
  for (i = 0; i < 4; i++) {
    double vtx[3] = {geom.vtx_x[i], geom.vtx_y[i], 0};
    mdlTopo.vertices.push_back(GR_createVertex(mdlTopo.region, vtx));
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
    auto startVert = mdlTopo.vertices.at(startVertIdx);
    auto endVert = mdlTopo.vertices.at(endVertIdx);
    GV_point(startVert, point0);
    GV_point(endVert, point1);
    linearCurve = SCurve_createLine(point0, point1);
    auto edge = GR_createEdge(mdlTopo.region, startVert, endVert, linearCurve, 1);
    mdlTopo.edges.push_back(edge);
    if (debug) {
      std::cout << "edge " << i << " (" << point0[0] << " , " << point0[1]
        << ")"
        << ",(" << point1[0] << " , " << point1[1] << ")\n";
    }
  }

}

void createFaces(ModelTopo& mdlTopo, GeomInfo& geom, bool debug) {
  auto planeBounds = getBoundingPlane(geom);
  // Now add the faces
  double corner[3], xPt[3], yPt[3]; // the points defining the surface of the face

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
  // the first four edges define the outer bounding rectangle
  for (int i = 0; i < 4; i++) {
    mdlTopo.faceDirs.push_back(faceDirectionFwd); // clockwise
    mdlTopo.faceEdges.push_back(mdlTopo.edges.at(i));
  }
  if (mdlTopo.edges.size() > 4) {
    // the remaining edges define the grounding line
    // TODO generalize loop creation
    int j = mdlTopo.edges.size() - 1;
    for (int i = 4; i < mdlTopo.edges.size(); i++) {
      mdlTopo.faceDirs.push_back(faceDirectionRev); // counter clockwise
      // all edges are input in counter clockwise order,
      // reverse the order so the face is on the left (simmetrix requirement)
      mdlTopo.faceEdges.push_back(mdlTopo.edges.at(j--));
    }

    int numLoopsOuterFace = 2;
    int loopFirstEdgeIdx[2] = {0, 4};
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
  } else {
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

  mdlTopo.faceEdges.clear();
  mdlTopo.faceDirs.clear();

  if (mdlTopo.edges.size() > 4) {
    // **************
    // Create the 'ice' face bounded by the grounding line
    // **************
    planarSurface = SSurface_createPlane(corner, xPt, yPt);
    const int numEdgesInnerFace = mdlTopo.edges.size() - 4;
    const int numLoopsInnerFace = 1;
    int loopFirstEdgeIdx[1] = {0};
    int j = 4;
    for (int i = 0; i < numEdgesInnerFace; i++) {
      mdlTopo.faceDirs.push_back(faceDirectionFwd); // clockwise
      mdlTopo.faceEdges.push_back(mdlTopo.edges.at(j++));
    }
    mdlTopo.faces.push_back(GR_createFace(mdlTopo.region, numEdgesInnerFace,
          mdlTopo.faceEdges.data(),
          mdlTopo.faceDirs.data(),
          numLoopsInnerFace, loopFirstEdgeIdx,
          planarSurface, sameNormal));
    if(debug) {
      std::cout << "faces[1] area: " << GF_area(mdlTopo.faces[1], 0.2) << "\n";
    }
    assert(GF_area(mdlTopo.faces[1], 0.2) > 0);
  }
}

void createMesh(ModelTopo mdlTopo, std::string& meshFileName, pProgress progress, bool debug) {
  pMesh mesh = M_new(0, mdlTopo.model);
  pACase meshCase = MS_newMeshCase(mdlTopo.model);

  pModelItem domain = GM_domain(mdlTopo.model);
  // find the smallest size of the geometric model edges
  auto minGEdgeLen = std::numeric_limits<double>::max();
  for (int i = 0; i < mdlTopo.edges.size(); i++) {
    auto len = GE_length(mdlTopo.edges.at(i));
    if (len < minGEdgeLen)
      minGEdgeLen = len;
  }
  const auto contourMeshSize = minGEdgeLen * 128;
  const auto globMeshSize = contourMeshSize * 128;
  if(debug) {
    std::cout << "Min geometric model edge length: " << minGEdgeLen << std::endl;
    std::cout << "Contour absolute mesh size target: " << contourMeshSize
      << std::endl;
    std::cout << "Global absolute mesh size target: " << globMeshSize
      << std::endl;
  }
  MS_setMeshSize(meshCase, domain, 1, globMeshSize, NULL);
  for (int i = 4; i < mdlTopo.edges.size(); i++)
    MS_setMeshSize(meshCase, mdlTopo.edges.at(i), 1, contourMeshSize, NULL);

  {
    GFIter fIter = GM_faceIter(mdlTopo.model);
    pGFace modelFace;
    while (modelFace = GFIter_next(fIter)) {
      const double area = GF_area(modelFace, 0.2);
      if(debug) {
        std::cout << "face area: " << area << "\n";
      }
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
}

