#include "landIceMeshGen.h"
#include "Quadtree.h"
#include <map>

namespace TC {
void normalize(double x, double y, double& nx, double& ny) {
    double magnitude = std::fabs(x) + std::fabs(y); // taxicab "magnitude"
    if (magnitude == 0) {
        nx = ny = 0;
    } else {
        nx = x / magnitude;
        ny = y / magnitude;
    }
}

// Taxicab angle from positive x-axis (counterclockwise)
// Mapping taxicab unit directions to "angles" in 0 to 4 range
double angle(double x, double y) {
    if (x >= 0 && y >= 0)
        return y / (x + y); // First quadrant
    else if (x < 0 && y >= 0)
        return 1 + (-x) / ((-x) + y); // Second
    else if (x < 0 && y < 0)
        return 2 + (-y) / ((-x) + (-y)); // Third
    else
        return 3 + x / (x + (-y)); // Fourth
}

// Compute taxicab angle between two vectors
double angleBetween(double x1, double y1, double x2, double y2) {
    double nx1, ny1, nx2, ny2;
    TC::normalize(x1, y1, nx1, ny1);
    TC::normalize(x2, y2, nx2, ny2);
    double angle1 = TC::angle(nx1, ny1);
    double angle2 = TC::angle(nx2, ny2);
    double diff = std::abs(angle1 - angle2);
    if (diff > 4.0)
        diff = 8.0 - diff; // Ensure smallest angular distance
    return diff;
}

double radiansTo(double rad) {
    const double x = std::cos(rad);
    const double y = std::sin(rad);
    return TC::angle(x, y);
}

double degreesTo(double deg) {
    const double rad = deg * (M_PI / 180);
    return TC::radiansTo(rad);
}

} //end namespace TC

struct PlaneBounds {
  double minX;
  double maxX;
  double minY;
  double maxY;
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
GeomInfo readVtkGeom(std::string fname, bool debug) {
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

GeomInfo readJigGeom(std::string fname, bool debug) {
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

double getLengthSquared(double ax, double ay, double bx, double by) {
  double xDelta = std::abs(ax - bx);
  double yDelta = std::abs(ay - by);
  double length = xDelta * xDelta + yDelta * yDelta;
  return length;
}

bool isPtCoincident(double ax, double ay, double bx, double by,
                    double tolSquared = 1) {
  const double lengthSquared = getLengthSquared(ax, ay, bx, by);
  return (lengthSquared < tolSquared);
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

void convertMetersToKm(GeomInfo &geom) {
  std::transform(geom.vtx_x.cbegin(), geom.vtx_x.cend(), geom.vtx_x.begin(), [](double v) { return v * 0.001; });
  std::transform(geom.vtx_y.cbegin(), geom.vtx_y.cend(), geom.vtx_y.begin(), [](double v) { return v * 0.001; });
}

quadtree::Box<double> makeBoxAroundPt(double x, double y, double pad) {
  const double left = x-pad;
  const double bottom = y-pad; //see https://github.com/pvigier/Quadtree/issues/10#issuecomment-3058271986
  const double width = pad*2;
  const double height = pad*2;
  return {left, bottom, width, height};
}

bool isNumEdgesBtwnPtsGreaterThanOne(size_t small, size_t large, size_t firstPt, size_t lastPt) {
  assert(small<=large);
  if(large == lastPt && small == firstPt) {
    return false; //difference is one
  } else {
    return (large-small > 1);
  }
}

std::vector<std::pair<int,int>> removeNestedSegments(std::map<int,int> longPairs, int firstPt, int lastPt) {
  //FIXME: n^2 brute force search
  std::vector<std::pair<int,int>> unNestedPairs;
  for( auto it = longPairs.begin(); it != longPairs.end(); it++) {
    auto [aLow,aHigh] = *it;
    bool isNested = false;
    for( auto sit = longPairs.begin(); sit != longPairs.end(); sit++) {
      if( it == sit ) continue; //don't check against self
      auto [bLow,bHigh] = *sit;
      //check if a is within b
      if( aLow >= bLow && aHigh <= bHigh ) {
        isNested = true;
        break;
      }
    }
    if( !isNested ) {
      unNestedPairs.push_back({aLow,aHigh});
    }
  }
  if(true) {
    std::cout << "number of un-nested pairs found (brute force) " << unNestedPairs.size() << "\n";
    std::cout << "a_min, a_max\n";
    for(auto& [a,b] : unNestedPairs) {
      std::cout << a << ", " << b << "\n";
    }
    std::cout << "done\n";
  }
  return unNestedPairs;
}


//find pairs of points that are not consecutative, but are within some length
//tolerance of each other - mark these points as model vertices to help prevent
//intersecting bsplines
std::map<int,int> findNarrowChannels(GeomInfo &g, double coincidentVtxToleranceSquared, int firstPt, bool debug=false) {
  assert(g.vtx_x.size() >= firstPt);

  //use a quadtree
  struct Node
  {
    quadtree::Box<double> box;
    std::size_t id;
  };
  auto n = std::size_t(g.vtx_x.size());
  auto getBox = [](Node* node)
  {
    return node->box;
  };
  const auto bbox = getBoundingPlane(g);
  auto domain = quadtree::Box<double>(bbox.minX, bbox.minY, bbox.maxX-bbox.minX, bbox.maxY-bbox.minY);
  auto quadtree = quadtree::Quadtree<Node*, decltype(getBox), std::equal_to<Node*>, double>(domain, getBox);

  double padding = std::sqrt(coincidentVtxToleranceSquared)/2;
  std::vector<Node> nodes;
  for(size_t i = 4; i < g.vtx_x.size(); i++) {
    auto box = makeBoxAroundPt(g.vtx_x.at(i),g.vtx_y.at(i),padding);
    nodes.push_back({box,i});
  }
  for(auto& node : nodes) {
    quadtree.add(&node);
  }
  auto intersections = quadtree.findAllIntersections();
  if(true) {
    std::cout << "number of point pairs within " << std::sqrt(coincidentVtxToleranceSquared) << "km found: " << intersections.size() << '\n';
    std::cout << "pt0_id, pt0_x, pt0_y, pt1_id, pt1_x, pt1_y\n";
    for(auto& [a,b] : intersections) {
      const int distance = std::abs(static_cast<int>(a->id)-static_cast<int>(b->id));
      std::cout << a->id << ", " << g.vtx_x.at(a->id) << ", " << g.vtx_y.at(a->id) << ", "
        << b->id << ", " << g.vtx_x.at(b->id) << ", " << g.vtx_y.at(b->id) << ", "
        << distance << "\n";
    }
    std::cout << "done\n";
  }
  //remove consecutative pairs
  std::map<int,int> longPairs;
  const int lastPt = g.vtx_x.size()-1;
  for(auto& [a,b] : intersections) {
    const auto small = std::min(a->id, b->id);
    const auto large = std::max(a->id, b->id);
    if(isNumEdgesBtwnPtsGreaterThanOne(small, large, firstPt, lastPt)) {
      assert(longPairs.count(small) == 0);
      longPairs.insert({small, large});
    }
  }
  if(true) {
    std::cout << "longPairs " << longPairs.size() << "\n";
    std::cout << "id, min, max\n";
    int i=0;
    for(auto& [a,b] : longPairs) {
      std::cout << i++ << ", " << a << ", " << b << "\n";
    }
    std::cout << "done\n";
  }
  return longPairs;
}

GeomInfo cleanGeom(GeomInfo &dirty, double coincidentVtxToleranceSquared,
                      bool debug) {
  assert(checkVertexUse(dirty));
  // trying to check the the dirty geom has a chain of edges
  assert(dirty.numEdges == dirty.numVtx);
  // the first four edges form a loop
  assert(dirty.edges[0][0] == dirty.edges[3][1]);
  // the remaining edges form a loop
  assert(dirty.edges[4][0] == dirty.edges[dirty.numEdges - 1][1]);

  int numPtsRemoved = 0;
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
      numPtsRemoved++;
      if (debug) {
        std::cout << "coincident pt " << i - 1 << " (" << dirty.vtx_x[i - 1]
                  << ", " << dirty.vtx_y[i - 1] << ") " << i << " ("
                  << dirty.vtx_x[i] << ", " << dirty.vtx_y[i] << ")\n";
      }
    }
  }
  clean.numVtx = clean.vtx_x.size();
  std::cout << "removed " << numPtsRemoved << " coincident points\n";

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

pGEdge fitCurveToContourSimInterp(pGRegion region, pGVertex first, pGVertex last,
                         std::vector<double>& pts, bool debug=false) {
  assert(pts.size() % 3 == 0); //pts must contain coordinates x1,y1,z1, x2,y2,z2, ...
  const int numPts = pts.size()/3;
  assert(numPts > 1);
  pCurve curve;
  if( numPts == 2 || numPts == 3) {
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

//similar to scorec/tomms @ 2f97d13 (simapis-mod branch)
int onCurve(double tc_m1, double tc, double tc_p1, double onCurveAngleTol) {
  const double deg_angle_lower = onCurveAngleTol;
  const double deg_angle_upper = -onCurveAngleTol;
  const double tc_angle_lower = TC::degreesTo(deg_angle_lower);
  const double tc_angle_upper = TC::degreesTo(deg_angle_upper);
  if ((tc_m1>tc_angle_lower) && (tc_m1<tc_angle_upper) &&
      (tc   >tc_angle_lower) && (tc   <tc_angle_upper) &&
      (tc_p1>tc_angle_lower) && (tc_p1<tc_angle_upper)) {
    return 1;
  } else {
    return 0;
  }
}

int findStartingMdlVtx(std::vector<int>& isMdlVtx, const int offset) {
  const int mdlVtx = 1;
  auto it = std::find(isMdlVtx.begin()+offset, isMdlVtx.end(), mdlVtx);
  if( it == isMdlVtx.end()) {
    exit(EXIT_FAILURE);
    return -1;
  } else {
    return it - isMdlVtx.begin();
  }
}

void createEdges(ModelTopo& mdlTopo, GeomInfo& geom, std::vector<int>& isPtOnCurve, std::vector<int>& isMdlVtx, const bool debug, const int firstContourPt) {
  if(geom.numVtx <= firstContourPt) { // no internal contour
    return;
  }
  enum class State {MdlVtx = 0, OnCurve = 1, NotOnCurve = 2};
  enum class Action {Init, Advance, Line, Curve};
  typedef std::pair<State,Action> psa; // next state, action
  using func=std::function<psa(int pt)>;

  pGVertex firstMdlVtx;
  int firstPtIdx;
  int startingCurvePtIdx;
  pGVertex startingMdlVtx;
  std::vector<int> ptsOnCurve;
  func createCurve = [&](int pt) {
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

    auto edge = fitCurveToContourSimInterp(mdlTopo.region, startingMdlVtx, endMdlVtx, pts, true);
    mdlTopo.edges.push_back(edge);

    if (debug) {
      std::cerr << "edge " << mdlTopo.edges.size()
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
  func createLine = [&](int pt) {
    assert(ptsOnCurve.size() == 1);
    ptsOnCurve.push_back(pt);
    auto ignored = createCurve(pt);
    return psa{State::MdlVtx,Action::Line};
  };
  func advance = [&](int pt) {
    ptsOnCurve.push_back(pt);
    return psa{State::OnCurve,Action::Advance};
  };
  func createLineAndStartCurve = [&](int pt) {
    assert(ptsOnCurve.size() == 1);
    ptsOnCurve.push_back(pt);
    auto ignored = createCurve(pt);
    return psa{State::OnCurve,Action::Line};
  };
  func createCurveFromPriorPt = [&](int pt) {
    if(ptsOnCurve.size() == 1 ) {
      return createLine(pt);
    } else if (ptsOnCurve.size() >= 2 ) {
      //we are not adding the current point yet, so there must be
      //at least two points in the list to form a curve
      auto ignored = createCurve(ptsOnCurve.back());
      auto ret = createLine(pt);
      return ret;
    } else {
      std::cerr << "createCurveFromPriorPt: no points on the curve.... exiting\n";
      exit(EXIT_FAILURE);
    }
  };
  func createCurveFromCurrentPt = [&](int pt) {
    ptsOnCurve.push_back(pt);
    auto ret = createCurve(pt);
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
    {{State::MdlVtx,State::OnCurve}, createLineAndStartCurve},
    {{State::MdlVtx,State::NotOnCurve}, createLine},
    {{State::OnCurve,State::MdlVtx}, createCurveFromCurrentPt},
    {{State::OnCurve,State::OnCurve}, advance},
    {{State::OnCurve,State::NotOnCurve}, createCurveFromPriorPt},
    {{State::NotOnCurve,State::MdlVtx}, fail},
    {{State::NotOnCurve,State::OnCurve}, fail},
    {{State::NotOnCurve,State::NotOnCurve}, fail}
  };

  firstPtIdx = startingCurvePtIdx = findStartingMdlVtx(isMdlVtx, firstContourPt);
  double vtx[3] = {geom.vtx_x[startingCurvePtIdx], geom.vtx_y[startingCurvePtIdx], 0};
  firstMdlVtx = startingMdlVtx = GR_createVertex(mdlTopo.region, vtx);
  mdlTopo.vertices.push_back(firstMdlVtx);
  ptsOnCurve.push_back(startingCurvePtIdx);

  State state = State::MdlVtx;
  int ptsVisited = 0; //don't count the first vertex until we close the loop
  int ptIdx = startingCurvePtIdx+1;
  while(ptsVisited < isMdlVtx.size()-firstContourPt) {
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
    if(ptIdx == isMdlVtx.size()-1) {
      ptIdx = firstContourPt; //wrap around
    } else {
      ptIdx++;
    }
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

void createMesh(ModelTopo mdlTopo, std::string& meshFileName, pProgress progress) {
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
  std::cout << "Min geometric model edge length: " << minGEdgeLen
    << std::endl;
  const auto contourMeshSize = minGEdgeLen * 128;
  const auto globMeshSize = contourMeshSize * 128;
  std::cout << "Contour absolute mesh size target: " << contourMeshSize
    << std::endl;
  std::cout << "Global absolute mesh size target: " << globMeshSize
    << std::endl;
  MS_setMeshSize(meshCase, domain, 1, globMeshSize, NULL);
  for (int i = 4; i < mdlTopo.edges.size(); i++)
    MS_setMeshSize(meshCase, mdlTopo.edges.at(i), 1, contourMeshSize, NULL);

  {
    GFIter fIter = GM_faceIter(mdlTopo.model);
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
}

std::tuple<std::vector<int>,std::vector<int>>
discoverTopology(GeomInfo& geom, double coincidentPtTolSquared, double angleTol, double onCurveAngleTol, bool debug) {
  if(geom.numVtx <= 4) { // no internal contour
    return {std::vector<int>(), std::vector<int>()};
  }
  if(debug) {
    std::cout << "tc(30) " << TC::degreesTo(30) << "\n";
    std::cout << "tc(60) " << TC::degreesTo(60) << "\n";
    std::cout << "tc(90) " << TC::degreesTo(90) << "\n";
    std::cout << "tc(120) " << TC::degreesTo(120) << "\n";
    std::cout << "tc(150) " << TC::degreesTo(150) << "\n";
    std::cout << "tc(180) " << TC::degreesTo(180) << "\n";
    std::cout << "tc(270) " << TC::degreesTo(270) << "\n";
    std::cout << "tc(-120) " << TC::degreesTo(-120) << "\n";
  }
  const double deg_angle_lower = angleTol;
  const double deg_angle_upper = -deg_angle_lower;
  const double tc_angle_lower = TC::degreesTo(deg_angle_lower);
  const double tc_angle_upper = TC::degreesTo(deg_angle_upper);
  std::cout << "deg_angle_lower " << deg_angle_lower <<
               " tc_angle_lower " << tc_angle_lower << "\n";
  std::cout << "deg_angle_upper " << deg_angle_upper <<
               " tc_angle_upper " << tc_angle_upper << "\n";
  std::cout << "numPts " << geom.numVtx-4 << " lastPt " << geom.numVtx << "\n";

  std::vector<double> angle;
  std::vector<int> isMdlVtx;
  std::vector<int> isPointOnCurve; //1: along a curve, 0: otherwise
  angle.reserve(geom.numVtx);
  isMdlVtx.reserve(geom.numVtx);
  isPointOnCurve.reserve(geom.numVtx);
  //hack: add data for the first four boundary verts so 'createEdges' indexing
  //matches the GeomInfo struct indexing
  for(int i=0; i<4; i++){
    angle.push_back(TC::degreesTo(90)); //hack - 90deg corners
    isMdlVtx.push_back(1); //hack - all model verts
    isPointOnCurve.push_back(0); //hack - not on curve
  }

  //compute angle and determine if each pt is a model vertex
  for(int i=4; i<geom.numVtx; i++) {
    int m1 = i-1;
    int p1 = i+1;
    if (i == 4) { //first point
      m1 = geom.numVtx - 1;
    } else if ( i == geom.numVtx-1 ) { //last point
      p1 = 4;
    }
    const double norm_prev_x = geom.vtx_x[m1] - geom.vtx_x[i];
    const double norm_prev_y = geom.vtx_y[m1] - geom.vtx_y[i];
    const double norm_next_x = geom.vtx_x[p1] - geom.vtx_x[i];
    const double norm_next_y = geom.vtx_y[p1] - geom.vtx_y[i];
    const double tc_angle = TC::angleBetween(norm_prev_x, norm_prev_y, norm_next_x, norm_next_y);
    angle.push_back(tc_angle);
    isMdlVtx.push_back(tc_angle < tc_angle_lower || tc_angle > tc_angle_upper);
  }

  const int firstPt = 4;
  auto narrowPtPairs = findNarrowChannels(geom, coincidentPtTolSquared, firstPt);

  //mark the narrow points as model vertices
  for(auto& [a,b] : narrowPtPairs) {
    isMdlVtx.at(a) = 1;
    isMdlVtx.at(b) = 1;
  }

  for (int j = 4;j < geom.numVtx; ++j) {
    int m1 = j-1;
    int p1 = j+1;
    if (j == 4) { //first point
      m1 = geom.numVtx - 1;
    } else if ( j == geom.numVtx-1 ) { //last point
      p1 = 4;
    }
    const double tc_m1 = angle.at(m1);
    const double tc = angle.at(j);
    const double tc_p1 = angle.at(p1);
    const auto on = onCurve(tc_m1, tc, tc_p1, onCurveAngleTol);
    isPointOnCurve.push_back(on);
  }

  if(debug) {
    std::cout << "x,y,z,isOnCurve,angle,isMdlVtx\n";
    for (int j = 0;j < isPointOnCurve.size(); j++) {
      std::cout << geom.vtx_x.at(j) << ", " << geom.vtx_y.at(j) << ", " << 0
        << ", " << isPointOnCurve.at(j) << ", " << angle.at(j)
        << ", " << isMdlVtx.at(j) << "\n";
    }
    std::cout << "done\n";
  }

  //find points marked as on a curve that have no
  // adjacent points that are also marked as on the curve
  //first point
  const int firstContourPt = 4;
  if( isPointOnCurve.back() == 0 &&
      isPointOnCurve.at(firstContourPt) == 1 &&
      isPointOnCurve.at(firstContourPt+1) == 0 ) {
    isPointOnCurve.at(firstContourPt) = 0;
  }
  //interior
  for (int j = firstContourPt+1;j < isPointOnCurve.size(); j++) {
    if( isPointOnCurve.at(j-1) == 0 &&
        isPointOnCurve.at(j) == 1 &&
        isPointOnCurve.at(j+1) == 0 ) {
      isPointOnCurve.at(j) = 0;
    }
  }
  //last point
  const int last = isPointOnCurve.size()-1;
  if( isPointOnCurve.at(last-1) == 0 &&
      isPointOnCurve.at(last) == 1 &&
      isPointOnCurve.at(firstContourPt) == 0 ) {
    isPointOnCurve.at(last) = 0;
  }

  if(debug) {
    std::cout << "x,y,z,isOnCurveMod,angle,isMdlVtx\n";
    for (int j = 0;j < isPointOnCurve.size(); j++) {
      std::cout << geom.vtx_x.at(j) << ", " << geom.vtx_y.at(j) << ", " << 0
        << ", " << isPointOnCurve.at(j) << ", " << angle.at(j)
        << ", " << isMdlVtx.at(j) << "\n";
    }
    std::cout << "doneMod\n";
  }
  return {isPointOnCurve,isMdlVtx};
}

void createFaces(ModelTopo& mdlTopo, GeomInfo& geom) {
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
    std::cout << "faces[0] area: " << GF_area(mdlTopo.faces[0], 0.2) << "\n";
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
    std::cout << "faces[0] area: " << GF_area(mdlTopo.faces[0], 0.2) << "\n";
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
    std::cout << "faces[1] area: " << GF_area(mdlTopo.faces[1], 0.2) << "\n";
    assert(GF_area(mdlTopo.faces[1], 0.2) > 0);
  }
}


