#include "modelGen2d.h"
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

typedef std::array<double, 2> Vec2d;

double crossProduct2d(const Vec2d& ab, const Vec2d& cd) {
  //a*d - c*b
  return ab[0]*cd[1]-cd[0]*ab[1];
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

void convertMetersToKm(GeomInfo &geom) {
  std::transform(geom.vtx_x.cbegin(), geom.vtx_x.cend(), geom.vtx_x.begin(), [](double v) { return v * 0.001; });
  std::transform(geom.vtx_y.cbegin(), geom.vtx_y.cend(), geom.vtx_y.begin(), [](double v) { return v * 0.001; });
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
  if(debug)
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
  clean.firstContourPt = 4;

  return clean;
}

OnCurve::OnCurve(double onCurveAngleTol, bool isDebug) :
  deg_angle_lower(onCurveAngleTol),
  deg_angle_upper(-onCurveAngleTol),
  tc_angle_lower(TC::degreesTo(deg_angle_lower)),
  tc_angle_upper(TC::degreesTo(deg_angle_upper)),
  debug(isDebug)
{
  if(debug) {
    std::cout << "OnCurve deg_angle_lower " << deg_angle_lower <<
      " tc_angle_lower " << tc_angle_lower << "\n";
    std::cout << "OnCurve deg_angle_upper " << deg_angle_upper <<
      " tc_angle_upper " << tc_angle_upper << "\n";
  }
}

//similar to scorec/tomms @ 2f97d13 (simapis-mod branch)
int OnCurve::operator()(double tc_m1, double tc, double tc_p1) {
  if ((tc_m1>tc_angle_lower) && (tc_m1<tc_angle_upper) &&
      (tc   >tc_angle_lower) && (tc   <tc_angle_upper) &&
      (tc_p1>tc_angle_lower) && (tc_p1<tc_angle_upper)) {
    return 1;
  } else {
    return 0;
  }
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

//find pairs of points that are not consecutative, but are within some length
//tolerance of each other - mark these points as model vertices to help prevent
//intersecting bsplines
std::map<int,int> findNarrowChannels(GeomInfo &g, double coincidentVtxToleranceSquared, bool debug=false) {
  assert(g.vtx_x.size() >= g.firstContourPt);

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
  if(debug) {
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
    if(isNumEdgesBtwnPtsGreaterThanOne(small, large, g.firstContourPt, lastPt)) {
      assert(longPairs.count(small) == 0);
      longPairs.insert({small, large});
    }
  }
  if(debug) {
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

void writeToCSV(std::string fname, GeomInfo& geom,
    std::vector<double>& angle,
    std::vector<int>& isPointOnCurve,
    std::vector<int>& isMdlVtx) {
  std::ofstream csv(fname);
  assert(csv.is_open());
  csv << "x,y,z,isOnCurve,angle,isMdlVtx\n";
  for (int j = 0;j < isPointOnCurve.size(); j++) {
    csv << geom.vtx_x.at(j) << ", " << geom.vtx_y.at(j) << ", " << 0
      << ", " << isPointOnCurve.at(j) << ", " << angle.at(j)
      << ", " << isMdlVtx.at(j) << "\n";
  }
  csv.close();
}

std::tuple<std::vector<int>,std::vector<int>>
discoverTopology(GeomInfo& geom, double coincidentPtTolSquared, double angleTol, double onCurveAngleTol, bool debug) {
  if(geom.numVtx <= geom.firstContourPt) { // no internal contour
    return {std::vector<int>(), std::vector<int>()};
  }
  const double deg_angle_lower = angleTol;
  const double deg_angle_upper = -deg_angle_lower;
  const double tc_angle_lower = TC::degreesTo(deg_angle_lower);
  const double tc_angle_upper = TC::degreesTo(deg_angle_upper);
  if(debug) {
    std::cout << "tc(30) " << TC::degreesTo(30) << "\n";
    std::cout << "tc(60) " << TC::degreesTo(60) << "\n";
    std::cout << "tc(90) " << TC::degreesTo(90) << "\n";
    std::cout << "tc(120) " << TC::degreesTo(120) << "\n";
    std::cout << "tc(150) " << TC::degreesTo(150) << "\n";
    std::cout << "tc(180) " << TC::degreesTo(180) << "\n";
    std::cout << "tc(270) " << TC::degreesTo(270) << "\n";
    std::cout << "tc(-120) " << TC::degreesTo(-120) << "\n";
    std::cout << "deg_angle_lower " << deg_angle_lower <<
                 " tc_angle_lower " << tc_angle_lower << "\n";
    std::cout << "deg_angle_upper " << deg_angle_upper <<
                 " tc_angle_upper " << tc_angle_upper << "\n";
    std::cout << "numPts " << geom.numVtx-geom.firstContourPt << " lastPt " << geom.numVtx << "\n";
  }

  std::vector<double> angle;
  std::vector<int> isMdlVtx;
  std::vector<int> isPointOnCurve; //1: along a curve, 0: otherwise
  angle.reserve(geom.numVtx);
  isMdlVtx.reserve(geom.numVtx);
  isPointOnCurve.reserve(geom.numVtx);
  //hack: add data for the first four boundary verts so 'createEdges' indexing
  //matches the GeomInfo struct indexing
  for(int i=0; i<geom.firstContourPt; i++){
    angle.push_back(TC::degreesTo(90)); //hack - 90deg corners
    isMdlVtx.push_back(1); //hack - all model verts
    isPointOnCurve.push_back(0); //hack - not on curve
  }

  //compute angle and determine if each pt is a model vertex
  for(int i=geom.firstContourPt; i<geom.numVtx; i++) {
    const int m1 = geom.getPrevPtIdx(i);
    const int p1 = geom.getNextPtIdx(i);
    const double norm_prev_x = geom.vtx_x[m1] - geom.vtx_x[i];
    const double norm_prev_y = geom.vtx_y[m1] - geom.vtx_y[i];
    const double norm_next_x = geom.vtx_x[p1] - geom.vtx_x[i];
    const double norm_next_y = geom.vtx_y[p1] - geom.vtx_y[i];
    const double tc_angle = TC::angleBetween(norm_prev_x, norm_prev_y, norm_next_x, norm_next_y);
    angle.push_back(tc_angle);
    isMdlVtx.push_back(tc_angle < tc_angle_lower || tc_angle > tc_angle_upper);
  }

  //mark points that are on smooth curves
  OnCurve onCurve(onCurveAngleTol);
  const double smoothAngle = (onCurve.getLowerTolTC()+onCurve.getUpperTolTC())/2;
  for (int j = geom.firstContourPt;j < geom.numVtx; ++j) {
    const int m1 = geom.getPrevPtIdx(j);
    const int p1 = geom.getNextPtIdx(j);
    const double tc_m1 = isMdlVtx.at(m1) ? smoothAngle : angle.at(m1); //ignore the point if it is a model vtx
    const double tc = angle.at(j);
    const double tc_p1 = isMdlVtx.at(p1) ? smoothAngle : angle.at(p1); //ignore the point if it is a model vtx
    const auto on = onCurve(tc_m1, tc, tc_p1);
    isPointOnCurve.push_back(on);
  }

  if(debug) {
    writeToCSV("init.csv", geom, angle, isPointOnCurve, isMdlVtx);
  }

  //mark pairs of points that are within a tolerance of each other as not on
  //smooth curves to force a linear spline through them
  auto narrowPtPairs = findNarrowChannels(geom, coincidentPtTolSquared);
  for(auto& [a,b] : narrowPtPairs) {
    isPointOnCurve.at(a) = 0;
    isPointOnCurve.at(b) = 0;
  }

  if(debug) {
    writeToCSV("narrowChannels.csv", geom, angle, isPointOnCurve, isMdlVtx);
  }

  return {isPointOnCurve, isMdlVtx};
}