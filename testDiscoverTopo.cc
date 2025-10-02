#include "modelGen2d.h"
#include "curveReader.h"
#include <numeric>

int main(int argc, char **argv) {
  const int numExpectedArgs = 2;
  if (argc != numExpectedArgs) {
    std::cerr << "Usage: <input csv file>\n";
    std::cerr << "input csv file with the following columns: "
                 "x,y,z,isOnCurve,angle,isMdlVtx\n";
    return 1;
  }
  assert(argc == numExpectedArgs);

  std::string filename = argv[1];
  std::cout << "input csv file: " << argv[1] << "\n";

  auto curveInfo = CurveReader::readCurveInfo(filename);
  CurveReader::printCurveInfo(curveInfo);
  ModelFeatures features;
  features.outer.numVtx = 4;
  for(int i=0; i < 4; i++) {
    features.outer.vtx_x.push_back(curveInfo.x[i]);
    features.outer.vtx_y.push_back(curveInfo.y[i]);
  }
  features.inner.numVtx = curveInfo.isMdlVtx.size()-4;
  for(int i=0, j=4; j < curveInfo.isMdlVtx.size(); j++, i++) {
    features.inner.vtx_x.push_back(curveInfo.x[j]);
    features.inner.vtx_y.push_back(curveInfo.y[j]);
  }

  auto coincidentPtTolSquared = 1.0;
  auto angleTol = 120.0;
  auto onCurveAngleTol = 40.0;
  auto debug = true;
  auto [isPointOnCurve, isMdlVtx] = discoverTopology(features.inner, coincidentPtTolSquared, angleTol, onCurveAngleTol, debug);
  auto numMdlVerts = std::accumulate(isMdlVtx.begin(), isMdlVtx.end(), 0);
  std::cout << "number of model vertices: " << numMdlVerts << std::endl;
  return 0;
}
