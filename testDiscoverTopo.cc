#include "landIceMeshGen.h"
#include "curveReader.h"

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
  GeomInfo geom;
  geom.numVtx = curveInfo.isMdlVtx.size();
  geom.vtx_x = curveInfo.x;
  geom.vtx_y = curveInfo.y;

  auto coincidentPtTolSquared = 1.0;
  auto angleTol = 120.0;
  auto onCurveAngleTol = 40.0;
  auto debug = true;
  auto [isPointOnCurve, isMdlVtx] = discoverTopology(geom, coincidentPtTolSquared, angleTol, onCurveAngleTol, debug);
  return 0;
}
