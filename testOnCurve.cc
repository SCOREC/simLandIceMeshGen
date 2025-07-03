#include "landIceMeshGen.h"
#include "curveReader.h"

int main(int argc, char **argv) {
  const int numExpectedArgs = 5;
  if (argc != numExpectedArgs) {
    std::cerr << "Usage: <input csv file> <ptIdx> <isOnCurve> <onCurveAngleTol>\n";
    std::cerr << "input csv file with the following columns: "
                 "x,y,z,isOnCurve,angle,isMdlVtx\n";
    std::cerr << "ptIdx: the point to check\n";
    std::cerr << "isOnCurve: the expected result for the point being checked\n";
    std::cerr << "onCurveAngleTol: the angle being passed to onCurve(...)\n";
    return 1;
  }
  assert(argc == numExpectedArgs);

  std::string filename = argv[1];
  const auto ptIdx = std::atoi(argv[2]);
  const auto isOnCurve = std::atoi(argv[3]);
  const auto onCurveAngleTol = std::atof(argv[4]);
  std::cout << "input csv file: " << argv[1] << " "
            << "onCurveAngleTol : " << onCurveAngleTol << "\n";

  const auto debug = true;

  auto curveInfo = CurveReader::readCurveInfo(filename);
  CurveReader::printCurveInfo(curveInfo);

  Point m2{curveInfo.x.at(ptIdx-2), curveInfo.y.at(ptIdx-2), 0.0};
  Point m1{curveInfo.x.at(ptIdx-1), curveInfo.y.at(ptIdx-1), 0.0};
  Point m0{curveInfo.x.at(ptIdx),   curveInfo.y.at(ptIdx),   0.0};
  Point p1{curveInfo.x.at(ptIdx+1), curveInfo.y.at(ptIdx+1), 0.0};
  Point p2{curveInfo.x.at(ptIdx+2), curveInfo.y.at(ptIdx+2), 0.0};
  const auto on = onCurve(m2, m1, m0, p1, p2, onCurveAngleTol);
  std::cerr << "pt " << ptIdx 
            << " (" << curveInfo.x.at(ptIdx) << ", " << curveInfo.y.at(ptIdx) << ")"
            << " isOnCurve = " << on << "\n";
  if( on != isOnCurve ) {
    return 1;
  } else {
    return 0;
  }
}

