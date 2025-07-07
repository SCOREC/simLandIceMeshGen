#ifndef CURVEREADER_H
#define CURVEREADER_H
#include <vector>

namespace CurveReader {
struct CurveInfo {
    std::vector<double> x; 
    std::vector<double> y;
    std::vector<int> isOnCurve;
    std::vector<int> isMdlVtx;
    std::vector<double> tcAngle;
};

CurveInfo readCurveInfo(const std::string& filename);
void printCurveInfo(CurveInfo& c);
};

#endif
