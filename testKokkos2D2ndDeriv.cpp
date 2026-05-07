#include <Kokkos_Core.hpp>
#include "BSplineKokkos2D.h"
#include "BSpline.h"
#include <vector>

#include "splineInterpolation.h"
#include "curveReader.h"
#include <cassert>
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <numeric>

double EPSILON = 1e-12;

int main(int argc, char* argv[]) {
  int retVal = 0;
  if (argc != 3) {
    std::cerr<< "Input arguments: <input csv file> <expected curve length>" << std::endl;
    std::cerr << "input csv need these columns: ";
    std::cerr << "coordinate x, coordinate y, coordinate z,isOnCurve,angle,isMdlVtx" << std::endl;
    return 1;
  }

  Kokkos::initialize(argc, argv);
  {
    #ifdef KOKKOS_ENABLE_CUDA
    #define MemSpace Kokkos::CudaSpace
    #endif
    #ifndef MemSpace
    #define MemSpace Kokkos::HostSpace
    #endif

    using ExecutionSpace = MemSpace::execution_space;

    std::string inputCSV = argv[1];
    int extensionPos = inputCSV.rfind(".");
    int slashPos = inputCSV.rfind("/");
    std::string fileNameNoExt = inputCSV.substr(slashPos+1, extensionPos);
    double expectedCurveLength = std::stod(argv[2]);
    auto curve = CurveReader::readCurveInfo(inputCSV);

    SplineInterp::BSpline2d serialBSP;
    if (curve.x.size() == 2) {
      serialBSP = SplineInterp::attach_piecewise_linear_curve(curve.x, curve.y);
    }
    else {
      serialBSP = SplineInterp::fitCubicSplineToPoints(curve.x, curve.y);
    }

    std::vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
    int order;
    serialBSP.x.getpara(order, ctrlPtsX, knots, weight);
    serialBSP.y.getpara(order, ctrlPtsY, knots, weight);

    BSplineKokkos2D<ExecutionSpace> kokkosBSP(order, ctrlPtsX, ctrlPtsY, knots);

    std::vector<double> evalAt = {0, 0.2, 0.41, 0.5, 0.66, 0.73, 0.75, 0.89, 0.94, 1};
    for (int i = 0; i < 10; i++) {
      double derivX = serialBSP.x.evalSecondDeriv(evalAt[i]);
      double derivY = serialBSP.y.evalSecondDeriv(evalAt[i]);

      Kokkos::View<double*, MemSpace> res("Result", 2);
      Kokkos::View<double*, MemSpace> xVals("paraCoor", 1);
      auto mvXVals = Kokkos::create_mirror_view(xVals);
      mvXVals(0) = evalAt[i];
      Kokkos::deep_copy(xVals, mvXVals);

      res = kokkosBSP.eval2ndDeriv(xVals, 0);
      auto mvRes = Kokkos::create_mirror_view(res);
      Kokkos::deep_copy(mvRes, res);

      double xDiff = std::fabs(derivX - mvRes(0));
      double yDiff = std::fabs(derivY - mvRes(1));

      if (xDiff > EPSILON || yDiff > EPSILON) {
        std::cout << "Test " << i+1 << " failed, eval at: " << evalAt[i] << std::endl;
        std::cout << "Difference: x = " << xDiff << " y = " << derivY << std::endl;
        std::cout << "SERIAL 2nd deriv: x = " << derivX << ", " << derivY << std::endl;
        std::cout << "KOKKOS 2nd deriv: x = " << mvRes(0) << ", " << mvRes(1) << std::endl;
        retVal = 1;
      }
    }
  }
  Kokkos::finalize();
  return retVal;
}
