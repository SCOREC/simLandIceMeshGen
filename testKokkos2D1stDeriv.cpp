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

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecutionSpace::memory_space;

int main(int argc, char* argv[]) {
  int retVal;
  if (argc != 3) {
    std::cerr<< "Input arguments: <input csv file> <expected curve length>" << std::endl;
    std::cerr << "input csv need these columns: ";
    std::cerr << "coordinate x, coordinate y, coordinate z,isOnCurve,angle,isMdlVtx" << std::endl;
    return 1;
  }

  Kokkos::initialize(argc, argv);
  {
    std::string inputCSV = argv[1];
    int extensionPos = inputCSV.rfind(".");
    int slashPos = inputCSV.rfind("/");
    std::string fileNameNoExt = inputCSV.substr(slashPos+1, extensionPos);
    double expectedCurveLength = std::stod(argv[2]);
    auto curve = CurveReader::readCurveInfo(inputCSV);

    //Construct BSpline2d object
    SplineInterp::BSpline2d serialBSP;
    if (curve.x.size() == 2) {
      serialBSP = SplineInterp::attach_piecewise_linear_curve(curve.x, curve.y);
    } else {
      serialBSP = SplineInterp::fitCubicSplineToPoints(curve.x, curve.y);
    }

    std::vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
    int order;

    //Get the info from serial spline, we will feed this to kokkos spline
    serialBSP.x.getpara(order, ctrlPtsX, knots, weight);
    serialBSP.y.getpara(order, ctrlPtsY, knots, weight);

    BSplineKokkos2D<ExecutionSpace> kokkosBSP(order, ctrlPtsX, ctrlPtsY, knots);
    
    std::vector<double> evalAt = {0, 0.2, 0.41, 0.5, 0.66, 0.73, 0.75, 0.89, 0.94, 1};

    //Create a derivCSR that will be used during deriv eval
    const size_t paraSize = evalAt.size();
    const size_t splineIdxSize = 3;
    BSplineKokkos2D<ExecutionSpace>::CSR kokkosBSPCSR(splineIdxSize, paraSize);
    auto valsMirror = Kokkos::create_mirror_view(kokkosBSPCSR.paraCoor);
    for (int i = 0; i < evalAt.size(); i++) {
      valsMirror(i) = evalAt[i];
    }
    auto splineIdxMirror = Kokkos::create_mirror_view(kokkosBSPCSR.splineIdx);
    for (int i = 0; i < 3; i++) {
      splineIdxMirror(i) = 0;
    }
    auto offsetMirror = Kokkos::create_mirror_view(kokkosBSPCSR.offset);

    offsetMirror(0) = 0;
    offsetMirror(1) = 3;
    offsetMirror(2) = 4;
    offsetMirror(3) = 10;
    
    Kokkos::deep_copy(kokkosBSPCSR.splineIdx,splineIdxMirror);
    Kokkos::deep_copy(kokkosBSPCSR.offset, offsetMirror);
    Kokkos::deep_copy(kokkosBSPCSR.paraCoor, valsMirror);
    
    //Calling our batch 1st derivative function
    Kokkos::View<double*[2], MemSpace> res("deriv result", offsetMirror(offsetMirror.extent(0)-1));
    res = kokkosBSP.eval1stDeriv(kokkosBSPCSR);
    auto mvRes = Kokkos::create_mirror_view(res);
    Kokkos::deep_copy(mvRes, res);
    
    for (int i = 0; i < 10; i++) {
      double derivX = serialBSP.x.evalFirstDeriv(evalAt[i]);
      double derivY = serialBSP.y.evalFirstDeriv(evalAt[i]);
    
      double xDiff = std::fabs(derivX) - std::fabs(mvRes(i, 0));
      double yDiff = std::fabs(derivY) - std::fabs(mvRes(i,1));
      
      //std::cout << "Serial x, y: " << derivX << "|" << derivY << std::endl;
      //std::cout << "Kokkos x, y: " << mvRes(i, 0) << "|" << mvRes(i, 1) << std::endl;

      if (xDiff > EPSILON || yDiff > EPSILON) {
        std::cout << "Test " << i+1 << " failed, eval at: " << evalAt[i] << std::endl;
        std::cout << "Difference: x = " << xDiff << " y = " << yDiff << std::endl;
        std::cout << "SERIAL 1st deriv: x = " << derivX << " y = " << derivY << std::endl;
        std::cout << "KOKKOS 1st deriv: x = " << mvRes(i, 0) << " y = " << mvRes(i, 1) << std::endl;
        retVal = 1;
      }
    }
  }
  Kokkos::finalize();
  return retVal;
}
