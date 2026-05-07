#include <Kokkos_Core.hpp>
#include "BSplineKokkos.h"
#include "BSpline.h"
#include <vector>

//These imports are from testSplineFitting.cc
#include "splineInterpolation.h"
#include "curveReader.h"
#include <cassert>
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <numeric>


template<typename MemSpace>
void printDView(Kokkos::View<double*, MemSpace>& dView) {
  auto host_dView = Kokkos::create_mirror_view(dView);
  Kokkos::deep_copy(host_dView, dView);
  for (int i = 0; i < host_dView.extent(0); i++) {
    std::cout << host_dView(i) << std::endl;
  }
}

template<typename MemSpace>
void printIView(Kokkos::View<int*, MemSpace>& iView) {
  auto host_iView = Kokkos::create_mirror_view(iView);
  Kokkos::deep_copy(host_iView, iView);
  for (int i = 0; i < host_iView.extent(0); i++) {
    std::cout << host_iView(i) << std::endl;
  }
}

template <typename MemSpace>
void printSpline(BSplineKokkos<MemSpace> spline) {
  std::cout << "Spline content" << std::endl;
  Kokkos::View<double*, MemSpace> dView = spline.getCtrlPts();
  std::cout << "CtrlPts" << std::endl;
  printDView(dView);
  Kokkos::View<int*, MemSpace> iView = spline.getCPOffset();
  std::cout << "CtrlPts Offset" << std::endl;
  printIView(iView);
  dView = spline.getKnots();
  std::cout << "Knots" << std::endl;
  printDView(dView);
  iView = spline.getKnotsOffset();
  std::cout << "Knots Offset" << std::endl;
  printIView(iView);
  iView = spline.getOrder();
  std::cout << "Order" << std::endl;
  printIView(iView);
}

void printVector(const std::string& name, std::vector<double>& vector) {
  std::cout << name << ": " << vector.size() << std::endl;
  for (int i = 0; i < vector.size(); i++) {
    std::cout << vector[i] << std::endl;
  }
}

double EPSILON = 1e-12;

int main(int argc, char* argv[]) {
  int retVal = 0;                                             
  if (argc != 3) {
    //Check for the arguments needed 
    std::cerr<< "Input arguments: <input csv file> <expected curve length>" << std::endl;
    std::cerr << "input csv need these columns: ";
    std::cerr << "coordinate x, coordinate y, coordinate z,isOnCurve,angle,isMdlVtx" << std::endl;
    return 1;
  }
  Kokkos::initialize(argc, argv);
  {
    //Check for what memory space we should use
    //Sticking to just cuda space or host space for now
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

    BSplineKokkos<ExecutionSpace> kokkosBSP(order, ctrlPtsX, ctrlPtsY, knots);
   //printSpline(kokkosBSP);
    std::vector<double> evalAt = {0, 0.2, 0.41, 0.5, 0.66, 0.73, 0.75, 0.89, 0.94, 1};

    for (int i = 0; i < 10; i++) {
      double derivX = serialBSP.x.evalFirstDeriv(evalAt[i]);
      double derivY = serialBSP.y.evalFirstDeriv(evalAt[i]);
      std::vector<double> kokkos1stDeriv = kokkosBSP.eval1stDeriv({evalAt[i]}, 0);
      double kokkosDerivX = kokkos1stDeriv[0];
      double kokkosDerivY = kokkos1stDeriv[1];

      double derivXDiff = std::fabs(derivX) - std::fabs(kokkosDerivX);
      double derivYDiff = std::fabs(derivY) - std::fabs(kokkosDerivY);

      if (derivXDiff > EPSILON || derivYDiff > EPSILON) {
        std::cout << "Test " << i+1 << " failed, eval at: " << evalAt[i] << std::endl;
        std::cout << "EPSILON = " << EPSILON << std::endl;
        std::cout << "SERIAL/KOKKOS DIFFERENCE: xcoor = " << derivXDiff << " ycoor = " << derivYDiff << std::endl;
        std::cout << "SERIAL 1st Derivative: x = " << derivX << " y = " << derivY << std::endl;
        std::cout << "KOKKOS 1st Derivative: x = " << kokkosDerivX << " y = " << kokkosDerivY << std::endl;
        retVal = 1;
      }
    }

  }
  Kokkos::finalize();
  return retVal;
}
