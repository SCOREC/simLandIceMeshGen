//Constructor test for new representation of Kokkos spline
//We will be comparing this against the serial version
#include <Kokkos_Core.hpp>
#include "BSplineKokkos2D.h"
#include "BSpline.h"
#include "splineInterpolation.h"
#include "curveReader.h"

#include <vector>
#include <iostream>
#include <string>
#include <cassert>
#include <math.h>
#include <fstream>
#include <numeric>

//For checking if the content of the splines are correct
//

double EPSILON = 1e-12;

int main(int argc, char* argv[]) {
  //We check how many arguments are given
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

    //Construct BSpline2d object
    SplineInterp::BSpline2d serialBSP;
    if (curve.x.size() == 2) {
      serialBSP = SplineInterp::attach_piecewise_linear_curve(curve.x, curve.y);
    } else {
      serialBSP = SplineInterp::fitCubicSplineToPoints(curve.x, curve.y);
    }

    std::vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
    int order;

    serialBSP.x.getpara(order, ctrlPtsX, knots, weight);
    serialBSP.y.getpara(order, ctrlPtsY, knots, weight);

    BSplineKokkos2D<ExecutionSpace> kokkosBSP(order, ctrlPtsX, ctrlPtsY, knots);

    auto intView = Kokkos::create_mirror_view(kokkosBSP.getOrder());
    Kokkos::deep_copy(intView, kokkosBSP.getOrder());
    //Testing the order initialization
    double diff = std::fabs(order - intView(0));
    if (diff > EPSILON) {
      std::cout << "Order difference : " << diff << std::endl;
      std::cout << "Serial: " << order << " Kokkos: " << intView(0) << std::endl;
      retVal = 1;
    }

    auto double2DView = Kokkos::create_mirror_view(kokkosBSP.getCtrlPts());
    Kokkos::deep_copy(double2DView, kokkosBSP.getCtrlPts());

    //Testing ctrlPts initialization
    double xDiff, yDiff;
    for (int i = 0; i < ctrlPtsX.size(); i++) {
      xDiff = std::fabs(ctrlPtsX[i] - double2DView(i, 0));
      yDiff = std::fabs(ctrlPtsY[i] - double2DView(i, 1));
      if (xDiff > EPSILON || yDiff > EPSILON) {
        std::cout << "CtrlPts difference: x = " << xDiff << " y = " << yDiff << std::endl;
        std::cout << "Serial: " << ctrlPtsX[i] << ", " << ctrlPtsY[i] << std::endl;
        std::cout << "Kokkos: " << double2DView(i, 0) << ", " << double2DView(i, 1) << std::endl;
        retVal = 1;
      }
    }

    //Test knots initialization
    auto doubleView = Kokkos::create_mirror_view(kokkosBSP.getKnots());
    Kokkos::deep_copy(doubleView, kokkosBSP.getKnots());
    for (int i = 0; i < knots.size(); i++) {
      diff = std::fabs(knots[i] - doubleView(i));
      if (diff > EPSILON) {
        std::cout << "Knots difference: " << diff << std::endl;
        std::cout << "Serial: " << knots[i] << std::endl;
        std::cout << "Kokkos: " << doubleView(i) << std::endl;
        retVal = 1;
      }
    }

    //Test 1st deriv coef initialization
    Kokkos::View<double*[2], MemSpace> coef ("getCPCoe2", kokkosBSP.getCP1stD().extent(0));
    coef = kokkosBSP.getCP1stD();
    auto mv_coef = Kokkos::create_mirror_view(coef);

    Kokkos::deep_copy(mv_coef, coef);
    std::vector<double> serialCoefX = serialBSP.x.getCtrlPts_1st();
    std::vector<double> serialCoefY = serialBSP.y.getCtrlPts_1st();

    for (int i = 0; i < serialCoefX.size(); i++) {
      double diffX = std::fabs(serialCoefX[i] - mv_coef(i, 0));
      double diffY = std::fabs(serialCoefY[i] - mv_coef(i, 1));
      if (diffX > EPSILON || diffY > EPSILON) {
        std::cout << "diffX: " << diffX << std::endl;
        std::cout << "diffY: " << diffY << std::endl;
        std::cout << "Serial 1st coef: " << serialCoefX[i] << ", " << serialCoefY[i] << std::endl;
        std::cout << "Kokkos 1st coef: " << mv_coef(i, 0)<< ", " << mv_coef(i, 1) << std::endl;
        retVal = 1;
      }
    }

    Kokkos::View<double*[2], MemSpace>coef2 ("getCPCoe2", kokkosBSP.getCP2ndD().extent(0));
    coef2 = kokkosBSP.getCP2ndD();
    mv_coef = Kokkos::create_mirror_view(coef2);

    Kokkos::deep_copy(mv_coef, coef2);

    serialCoefX = serialBSP.x.getCtrlPts_2nd();
    serialCoefY = serialBSP.y.getCtrlPts_2nd();
    for (int i = 0; i < serialCoefX.size(); i++) {
      double diffX = std::fabs(serialCoefX[i] - mv_coef(i, 0));
      double diffY = std::fabs(serialCoefY[i] - mv_coef(i, 1));
      if (diffX > EPSILON || diffY > EPSILON) {
        std::cout << "diffX: " << diffX << std::endl;
        std::cout << "diffY: " << diffY << std::endl;
        std::cout << "Serial 2nd coef: " << serialCoefX[i] << ", " << serialCoefY[i] << std::endl;
        std::cout << "Kokkos 2nd coef: " << mv_coef(i, 0)<< ", " << mv_coef(i, 1) << std::endl;
        retVal = 1;
      }
    }

  }
  Kokkos::finalize();
  return retVal;
}

