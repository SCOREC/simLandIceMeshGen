#include "BSpline.h"
#include "BSplineKokkos2D.h"
#include "makeScalingTestData.h"
#include <Kokkos_Core.hpp>
#include <vector>

#include "curveReader.h"
#include "splineInterpolation.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>

double EPSILON = 1e-11;

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecutionSpace::memory_space;

int main(int argc, char *argv[]) {
  int retVal;
  if (argc != 4) {
    std::cout << "Input arguments: <number of splines> <number of points per "
                 "spline> <number of para coords to evaluate>"
              << std::endl;
    retVal = 1;
    return retVal;
  }
  Kokkos::initialize(argc, argv);
  {
    // Creating Scaling Data
    // Using Kokkos timer to keep track of time elapsed
    Kokkos::Timer timer;
    const int numSplines = std::atoi(argv[1]);
    const int ptsPerSpline = std::atoi(argv[2]);
    Kokkos::View<double *[2], MemSpace> pts("PointsOnCircle",
                                            numSplines * ptsPerSpline);
    pts = makeCircle(1000.0, 2000.0, 500.0, numSplines, ptsPerSpline);
    double dataGenTime = timer.seconds();

    timer.reset();
    // Copy the result to vectors for serial initialization
    auto ptsMirror = Kokkos::create_mirror_view(pts);
    Kokkos::deep_copy(ptsMirror, pts);
    std::vector<double> ptsX(numSplines * ptsPerSpline),
        ptsY(numSplines * ptsPerSpline);
    for (int i = 0; i < numSplines * ptsPerSpline; i++) {
      ptsX[i] = ptsMirror(i, 0);
      ptsY[i] = ptsMirror(i, 1);
    }
    double dataCopyTime = timer.seconds();

    timer.reset();
    // Initializing Serial BSpline
    SplineInterp::BSpline2d serialBSP;
    if (ptsX.size() == 2) {
      serialBSP = SplineInterp::attach_piecewise_linear_curve(ptsX, ptsY);
    } else {
      serialBSP = SplineInterp::fitCubicSplineToPoints(ptsX, ptsY);
    }
    double serialSplineCreationTime = timer.seconds();

    // Getting the data from serial spline, prepare for BSPlineKokkos2D
    // generation
    int order;
    std::vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
    serialBSP.x.getpara(order, ctrlPtsX, knots, weight);
    serialBSP.y.getpara(order, ctrlPtsY, knots, weight);

    // Initializing BSplineKokkos2D object
    timer.reset();
    BSplineKokkos2D<ExecutionSpace> kokkosBSP(order, ctrlPtsX, ctrlPtsY, knots);
    double kokkosSplineCreationTime = timer.seconds();

    /*-------- Starting 1st derivative test --------*/
    // Creating parametric coordinates to evaluate at
    int paraCoords = std::atoi(argv[3]);
    double incr = 1.0 / paraCoords;
    double val = 0.0;
    std::vector<double> evalAt(paraCoords);
    for (int i = 0; i < paraCoords; i++) {
      evalAt[i] = val;
      val += incr;
    }

    // Serial evaluation of 1st derivative
    std::vector<double> serialResX(evalAt.size());
    std::vector<double> serialResY(evalAt.size());

    timer.reset();

    for (int i = 0; i < evalAt.size(); i++) {
      serialResX[i] = serialBSP.x.evalFirstDeriv(evalAt[i]);
      serialResY[i] = serialBSP.y.evalFirstDeriv(evalAt[i]);
    }
    double serial1stDerivTime = timer.seconds();

    // Kokkos evaluation of 1st derivative
    // We will split this spline into 10 chunks to evaluate
    const size_t paraSize = evalAt.size();
    const size_t splineIdxSize = 10;

    timer.reset(); // Recording the time needed for CSR initialization
    BSplineKokkos2D<ExecutionSpace>::CSR kokkosBSPCSR(splineIdxSize, paraSize);
    auto valsMirror = Kokkos::create_mirror_view(kokkosBSPCSR.paraCoor);
    for (int i = 0; i < evalAt.size(); i++) {
      valsMirror(i) = evalAt[i];
    }
    auto splineIdxMirror = Kokkos::create_mirror_view(kokkosBSPCSR.splineIdx);
    for (int i = 0; i < splineIdxSize; i++) {
      splineIdxMirror(i) = 0;
    }
    auto offsetMirror = Kokkos::create_mirror_view(kokkosBSPCSR.offset);
    int offsetSize = paraSize / splineIdxSize;
    int curOffset = 0;
    for (int i = 0; i < splineIdxSize; i++) {
      offsetMirror(i) = curOffset;
      curOffset += offsetSize;
    }
    offsetMirror(offsetMirror.extent(0) - 1) = evalAt.size();
    Kokkos::deep_copy(kokkosBSPCSR.splineIdx, splineIdxMirror);
    Kokkos::deep_copy(kokkosBSPCSR.offset, offsetMirror);
    Kokkos::deep_copy(kokkosBSPCSR.paraCoor, valsMirror);
    double csrInitTime = timer.seconds();
    // Calling batch 1st derivative function
    Kokkos::View<double *[2], MemSpace> res("derivResult", evalAt.size());
    timer.reset();
    res = kokkosBSP.eval1stDeriv(kokkosBSPCSR);
    double kokkos1stDerivTime = timer.seconds();
    auto resMirror = Kokkos::create_mirror_view(res);
    Kokkos::deep_copy(resMirror, res);

    // Compare the serial and kokkos result
    timer.reset();
    for (int i = 0; i < evalAt.size(); i++) {
      double xDiff = std::fabs(resMirror(i, 0)) - std::fabs(serialResX[i]);
      double yDiff = std::fabs(resMirror(i, 1)) - std::fabs(serialResY[i]);
      if (xDiff > EPSILON || yDiff > EPSILON) {
        std::cout << "1st Deriv Test " << i + 1
                  << " failed, eval at: " << evalAt[i] << std::endl;
        std::cout << "Difference: x = " << xDiff << " y = " << yDiff
                  << std::endl;
        std::cout << "SERIAL 1st deriv: x = " << serialResX[i]
                  << " y = " << serialResY[i] << std::endl;
        std::cout << "KOKKOS 1st deriv: x = " << resMirror(i, 0)
                  << " y = " << resMirror(i, 1) << std::endl;
        retVal = 1;
      }
    }
    double verifyTime1stDeriv = timer.seconds();

    /*-------- End of 1st Deriv Test --------*/
    /*-------- Start of 2nd Deriv Test --------*/
    // Serial evaluation

    timer.reset();
    for (int i = 0; i < evalAt.size(); i++) {
      serialResX[i] = serialBSP.x.evalSecondDeriv(evalAt[i]);
      serialResY[i] = serialBSP.y.evalSecondDeriv(evalAt[i]);
    }
    double serial2ndDerivTime = timer.seconds();
    // Kokkos 2nd derivative function, using the existing csr
    timer.reset();
    res = kokkosBSP.eval2ndDeriv(kokkosBSPCSR);
    double kokkos2ndDerivTime = timer.seconds();
    Kokkos::deep_copy(resMirror, res);
    // Compare serial and kokkos result
    timer.reset();
    for (int i = 0; i < evalAt.size(); i++) {
      double xDiff = std::fabs(resMirror(i, 0)) - std::fabs(serialResX[i]);
      double yDiff = std::fabs(resMirror(i, 1)) - std::fabs(serialResY[i]);
      if (xDiff > EPSILON || yDiff > EPSILON) {
        std::cout << "2nd Deriv Test " << i + 1
                  << "  failed, eval at: " << evalAt[i] << std::endl;
        std::cout << "Difference: x = " << xDiff << " y = " << yDiff
                  << std::endl;
        std::cout << "SERIAL 2nd deriv: x = " << serialResX[i]
                  << " y = " << serialResY[i] << std::endl;
        std::cout << "KOKKOS 2nd deriv: x = " << resMirror(i, 0)
                  << " y = " << resMirror(i, 1) << std::endl;
        retVal = 1;
      }
    }
    double verifyTime2ndDeriv = timer.seconds();

    /*-------- End of 2nd Deriv Test --------*/
    // Outputting the metrics
    std::cout << "-------- Duration Measured --------" << std::endl;
    std::cout << "Large scale data generation time: " << dataGenTime
              << std::endl;
    std::cout << "Data copy time: " << dataCopyTime << std::endl;
    std::cout << "Serial spline initialization time: "
              << serialSplineCreationTime << std::endl;
    std::cout << "Kokkos spline initialization time: "
              << kokkosSplineCreationTime << std::endl;
    std::cout << "Kokkos csr initialization time: " << csrInitTime << std::endl;
    std::cout << "Serial 1st derivative time: " << serial1stDerivTime
              << std::endl;
    std::cout << "Kokkos 1st derivative time: " << kokkos1stDerivTime
              << std::endl;
    std::cout << "1st deriv verify time: " << verifyTime1stDeriv << std::endl;
    std::cout << "Serial 2nd derivative time: " << serial2ndDerivTime
              << std::endl;
    std::cout << "Kokkos 2nd derivative time: " << kokkos2ndDerivTime
              << std::endl;
    std::cout << "2nd deriv verify time: " << verifyTime2ndDeriv << std::endl;
  }
  Kokkos::finalize();
  return retVal;
}
