#include "BSplineKokkos2D.h"
#include "makeScalingTestData.h"
#include <Kokkos_Core.hpp>
#include <vector>

#include "curveReader.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>

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
    const double EPSILON = 1e-11;
    // Creating Scaling Data
    // Using Kokkos timer to keep track of time elapsed
    Kokkos::Timer timer;
    const int numSplines = std::atoi(argv[1]);
    const int ptsPerSpline = std::atoi(argv[2]);
    Kokkos::View<double *[2], MemSpace> pts("PointsOnCircle",
                                            numSplines * ptsPerSpline);
    pts = makeCircle(1000.0, 2000.0, 500.0, numSplines, ptsPerSpline);
    double dataGenTime = timer.seconds();

    // Copy the result to vectors for serial initialization
    auto ptsMirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pts);
    Kokkos::deep_copy(ptsMirror, pts);
    std::vector<double> ptsX(numSplines * ptsPerSpline),
        ptsY(numSplines * ptsPerSpline);
    for (int i = 0; i < numSplines * ptsPerSpline; i++) {
      ptsX[i] = ptsMirror(i, 0);
      ptsY[i] = ptsMirror(i, 1);
    }
    double dataCopyTime = timer.seconds();
    timer.reset();

    // Initializing Serial BSpline, divide the given points into correct number
    // of splines to initialize A vector that holds all the BSpline2d created
    std::vector<SplineInterp::BSpline2d> allSerial(numSplines);
    std::vector<double> subPtsX(ptsPerSpline), subPtsY(ptsPerSpline);
    for (int i = 0; i < numSplines; i++) {
      int start = i * ptsPerSpline;
      for (int j = start; j < start + ptsPerSpline; j++) {
        subPtsX[j - start] = ptsX[j];
        subPtsY[j - start] = ptsY[j];
      }
      if (ptsPerSpline == 2) {
        allSerial[i] =
            SplineInterp::attach_piecewise_linear_curve(subPtsX, subPtsY);
      } else {
        allSerial[i] = SplineInterp::fitCubicSplineToPoints(subPtsX, subPtsY);
      }
    }
    double serialSplineCreationTime = timer.seconds();

    // BSplineKokkos2D object initialization based on serial spline
    timer.reset();
    BSplineKokkos2D<ExecutionSpace> kokkosBSP(allSerial);
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
    // Each of the splines will evaluate at all the para coords
    std::vector<double> serialResX(evalAt.size() * numSplines);
    std::vector<double> serialResY(evalAt.size() * numSplines);
    int idx = 0;
    timer.reset();
    for (int i = 0; i < numSplines; i++) {
      for (int j = 0; j < evalAt.size(); j++) {
        serialResX[idx + j] = allSerial[i].x.evalFirstDeriv(evalAt[j]);
        serialResY[idx + j] = allSerial[i].y.evalFirstDeriv(evalAt[j]);
      }
      idx += paraCoords;
    }
    double serial1stDerivTime = timer.seconds();

    // Kokkos evaluation of 1st derivative
    const size_t paraSize = evalAt.size();

    timer.reset(); // Recording the time needed for CSR initialization
    BSplineKokkos2D<ExecutionSpace>::CSR kokkosBSPCSR(allSerial.size(),
                                                      paraSize);

    auto valsMirror = Kokkos::create_mirror_view(kokkosBSPCSR.paraCoor);
    for (int i = 0; i < evalAt.size(); i++) {
      valsMirror(i) = evalAt[i];
    }
    auto splineIdxMirror = Kokkos::create_mirror_view(kokkosBSPCSR.splineIdx);
    for (int i = 0; i < allSerial.size(); i++) {
      splineIdxMirror(i) = i;
    }
    auto offsetMirror = Kokkos::create_mirror_view(kokkosBSPCSR.offset);
    int curOffset = 0;
    for (int i = 0; i < allSerial.size(); i++) {
      offsetMirror(i) = curOffset;
      curOffset += paraCoords;
    }
    offsetMirror(offsetMirror.extent(0) - 1) = curOffset;
    Kokkos::deep_copy(kokkosBSPCSR.splineIdx, splineIdxMirror);
    Kokkos::deep_copy(kokkosBSPCSR.offset, offsetMirror);
    Kokkos::deep_copy(kokkosBSPCSR.paraCoor, valsMirror);
    double csrInitTime = timer.seconds();
    // Calling batch 1st derivative function
    Kokkos::View<double *[2], MemSpace> res("derivResult",
                                            evalAt.size() * numSplines);
    timer.reset();
    res = kokkosBSP.eval1stDeriv(kokkosBSPCSR);
    double kokkos1stDerivTime = timer.seconds();
    auto resMirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), res);

    // Compare the serial and kokkos result
    timer.reset();
    for (int i = 0; i < evalAt.size() * numSplines; i++) {
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
    idx = 0;
    timer.reset();
    for (int i = 0; i < numSplines; i++) {
      for (int j = 0; j < evalAt.size(); j++) {
        serialResX[idx + j] = allSerial[i].x.evalSecondDeriv(evalAt[j]);
        serialResY[idx + j] = allSerial[i].y.evalSecondDeriv(evalAt[j]);
      }
      idx += paraCoords;
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
