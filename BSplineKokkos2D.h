#ifndef BSPLINEKOKKOS2D_H
#define BSPLINEKOKKOS2D_H

#include <Kokkos_Core.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "splineInterpolation.h"

const int MAX_DEGREE = 3;

template <typename ExecutionSpace> class BSplineKokkos2D {
public:
  using MemSpace = typename ExecutionSpace::memory_space;

  // Custom CSR
  struct CSR {
    CSR(const size_t splineIdxSize, const size_t paraCoorSize)
        : splineIdx("splineIdx", splineIdxSize),
          offset("offset", splineIdxSize + 1),
          paraCoor("paraCoor", paraCoorSize){};
    Kokkos::View<int *, MemSpace> splineIdx;
    Kokkos::View<int *, MemSpace> offset;
    Kokkos::View<double *, MemSpace> paraCoor;
  };

  BSplineKokkos2D(int order_p, std::vector<double> &ctrlPts_x,
                  std::vector<double> &ctrlPts_y, std::vector<double> &knotsI) {
    Kokkos::View<int *, MemSpace> orderV("Orders", 1);
    auto host_orderV = Kokkos::create_mirror_view(orderV);
    host_orderV(0) = order_p;
    order = orderV;
    Kokkos::deep_copy(order, host_orderV);

    Kokkos::View<double *[2], MemSpace> ctrlPtsV("ctrlPts", ctrlPts_x.size());
    auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
    for (int i = 0; i < ctrlPts_x.size(); i++) {
      host_ctrlPtsV(i, 0) = ctrlPts_x[i];
      host_ctrlPtsV(i, 1) = ctrlPts_y[i];
    }
    ctrlPts = ctrlPtsV;
    Kokkos::deep_copy(ctrlPts, host_ctrlPtsV);

    Kokkos::View<int *, MemSpace> cpOffsetV("cpOffset", 2);
    auto host_cpOffsetV = Kokkos::create_mirror_view(cpOffsetV);
    host_cpOffsetV(0) = 0;
    host_cpOffsetV(1) = ctrlPts_x.size();

    cpOffset = cpOffsetV;
    Kokkos::deep_copy(cpOffset, host_cpOffsetV);

    Kokkos::View<double *, MemSpace> knotsV("knots", knotsI.size());
    auto host_knotsV = Kokkos::create_mirror_view(knotsV);
    for (int i = 0; i < knotsI.size(); i++) {
      host_knotsV(i) = knotsI[i];
    }
    knots = knotsV;
    Kokkos::deep_copy(knots, host_knotsV);

    Kokkos::View<int *, MemSpace> knotsOffsetV("knotsOffset", 2);
    auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
    host_knotsOffsetV(0) = 0;
    host_knotsOffsetV(1) = knotsI.size();
    knotsOffset = knotsOffsetV;
    Kokkos::deep_copy(knotsOffset, host_knotsOffsetV);

    calculateDerivCoeff();
  }

  BSplineKokkos2D(std::vector<SplineInterp::BSpline2d> serialBSP) {
    //Initializing order view
    Kokkos::View<int*, MemSpace> orderV("Orders", serialBSP.size());
    auto host_orderV = Kokkos::create_mirror_view(orderV);
    //Initializing ctrlPtsOffset
    Kokkos::View<int*, MemSpace> cpOffsetV("cpOffset", serialBSP.size()+1);
    auto host_cpOffsetV = Kokkos::create_mirror_view(cpOffsetV);
    host_cpOffsetV(0) = 0;
    //Initializing knotsOffset
    Kokkos::View<int*, MemSpace> knotsOffsetV("knotsOffset", serialBSP.size()+1);
    auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
    host_knotsOffsetV(0) = 0;
    
    //Obtain order, ctrlPtsOffset, knotsOffset values
    int splineOrder;
    std::vector<double> splineX, splineY, splineKnots, weights;
    for (int i = 0; i < serialBSP.size(); i++) {
      serialBSP[i].x.getpara(splineOrder, splineX, splineKnots, weights);
      host_orderV(i) = splineOrder;
      host_cpOffsetV(i+1) = host_cpOffsetV(i) + splineX.size();
      host_knotsOffsetV(i+1) = host_knotsOffsetV(i) + splineKnots.size();
    }
    order = orderV;
    Kokkos::deep_copy(order, host_orderV);
    cpOffset = cpOffsetV;
    Kokkos::deep_copy(cpOffset, host_cpOffsetV);
    knotsOffset = knotsOffsetV;
    Kokkos::deep_copy(knotsOffset, host_knotsOffsetV);

    //Initializing ctrlPts and knots
    Kokkos::View<double*[2], MemSpace> ctrlPtsV("ctrlPts", host_cpOffsetV(host_cpOffsetV.extent(0)-1));
    auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
    Kokkos::View<double*, MemSpace> knotsV("knots", host_knotsOffsetV(host_knotsOffsetV.extent(0)-1));
    auto host_knotsV = Kokkos::create_mirror_view(knotsV);

    for (int i = 0; i < serialBSP.size(); i++) {
      serialBSP[i].x.getpara(splineOrder, splineX, splineKnots, weights);
      serialBSP[i].y.getpara(splineOrder, splineY, splineKnots, weights);
      for (int j = 0; j < splineKnots.size(); j++) {
        host_knotsV(host_knotsOffsetV(i) + j) = splineKnots[j];
      }
      for (int j = 0; j < splineX.size(); j++) {
        host_ctrlPtsV(host_cpOffsetV(i) + j, 0) = splineX[j];
        host_ctrlPtsV(host_cpOffsetV(i) + j, 1) = splineY[j];
      }
    }
    //Copy result to device
    ctrlPts = ctrlPtsV;
    Kokkos::deep_copy(ctrlPts, host_ctrlPtsV);
    knots = knotsV;
    Kokkos::deep_copy(knots, host_knotsV); 
    
    //Calculating derivative coefficient
    calculateDerivCoeff();
  }

  void calculateDerivCoeff() {
    // Allocate the views we need
    Kokkos::View<double *[2], MemSpace> ctrlPts1stDV(
        "ctrlPts1stDeriv", ctrlPts.extent(0) - (cpOffset.extent(0)-1));
    Kokkos::View<int *, MemSpace> cP1stDOffsetV("cP1stDOffset",
                                                cpOffset.extent(0));
    Kokkos::View<int *, MemSpace> cP2ndDOffsetV("cP2ndDOffset",
                                                cpOffset.extent(0));

    // Set up the offset views
    auto mvCP1stDOffsetV = Kokkos::create_mirror_view(cP1stDOffsetV);
    auto mvCP2ndDOffsetV = Kokkos::create_mirror_view(cP2ndDOffsetV);
    mvCP1stDOffsetV(0) = 0;
    mvCP2ndDOffsetV(0) = 0;

    auto mvCPOffset = Kokkos::create_mirror_view(cpOffset);
    Kokkos::deep_copy(mvCPOffset, cpOffset);

    for (int i = 1; i < mvCPOffset.extent(0); i++) {
      mvCP1stDOffsetV(i) = mvCPOffset(i) - 1;
      mvCP2ndDOffsetV(i) = mvCPOffset(i) - 2;
    }
    
    auto mvOrder = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), order);
    auto mvKnots = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), knots);
    auto mvCtrlPts = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ctrlPts);
    auto mvCtrlPts1stDV = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ctrlPts1stDV);
    auto mvKnotsOffset = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), knotsOffset);

    // Calculate 1st derivative coef
    // TO DELETE: Knots view is initialized correctly, indexing is off
    int offidx = 1; // Offset index
    int oidx = 0;   // Order index
    int kOffset = 0;   // KnotsOffset index
    int kIdx = 1;
    int ctrlOffset = 0;
    for (int i = 1; i < mvCtrlPts1stDV.extent(0) + 1; i++) {
      //Check whether we are on border for next spline
      if (i+oidx == mvCP1stDOffsetV(offidx)+1) {
        ctrlOffset = mvCP1stDOffsetV(offidx)+1;
        oidx++;
        offidx++;
        kIdx = 1;
        kOffset = mvKnotsOffset(oidx);
      }
      
      double delta = double(mvOrder(oidx) - 1) /
                     (mvKnots(kOffset + kIdx + mvOrder(oidx) - 1) - mvKnots(kOffset + kIdx));
      mvCtrlPts1stDV(i - 1, 0) =
          (mvCtrlPts(ctrlOffset + kIdx, 0) - mvCtrlPts(ctrlOffset + kIdx - 1, 0)) * delta;
      mvCtrlPts1stDV(i - 1, 1) =
          (mvCtrlPts(ctrlOffset + kIdx, 1) - mvCtrlPts(ctrlOffset + kIdx - 1, 1)) * delta;
      kIdx++;
    }
    //TO DELETE: 1stDerivative coeff calculation corrected

    // Calculate 2nd derivative coef
    Kokkos::View<double *[2], MemSpace> ctrlPts2ndDV(
        "ctrlPts2ndDeriv", ctrlPts.extent(0) - (2 * (cpOffset.extent(0)-1)));
    auto mvCtrlPts2ndDV = Kokkos::create_mirror_view(ctrlPts2ndDV);

    offidx = 1;
    oidx = 0;
    for (int i = 1; i < ctrlPts1stDV.extent(0); i++) {
      if (i == mvCP2ndDOffsetV(offidx) + 2) {
        oidx++;
        offidx++;
        continue;
      }
      double delta = double(mvOrder(oidx) - 2) /
                     (mvKnots(i + mvOrder(oidx) - 1) - mvKnots(i + 1));

      mvCtrlPts2ndDV(i - 1, 0) =
          (mvCtrlPts1stDV(i, 0) - mvCtrlPts1stDV(i - 1, 0)) * delta;
      mvCtrlPts2ndDV(i - 1, 1) =
          (mvCtrlPts1stDV(i, 1) - mvCtrlPts1stDV(i - 1, 1)) * delta;
    }
    Kokkos::deep_copy(cP1stDOffsetV, mvCP1stDOffsetV);
    Kokkos::deep_copy(cP2ndDOffsetV, mvCP2ndDOffsetV);
    Kokkos::deep_copy(ctrlPts1stDV, mvCtrlPts1stDV);
    Kokkos::deep_copy(ctrlPts2ndDV, mvCtrlPts2ndDV);

    cp1stDOffset = cP1stDOffsetV;
    cp2ndDOffset = cP2ndDOffsetV;
    ctrlPts1stD = ctrlPts1stDV;
    ctrlPts2ndD = ctrlPts2ndDV;
  }

  KOKKOS_FUNCTION void eval1stDerivDeBoor(
      const double x, int offset, Kokkos::View<double *[2], MemSpace> result,
      const int order, const Kokkos::View<double *, MemSpace> knots,
      const Kokkos::View<double *[2], MemSpace> ctrlPts_1stD) const {
    // DeBoor's algorithm for BSpline 1st deriv calculation
    int lKnot = order;
    lKnot--;
    int resultOrder = order - 1;
    int leftPt = 0;

    while (x > knots(lKnot + 1)) {
      lKnot++;
      leftPt++;
    }

    // Allocate temporary points(MAX_DEGREE+1 is sufficient)
    double ptsX[MAX_DEGREE + 1];
    double ptsY[MAX_DEGREE + 1];

    int idx = 0;
    for (int i = leftPt; i < leftPt + resultOrder; i++) {
      ptsX[idx] = ctrlPts_1stD(i, 0);
      ptsY[idx] = ctrlPts_1stD(i, 1);
      idx++;
    }

    auto localKnots =
        Kokkos::subview(knots, Kokkos::pair<int, int>(lKnot - resultOrder + 2,
                                                      lKnot + resultOrder));

    // Calculation loop
    for (int r = 1; r <= resultOrder; r++) {
      for (int i = resultOrder - 1; i >= r; i--) {
        double alpha;
        if (localKnots(i + resultOrder - r - 1) - localKnots(i - 1) < 1e-12) {
          alpha = 0.;
        } else {
          alpha = (x - localKnots(i - 1)) /
                  (localKnots(i + resultOrder - r - 1) - localKnots(i - 1));
        }
        ptsX[i] = (1. - alpha) * ptsX[i - 1] + alpha * ptsX[i];
        ptsY[i] = (1. - alpha) * ptsY[i - 1] + alpha * ptsY[i];
      }
    }
    result(offset, 0) = ptsX[resultOrder - 1];
    result(offset, 1) = ptsY[resultOrder - 1];
  }

  Kokkos::View<double *[2], MemSpace>
  eval1stDeriv(const CSR &derivSplines) const {
    // Create the result view based on the offset in the derivCSR
    int resultSize;
    Kokkos::deep_copy(resultSize,
                      Kokkos::subview(derivSplines.offset,
                                      derivSplines.offset.extent(0) - 1));
    Kokkos::View<double *[2], MemSpace> res("result", resultSize);
    int offsetSize = derivSplines.offset.extent(0);
    Kokkos::parallel_for(
        "parallel De Boor's for 1st derivative", offsetSize - 1,
        KOKKOS_CLASS_LAMBDA(const int i) {
          int knotsStart = knotsOffset(i);
          int knotsEnd = knotsOffset(i+1);
          int cp1Start = cp1stDOffset(i);
          int cp1End = cp1stDOffset(i+1);
          //std::cout << "knotsStart: " << knotsStart << " knotsEnd: " << knotsEnd << std::endl;
          //std::cout << "cp1Start: " << cp1Start << " cp1End: " << cp1End << std::endl;
          for (int j = derivSplines.offset(i); j < derivSplines.offset(i + 1);
               j++) {
            auto x = Kokkos::subview(derivSplines.paraCoor, j - derivSplines.offset(i));
            auto splineIdx = Kokkos::subview(derivSplines.splineIdx, i);
            auto splineOrder = Kokkos::subview(order, splineIdx());
            eval1stDerivDeBoor(x(), j, res, splineOrder(), Kokkos::subview(knots, std::pair<int, int>(knotsStart, knotsEnd)), Kokkos::subview(ctrlPts1stD, std::pair<int, int>(cp1Start, cp1End), Kokkos::ALL()));
          }
        });
    return res;
  }

  // Second derivative calculation
  KOKKOS_FUNCTION void eval2ndDerivDeBoor(
      const double x, int offset, Kokkos::View<double *[2], MemSpace> result,
      const int order, const Kokkos::View<double *, MemSpace> knots,
      const Kokkos::View<double *[2], MemSpace> ctrlPts_2ndD) const {
    if (order == 2) {
      result(offset, 0) = 0;
      result(offset, 1) = 0;
      return;
    }
    int lKnot = order - 1;
    int resultOrder = order - 2;
    int leftPt = 0;
    while (knots(lKnot + 1) < x) {
      lKnot++;
      leftPt++;
    }

    double ptsX[MAX_DEGREE + 1];
    double ptsY[MAX_DEGREE + 1];

    int idx = 0;
    for (int i = leftPt; i < leftPt + resultOrder; i++) {
      ptsX[idx] = ctrlPts_2ndD(i, 0);
      ptsY[idx] = ctrlPts_2ndD(i, 1);
      idx++;
    }

    auto localKnots =
        Kokkos::subview(knots, Kokkos::pair<int, int>(lKnot - resultOrder + 2,
                                                      lKnot + resultOrder));

    for (int r = 1; r <= resultOrder; r++) {
      for (int i = resultOrder - 1; i >= r; i--) {
        double aLeft = localKnots(i - 1);
        double aRight = localKnots(i + resultOrder - r - 1);
        double alpha;
        if (aRight - aLeft < 1e-12) {
          alpha = 0.;
        } else {
          alpha = (x - aLeft) / (aRight - aLeft);
        }
        ptsX[i] = (1. - alpha) * ptsX[i - 1] + alpha * ptsX[i];
        ptsY[i] = (1. - alpha) * ptsY[i - 1] + alpha * ptsY[i];
      }
    }
    result(offset, 0) = ptsX[resultOrder - 1];
    result(offset, 1) = ptsY[resultOrder - 1];
  }

  Kokkos::View<double *[2], MemSpace>
  eval2ndDeriv(const CSR &derivSplines) const {
    int resultSize;
    Kokkos::deep_copy(resultSize,
                      Kokkos::subview(derivSplines.offset,
                                      derivSplines.offset.extent(0) - 1));
    Kokkos::View<double *[2], MemSpace> res("result", resultSize);
    int offsetSize = derivSplines.offset.extent(0);
    Kokkos::parallel_for(
        "parallel De Boor's for 2nd derivative", offsetSize - 1,
        KOKKOS_CLASS_LAMBDA(const int i) {
          for (int j = derivSplines.offset(i); j < derivSplines.offset(i + 1);
               j++) {
            auto x = Kokkos::subview(derivSplines.paraCoor, j);
            auto splineIdxSub = Kokkos::subview(derivSplines.splineIdx, i);
            auto splineOrder = Kokkos::subview(order, splineIdxSub());
            eval2ndDerivDeBoor(x(), j, res, splineOrder(), knots, ctrlPts2ndD);
          }
        });
    return res;
  }

  // Accessors
  Kokkos::View<int *, MemSpace> getOrder() const {
    return order;
  }
  Kokkos::View<double *[2], MemSpace> getCtrlPts() const {
    return ctrlPts;
  } Kokkos::View<int *, MemSpace> getCPOffset() const {
    return cpOffset;
  }
  Kokkos::View<double *, MemSpace> getKnots() const { return knots; }
  Kokkos::View<int *, MemSpace> getKnotsOffset() const { return knotsOffset; }
  Kokkos::View<double *[2], MemSpace> getCP1stD() const {
    return ctrlPts1stD;
  } Kokkos::View<int *, MemSpace> getCP1stDOffset() const {
    return cp1stDOffset;
  }
  Kokkos::View<double *[2], MemSpace> getCP2ndD() const {
    return ctrlPts2ndD;
  } Kokkos::View<int *, MemSpace> getCP2ndDOffset() const {
    return cp2ndDOffset;
  }

private:
  Kokkos::View<int *, MemSpace> order;
  Kokkos::View<double *[2], MemSpace> ctrlPts;
  Kokkos::View<int *, MemSpace> cpOffset;
  Kokkos::View<double *, MemSpace> knots;
  Kokkos::View<int *, MemSpace> knotsOffset;
  Kokkos::View<double *[2], MemSpace> ctrlPts1stD;
  Kokkos::View<double *[2], MemSpace> ctrlPts2ndD;
  Kokkos::View<int *, MemSpace> cp1stDOffset;
  Kokkos::View<int *, MemSpace> cp2ndDOffset;
};

#endif
