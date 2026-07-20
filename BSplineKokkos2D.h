#ifndef BSPLINEKOKKOS2D_H
#define BSPLINEKOKKOS2D_H

#include<Kokkos_Core.hpp>
#include<vector>
#include<iostream>
#include<string>
template<typename ExecutionSpace>

class BSplineKokkos2D {
public:
  using MemSpace = typename ExecutionSpace::memory_space;

  BSplineKokkos2D(int order_p, std::vector<double>& ctrlPts_x, std::vector<double>& ctrlPts_y, std::vector<double>& knotsI) {
    Kokkos::View<int*, MemSpace> orderV("Orders", 1);
    auto host_orderV = Kokkos::create_mirror_view(orderV);
    host_orderV(0) = order_p;
    order = orderV;
    Kokkos::deep_copy(order, host_orderV);
  
    Kokkos::View<double*[2], MemSpace> ctrlPtsV("ctrlPts", ctrlPts_x.size());
    auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
    for (int i = 0; i < ctrlPts_x.size(); i++) {
      host_ctrlPtsV(i, 0) = ctrlPts_x[i];
      host_ctrlPtsV(i, 1) = ctrlPts_y[i];
    }
    ctrlPts = ctrlPtsV;
    Kokkos::deep_copy(ctrlPts, host_ctrlPtsV);

    Kokkos::View<int*, MemSpace> cpOffsetV("cpOffset", 1);
    auto host_cpOffsetV = Kokkos::create_mirror_view(cpOffsetV);
    host_cpOffsetV(0) = ctrlPts_x.size();
    host_cpOffsetV(0) = ctrlPts_x.size();

    cpOffset = cpOffsetV;
    Kokkos::deep_copy(cpOffset, host_cpOffsetV);

    Kokkos::View<double*, MemSpace> knotsV("knots", knotsI.size());
    auto host_knotsV = Kokkos::create_mirror_view(knotsV);
    for (int i = 0; i < knotsI.size(); i++) {
      host_knotsV(i) = knotsI[i];
    }
    knots = knotsV;
    Kokkos::deep_copy(knots, host_knotsV);
    
    Kokkos::View<int*, MemSpace> knotsOffsetV("knotsOffset", 1);
    auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
    host_knotsOffsetV(0) = knotsI.size();
    knotsOffset = knotsOffsetV;
    Kokkos::deep_copy(knotsOffset, host_knotsOffsetV);

    calculateDerivCoeff();
  }

  void calculateDerivCoeff() {
    //Allocate the views we need
    Kokkos::View<double*[2], MemSpace> ctrlPts1stDV("ctrlPts1stDeriv", ctrlPts.extent(0)-cpOffset.extent(0));
    Kokkos::View<int*, MemSpace> cP1stDOffsetV("cP1stDOffset", cpOffset.extent(0));
    Kokkos::View<int*, MemSpace> cP2ndDOffsetV("cP2ndDOffset", cpOffset.extent(0));

    //Set up the offset views
    auto mvCP1stDOffsetV = Kokkos::create_mirror_view(cP1stDOffsetV);
    auto mvCP2ndDOffsetV = Kokkos::create_mirror_view(cP2ndDOffsetV);

    auto mvCPOffset = Kokkos::create_mirror_view(cpOffset);
    
    if (mvCPOffset.extent(0) > 1) {
      Kokkos::deep_copy(mvCPOffset, cpOffset);
      for (int i = 1; i < mvCPOffset.extent(0); i++) {
        mvCP1stDOffsetV(i) = mvCPOffset(i)-1;
        mvCP2ndDOffsetV(i) = mvCPOffset(i)-2;
      }
    }
    else {
      mvCP1stDOffsetV(0) = mvCPOffset(0)-1;
      mvCP2ndDOffsetV(0) = mvCPOffset(0)-2;
    }

    auto mvOrder = Kokkos::create_mirror_view(order);
    auto mvKnots = Kokkos::create_mirror_view(knots);
    auto mvCtrlPts = Kokkos::create_mirror_view(ctrlPts);
    auto mvCtrlPts1stDV = Kokkos::create_mirror_view(ctrlPts1stDV);
    Kokkos::deep_copy(mvOrder, order);
    Kokkos::deep_copy(mvKnots, knots);
    Kokkos::deep_copy(mvCtrlPts, ctrlPts);

    //Calculate 1st derivative coef
    int offidx = 0;  //Offset index
    int oidx = 0;  //Order index
    for (int i = 1; i < mvCtrlPts1stDV.extent(0)+1; i++) {
      if (i == mvCP1stDOffsetV(offidx)+1) {
        oidx++;
        offidx++;
        continue;
      }
      //We need to check whether we are on the border for the next spline in our structure
      double delta = double(mvOrder(oidx)-1) / (mvKnots(i+mvOrder(oidx)-1) - mvKnots(i));

      mvCtrlPts1stDV(i-1, 0) = (mvCtrlPts(i, 0) - mvCtrlPts(i-1, 0))*delta;
      mvCtrlPts1stDV(i-1, 1) = (mvCtrlPts(i, 1) - mvCtrlPts(i-1, 1))*delta;
    }

    //Calculate 2nd derivative coef
    Kokkos::View<double*[2], MemSpace> ctrlPts2ndDV("ctrlPts2ndDeriv", ctrlPts.extent(0)-(2*cpOffset.extent(0)));
    auto mvCtrlPts2ndDV = Kokkos::create_mirror_view(ctrlPts2ndDV);

    offidx = 0;
    oidx = 0;
    for (int i = 1; i < ctrlPts1stDV.extent(0); i++) {
      if (i == mvCP2ndDOffsetV(offidx)+2) {
        oidx++;
        offidx++;
        continue;
      }
      double delta = double(mvOrder(oidx)-2) / (mvKnots(i+mvOrder(oidx)-1) - mvKnots(i+1));

      mvCtrlPts2ndDV(i-1, 0) = (mvCtrlPts1stDV(i, 0) - mvCtrlPts1stDV(i-1, 0)) * delta;
      mvCtrlPts2ndDV(i-1, 1) = (mvCtrlPts1stDV(i, 1) - mvCtrlPts1stDV(i-1, 1)) * delta;
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

  KOKKOS_FUNCTION void eval1stDerivDeBoor(double x, int splineIdx, Kokkos::View<double*, MemSpace> result, const Kokkos::View<int*, MemSpace> order, const Kokkos::View<double*, MemSpace> knots, const Kokkos::View<double*[2], MemSpace> ctrlPts_1stD) const {
    //DeBoor's algorithm for BSpline 1st deriv calculation
    const int MAX_DEGREE = 3;
    int lKnot = order(splineIdx);
    lKnot--;
    int resultOrder = lKnot;
    int leftPt = 0;

    while (x > knots(lKnot+1)) {
      lKnot++;
      leftPt++;
    }

    //Allocate temporary points
    double ptsX[MAX_DEGREE+1];
    double ptsY[MAX_DEGREE+1];

    int idx = 0;
    for (int i = leftPt; i < leftPt + resultOrder; i++) {
      ptsX[idx] = ctrlPts_1stD(i, 0);
      ptsY[idx] = ctrlPts_1stD(i, 1);
      idx++;
    }

    auto localKnots = Kokkos::subview(knots, Kokkos::pair<int, int>(lKnot-resultOrder+2, lKnot+resultOrder));

    //Calculation loop
    for (int r = 1; r <= resultOrder; r++) {
      for(int i = resultOrder-1; i >= r; i--) {
        double alpha;
        if (localKnots(i-1) - localKnots(i+resultOrder-r-1) > 1e-12) {
          alpha = 0;
        }
        else {
          alpha = (x - localKnots(i-1)) / (knots(i+resultOrder-r-1) - knots(i-1));
        }
        ptsX[i] = (1. - alpha) * ptsX[i-1]+alpha * ptsX[i];
        ptsY[i] = (1. - alpha) * ptsY[i-1]+alpha * ptsY[i];
      }

      result(0) = ptsX[resultOrder-1];
      result(1) = ptsY[resultOrder-1];
    }
  }

  Kokkos::View<double*,MemSpace> eval1stDeriv(Kokkos::View<double*, MemSpace> xVals, int splineIdx) const {
    Kokkos::View<double*, MemSpace> res("result", 2);
    Kokkos::parallel_for("parallel De Boor's for 1st derivative", xVals.size(), KOKKOS_CLASS_LAMBDA(int i){
      eval1stDerivDeBoor(xVals(i), splineIdx, res, order, knots, ctrlPts1stD);
    });
    return res;
  }

  //Second derivative calculation
  KOKKOS_FUNCTION void eval2ndDerivDeBoor(double x, int splineIdx, Kokkos::View<double*, MemSpace> result, const Kokkos::View<int*, MemSpace> order, const Kokkos::View<double*, MemSpace> knots, const Kokkos::View<double*[2], MemSpace> ctrlPts_2ndD) const {
    if (order(splineIdx) == 2) {
      result(0) = 0;
      result(1) = 0;
      return;
    }
    const int MAX_DEGREE = 3;
    int lKnot = order(splineIdx);
    lKnot--;

    int resultOrder = lKnot-1;
    int leftPt = 0;

    while (knots(lKnot+1) < x) {
      lKnot++;
      leftPt++;
    }

    double ptsX[MAX_DEGREE+1];
    double ptsY[MAX_DEGREE+1];

    int idx = 0;
    for (int i = leftPt; i < leftPt+resultOrder; i++) {
      ptsX[idx] = ctrlPts_2ndD(i, 0);
      ptsY[idx] = ctrlPts_2ndD(i, 1);
      idx++;
    }

    auto localKnots = Kokkos::subview(knots, Kokkos::pair<int, int>(lKnot-resultOrder+2, lKnot+resultOrder));

    for (int r = 1; r <= resultOrder; r++) {
      for (int i = resultOrder-1; i >= r; i--) {
        double aLeft = localKnots(i-1);
        double aRight = localKnots(i+resultOrder-r-1);
        double alpha;
        if (aLeft == aRight) {
          alpha = 0.;
        }
        else {
          alpha = (x - aLeft) / (aRight - aLeft);
        }
        ptsX[i] = (1. - alpha) * ptsX[i-1] + alpha * ptsX[i];
        ptsY[i] = (1. - alpha) * ptsY[i-1] + alpha * ptsY[i];
      }
    }

    result(0) = ptsX[resultOrder-1];
    result(1) = ptsY[resultOrder-1];

  }

  Kokkos::View<double*, MemSpace> eval2ndDeriv(Kokkos::View<double*, MemSpace> xVals, int splineIdx) const {
    Kokkos::View<double*, MemSpace> res("result", 2);
    Kokkos::parallel_for("parallel De Boor's for 2nd derivative", xVals.size(), KOKKOS_CLASS_LAMBDA(int i) {
      eval2ndDerivDeBoor(xVals(i), splineIdx, res, order, knots, ctrlPts2ndD); 
    });
    return res;
  }

  //Accessors
  Kokkos::View<int*, MemSpace> getOrder() const {return order;}
  Kokkos::View<double*[2], MemSpace> getCtrlPts() const {return ctrlPts;}
  Kokkos::View<int*, MemSpace> getCPOffset() const {return cpOffset;}
  Kokkos::View<double*, MemSpace> getKnots() const {return knots;}
  Kokkos::View<int*, MemSpace> getKnotsOffset() const {return knotsOffset;}
  Kokkos::View<double*[2], MemSpace> getCP1stD() const {return ctrlPts1stD;}
  Kokkos::View<int*, MemSpace> getCP1stDOffset() const {return cp1stDOffset;}
  Kokkos::View<double*[2], MemSpace> getCP2ndD() const {return ctrlPts2ndD;}
  Kokkos::View<int*, MemSpace> getCP2ndDOffset() const {return cp2ndDOffset;}


private:
  Kokkos::View<int*, MemSpace> order;
  Kokkos::View<double*[2], MemSpace> ctrlPts;
  Kokkos::View<int*, MemSpace> cpOffset;
  Kokkos::View<double*, MemSpace> knots;
  Kokkos::View<int*, MemSpace> knotsOffset;
  Kokkos::View<double*[2], MemSpace> ctrlPts1stD;
  Kokkos::View<double*[2], MemSpace> ctrlPts2ndD;
  Kokkos::View<int*, MemSpace> cp1stDOffset;
  Kokkos::View<int*, MemSpace> cp2ndDOffset;

};

#endif
