#include "splineInterpolation.h"
#include <cassert>
#include <vector>

using std::vector;

/* DGESV prototype */
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
                       double *b, int *ldb, int *info);

namespace SplineInterp {

bool curveOrientation(std::vector<double> &curvePts) {
  int numPts = curvePts.size() / 3;
  double sumEdges = 0.0;
  for (int i = 1; i < numPts; ++i) {
    double x1 = curvePts[(i - 1) * 3];
    double y1 = curvePts[(i - 1) * 3 + 1];
    double x2 = curvePts[i * 3];
    double y2 = curvePts[i * 3 + 1];

    double areaUnderEdge = (x2 - x1) * (y2 + y1);
    sumEdges = sumEdges + areaUnderEdge;
  }

  if (sumEdges > 0)
    return true;

  return false;
}

void interpolateCubicBSpline(vector<double> &points, vector<double> &knots,
                             vector<double> &ctrlPoints, int bc) {
  int numPts = points.size();
  int order_p = 4;
  assert(numPts > 1);
  ctrlPoints.resize(numPts + 2);
  vector<double> coeffs((numPts + 2) * (numPts + 2),
                        0.); // numPts+2 * numPts+2 linear system
  // first find the constraint to coninside with points
  // 2 natural cubic spline at the boundary
  for (int i = 0; i < numPts + 2; i++) {
    vector<double> points_tmp(numPts + 2, 0.);
    points_tmp.at(i) = 1.0;
    vector<double> weight_p;
    Spline::BSpline basis(order_p, points_tmp, knots, weight_p);
    for (int j = 0; j < numPts; j++) {
      double para = knots.at(order_p + j - 1);
      double res = basis.eval(para);
      double secondDeriv0 = basis.evalSecondDeriv(0);
      double secondDeriv1 = basis.evalSecondDeriv(1);
      coeffs.at(i * (numPts + 2) + j) = res;
    }
    double secondDeriv0 = basis.evalSecondDeriv(0);
    double secondDeriv1 = basis.evalSecondDeriv(1);
    if (bc == 0) // natural
    {
      coeffs.at(i * (numPts + 2) + numPts) = secondDeriv0;
      coeffs.at(i * (numPts + 2) + numPts + 1) = secondDeriv1;
    } else // periodic
    {
      double firstDeriv0 = basis.evalFirstDeriv(0);
      double firstDeriv1 = basis.evalFirstDeriv(1);
      coeffs.at(i * (numPts + 2) + numPts) = firstDeriv0 - firstDeriv1;
      coeffs.at(i * (numPts + 2) + numPts + 1) = secondDeriv0 - secondDeriv1;
    }
  }

  // set up the linear system and solve
  vector<double> rhs(numPts + 2, 0.0);
  for (int i = 0; i < numPts; i++)
    rhs.at(i) = points.at(i);

  int info, one = 1, dim = numPts + 2;
  vector<int> ipiv(dim, 0);
  dgesv_(&dim, &one, &(coeffs.at(0)), &dim, &(ipiv.at(0)), &(rhs.at(0)), &dim,
         &info);
  assert(info == 0);
  for (int i = 0; i < numPts + 2; i++)
    ctrlPoints.at(i) = rhs.at(i);
}

BSpline2d fitCubicSplineToPoints(std::vector<double> xpts,
                                 std::vector<double> ypts) {
  assert(xpts.size() == ypts.size());
  const auto numPts = xpts.size();
  int order_p = 4;
  int knotsize = 2 * order_p + numPts - 2;
  vector<double> knots(knotsize, 0.);
  vector<double> ctrlPointsX(numPts + 2), ctrlPointsY(numPts + 2), weight;
  for (int i = 0; i < order_p; i++) {
    knots.at(knotsize - i - 1) = 1.0;
  }
  double increment = 1.0 / (numPts - 1);
  for (int i = 0; i < numPts - 2; i++) {
    // double increment=inter_len.at(i)/len;
    knots.at(order_p + i) = knots.at(order_p + i - 1) + increment;
  }
  interpolateCubicBSpline(xpts, knots, ctrlPointsX, 0);
  interpolateCubicBSpline(ypts, knots, ctrlPointsY, 0);
  Spline::BSpline xSpline(order_p, ctrlPointsX, knots, weight);
  Spline::BSpline ySpline(order_p, ctrlPointsY, knots, weight);
  return {xSpline, ySpline};
}

} // namespace SplineInterp
