#include "modelGen2d.h"
#include "Omega_h_file.hpp"
#include <cassert>
#include <cmath> //sqrt
#include <vector>

using std::vector;

/* DGESV prototype */
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
                       double *b, int *ldb, int *info);

namespace SplineInterp {

void SplineInfo::writeToOsh(std::string filename) {
    std::ofstream file(filename);
    assert(file.is_open());
    Omega_h::HostWrite<Omega_h::LO> splineToCtrlPts_h(splines.size()+1, "splineToCtrlPts_h");
    Omega_h::HostWrite<Omega_h::LO> splineToKnots_h(splines.size()+1, "splineToKnots_h");
    Omega_h::HostWrite<Omega_h::LO> order_h(splines.size(), "order_h");
    int i=0;
    for(auto& spline : splines) {
      assert(spline.x.getOrder() == spline.y.getOrder());
      order_h[i] = spline.x.getOrder();
    }
    bool compressed = true;
    bool needs_swapping = false;
    //Omega_h::write_array(file, splineToCtrlPts);
    //Omega_h::write_array(file, ctrlPts);
    //Omega_h::write_array(file, splineToKnots);
    //Omega_h::write_array(file, knots);
    Omega_h::binary::write_array(file, Omega_h::read(Omega_h::Write(order_h)), compressed, needs_swapping);
  }


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

BSpline2d attach_piecewise_linear_curve(std::vector<double> points) {
  const int dim = 3;
  assert(points.size() % dim == 0);
  const auto numPts = points.size() / dim;
  const int order_p=2;
  const int knotsize=2*order_p+numPts-2;
  vector<double> knots(knotsize,0.);
  vector<double> ctrlPointsX(numPts),ctrlPointsY(numPts),weight;
  for( int i=0; i<order_p; i++) {
    knots.at(knotsize-i-1)=1.0;
  }
  double totalLength = 0.0; 
  std::vector <double> lengthVector;
  for( int i=1; i<numPts; i++) {
    double l = getLengthSquared(points.at(dim*i-2), points.at(dim*i-1),
                                points.at(dim*i), points.at(dim*i+1));
    double length = std::sqrt(l);
    totalLength += length;
    lengthVector.push_back(totalLength);
  }
  for (int i=0; i<numPts-2; i++) {
    double par = lengthVector[i]/totalLength;
    knots.at(order_p+i)=par;
  }
  for( int i=0; i<numPts; i++) {
    ctrlPointsX.at(i)=points[dim*i];
    ctrlPointsY.at(i)=points[dim*i+1];
  }
  Spline::BSpline xSpline(order_p, ctrlPointsX, knots, weight);
  Spline::BSpline ySpline(order_p, ctrlPointsY, knots, weight);
  return {xSpline, ySpline};
}


void interpolateCubicBSpline(std::vector<double> &points, vector<double> &knots,
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

BSpline2d fitCubicSplineToPoints(std::vector<double> pts) {
  assert(pts.size() % 3 == 0);
  const auto numPts = pts.size() / 3;
  //TODO replace with mdspan
  std::vector<double> xpts;
  xpts.reserve(numPts);
  std::vector<double> ypts;
  ypts.reserve(numPts);
  for(int i=0; i<numPts; i++) {
     xpts.push_back(pts.at(i*3));
     ypts.push_back(pts.at(i*3+1));
  }
  return fitCubicSplineToPoints(xpts,ypts);
}

} // namespace SplineInterp
