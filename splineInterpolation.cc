#include "modelGen2d.h"
#include "Omega_h_file.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_for.hpp" //parallel_for
#include <cassert>
#include <cmath> //sqrt
#include <vector>

using std::vector;

/* DGESV prototype */
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
                       double *b, int *ldb, int *info);

template<typename T>
auto ohHostWriteToRead(Omega_h::HostWrite<T> hw) {
  return Omega_h::read(Omega_h::Write(hw));
}

bool areKnotsIncreasing(Omega_h::LOs splineToKnots, Omega_h::Reals x, Omega_h::Reals y) {
  assert(x.size() == y.size());
  Omega_h::parallel_for(splineToKnots.size()-1, OMEGA_H_LAMBDA(Omega_h::LO i) {
    for (auto j = splineToKnots[i]; j < splineToKnots[i+1]-1; j++) {
      OMEGA_H_CHECK(x[j] <= x[j+1]);
      OMEGA_H_CHECK(y[j] <= y[j+1]);
    }
  });
  return true;
}

namespace SplineInterp {

void SplineInfo::writeSamplesToCsv(std::string filename) {
  std::ofstream file(filename);
  assert(file.is_open());
  file << "splineId, x, y\n";
  int id=0;
  for(auto& spline : splines) {
    auto numSamples = 4;
    auto step = 1.0/(numSamples-1);
    for(int i = 0; i < numSamples; ++i) {
      auto t = step * i;
      double x = spline.x.eval(t);
      double y = spline.y.eval(t);
      file << id << ", " << x << ", " << y << "\n";
    }
    id++;
  }
  file.close();
}

void SplineInfo::writeToOsh(std::string filename) {
    std::ofstream file(filename);
    assert(file.is_open());

    Omega_h::HostWrite<Omega_h::LO> numCtrlPts_h(splines.size());
    Omega_h::HostWrite<Omega_h::LO> numKnots_h(splines.size());
    Omega_h::HostWrite<Omega_h::LO> order_h(splines.size());
    int i=0;
    for(auto& spline : splines) {
      assert(spline.x.getOrder() == spline.y.getOrder());
      order_h[i] = spline.x.getOrder();
      assert(spline.x.getNumCtrlPts() == spline.y.getNumCtrlPts());
      assert(spline.x.getNumKnots() == spline.y.getNumKnots());
      //store the number of ctrlPts/knots for each spline then build the offset array
      numCtrlPts_h[i] = spline.x.getNumCtrlPts();
      numKnots_h[i] = spline.x.getNumKnots();
      i++;
    }
    auto splineToCtrlPts = Omega_h::offset_scan( ohHostWriteToRead(numCtrlPts_h) );
    auto splineToKnots = Omega_h::offset_scan( ohHostWriteToRead(numKnots_h) );

    //fill in the arrays of all ctrlPts and knots
    const auto totNumCtrlPts = splineToCtrlPts.last();
    const auto totNumKnots = splineToKnots.last();
    Omega_h::HostWrite<Omega_h::Real> ctrlPts_x(totNumCtrlPts);
    Omega_h::HostWrite<Omega_h::Real> ctrlPts_y(totNumCtrlPts);
    Omega_h::HostWrite<Omega_h::Real> knots_x(totNumKnots);
    Omega_h::HostWrite<Omega_h::Real> knots_y(totNumKnots);
    int cpIdx = 0;
    int knotIdx = 0;
    for(auto& spline : splines) {
      for(std::size_t j=0;  j<spline.x.getNumCtrlPts(); j++) {
        ctrlPts_x[cpIdx] = spline.x.getCtrlPt(j);
        ctrlPts_y[cpIdx] = spline.y.getCtrlPt(j);
        cpIdx++;
      }
      for(std::size_t j=0;  j<spline.x.getNumKnots(); j++) {
        knots_x[knotIdx] = spline.x.getKnot(j);
        knots_y[knotIdx] = spline.y.getKnot(j);
        knotIdx++;
      }
    }
    assert(totNumCtrlPts == cpIdx);
    assert(totNumKnots == knotIdx);

    areKnotsIncreasing(splineToKnots, ohHostWriteToRead(knots_x), ohHostWriteToRead(knots_y));

    const int compressed = 0;
    //the following is from src/Omega_h_file.cpp write(...)
    unsigned char const magic[2] = {0xa1, 0x1a};
    file.write(reinterpret_cast<const char*>(magic), sizeof(magic));
    bool needs_swapping = !Omega_h::is_little_endian_cpu();
    Omega_h::binary::write_value(file, compressed, needs_swapping);

    Omega_h::binary::write_array(file, splineToCtrlPts, compressed, needs_swapping);
    Omega_h::binary::write_array(file, splineToKnots, compressed, needs_swapping);
    Omega_h::binary::write_array(file, ohHostWriteToRead(ctrlPts_x), compressed, needs_swapping);
    Omega_h::binary::write_array(file, ohHostWriteToRead(ctrlPts_y), compressed, needs_swapping);
    Omega_h::binary::write_array(file, ohHostWriteToRead(knots_x), compressed, needs_swapping);
    Omega_h::binary::write_array(file, ohHostWriteToRead(knots_y), compressed, needs_swapping);
    Omega_h::binary::write_array(file, ohHostWriteToRead(order_h), compressed, needs_swapping);

    file.close();
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
BSpline2d attach_piecewise_linear_curve(std::vector<double> xpts, std::vector<double> ypts) {
  assert(xpts.size() == ypts.size());
  std::vector<double> pts;
  pts.reserve(xpts.size()*3);
  for(int i=0; i<xpts.size(); i++) {
    pts.push_back(xpts.at(i));
    pts.push_back(ypts.at(i));
    pts.push_back(0);
  }
  return attach_piecewise_linear_curve(pts);
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
