/******************************************************************************

  (c) 2005-2016 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "BSpline.h"
#include <assert.h>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;
using std::vector;
using namespace Spline;

namespace {
  double norm(double x, double y) {
    return std::sqrt(x * x + y * y);
  }
}

BSpline::BSpline(int order_p, vector<double> &ctrlPts_p,
                 vector<double> &knots_p, vector<double> &weight_p) {
  assert(order_p > 1);
  order = order_p;
  ctrlPts = ctrlPts_p;
  knots = knots_p;
  // assert(knots_p.size()-ctrlPts_p.size()==order_p);
  weight = weight_p;
  calcuDerivCoeff();
}
void BSpline::calcuDerivCoeff() {
  // caculate coeffs of first deriv
  for (int i = 1; i < ctrlPts.size(); i++) {
    double delta =
        double((order - 1)) / (knots.at(i + order - 1) - knots.at(i));
    ctrlPts_1st.push_back((ctrlPts.at(i) - ctrlPts.at(i - 1)) * delta);
  }
  assert(ctrlPts_1st.size() == ctrlPts.size() - 1);
  // caculate coeffs of second deriv
  for (int i = 1; i < ctrlPts_1st.size(); i++) {
    double delta =
        double((order - 2)) / (knots.at(i + order - 1) - knots.at(i + 1));
    ctrlPts_2nd.push_back((ctrlPts_1st.at(i) - ctrlPts_1st.at(i - 1)) * delta);
  }
  if (order > 1)
    assert(ctrlPts_2nd.size() == ctrlPts.size() - 2);
}
double BSpline::eval(double x, bool debug) const {
  // first find the interval of x in knots
  int leftKnot = order - 1;
  int leftPt = 0;
  while (knots.at(leftKnot + 1) < x) {
    leftKnot++;
    leftPt++;
    if (leftKnot == knots.size() - 1)
      break;
  }

  if(debug) {
    for(int i=0; i<knots.size(); i++)
      printf("knots(%d) %f\n", i, knots[i]);
    for(int i=0; i<ctrlPts.size(); i++)
      printf("ctrlPts(%d) %f\n", i, ctrlPts[i]);
    printf("order %d coord %.2f leftPt %d leftKnot %d\n",
        order, x, leftPt, leftKnot);
  }
  vector<double> pts(&(ctrlPts[leftPt]), &ctrlPts[leftPt + order]);
  vector<double> localKnots(&(knots[leftKnot - order + 2]),
                            &(knots[leftKnot + order]));
  for (int r = 1; r <= order; r++) {
    // from bottom to top to save a buff
    for (int i = order - 1; i >= r; i--) {
      double a_left = localKnots.at(i - 1);
      double a_right = localKnots.at(i + order - r - 1);
      double alpha;
      if (a_right == a_left)
        alpha = 0.; // not sure??
      else
        alpha = (x - a_left) / (a_right - a_left);
      pts.at(i) = (1. - alpha) * pts.at(i - 1) + alpha * pts.at(i);
    }
  }
  return pts.at(order - 1);
}
double BSpline::evalFirstDeriv(double x) const {
  // first find the interval of x in knots
  int leftKnot = order - 1;
  int leftPt = 0;
  while (knots.at(leftKnot + 1) < x) {
    leftKnot++;
    leftPt++;
  }
  int order_t = order - 1;
  vector<double> pts(&(ctrlPts_1st.at(leftPt)),
                     &(ctrlPts_1st[leftPt + order_t]));
  vector<double> localKnots(&(knots.at(leftKnot - order_t + 2)),
                            &(knots[leftKnot + order_t]));
  for (int r = 1; r <= order_t; r++) {
    // from bottom to top to save a buff
    for (int i = order_t - 1; i >= r; i--) {
      double a_left = localKnots.at(i - 1);
      double a_right = localKnots.at(i + order_t - r - 1);
      double alpha;
      if (a_right == a_left)
        alpha = 0.; // not sure??
      else
        alpha = (x - a_left) / (a_right - a_left);
      pts.at(i) = (1. - alpha) * pts.at(i - 1) + alpha * pts.at(i);
    }
  }
  return pts.at(order_t - 1);
}
double BSpline::evalSecondDeriv(double x) const {
  if (order == 2)
    return 0;
  // first find the interval of x in knots
  int leftKnot = order - 1;
  int leftPt = 0;
  while (knots.at(leftKnot + 1) < x) {
    leftKnot++;
    leftPt++;
  }
  int order_t = order - 2;
  vector<double> pts(&(ctrlPts_2nd.at(leftPt)),
                     &(ctrlPts_2nd[leftPt + order_t]));
  vector<double> localKnots(&(knots.at(leftKnot - order_t + 2)),
                            &(knots[leftKnot + order_t]));
  for (int r = 1; r <= order_t; r++) {
    // from bottom to top to save a buff
    for (int i = order_t - 1; i >= r; i--) {
      double a_left = localKnots.at(i - 1);
      double a_right = localKnots.at(i + order_t - r - 1);
      double alpha;
      if (a_right == a_left)
        alpha = 0.; // not sure??
      else
        alpha = (x - a_left) / (a_right - a_left);
      pts.at(i) = (1. - alpha) * pts.at(i - 1) + alpha * pts.at(i);
    }
  }
  return pts.at(order_t - 1);
}
void BSpline::print() {
  cout << " ctrlPts " << ctrlPts.size() << endl;
  for (int i = 0; i < ctrlPts.size(); i++)
    cout << ctrlPts.at(i) << " ";
  cout << endl;
  cout << " knots " << knots.size() << endl;
  for (int i = 0; i < knots.size(); i++)
    cout << knots.at(i) << " ";
  cout << endl;
}
void BSpline::getpara(int &order_p, vector<double> &ctrlPts_p,
                      vector<double> &knots_p, vector<double> &weight_p) {
  order_p = order;
  ctrlPts_p = ctrlPts;
  knots_p = knots;
  weight_p = weight;
}
//draft provided by Claude Sonnet 4 from prompt:
// Given a 2d point along a bspline, specified in cartesian coordinates,
// how do I determine its parametric coordinate on the spline?
/**
 * Find parametric location of a given point using Newton-Raphson method
 * @param targetPt (in) The cartesian point to find parameter for
 * @param initialGuess (in) initial guess of parametric location
 * @param tolerance (in) Convergence tolerance
 * @param maxIterations (in) Maximum number of iterations
 * @param tMin (in) Minimum valid parameter value
 * @param tMax (in) Maximum valid parameter value
 * @return parametric location, or NaN if failed to converge
 */
double BSpline::newtonRaphson(const double targetPt,
                              const double initialGuess,
                              const double tolerance,
                              const int maxIterations,
                              const double tMin,
                              const double tMax) const {
  double t = initialGuess;

  for (int i = 0; i < maxIterations; ++i) {
    // Evaluate spline and derivative at current t
    double currentPt = eval(t);
    double derivative = evalFirstDeriv(t);

    // Calculate error
    double error = currentPt - targetPt;

    // Check convergence
    if (std::abs(error) < tolerance) {
      return t;
    }

    // Check if derivative is too small (avoid division by zero)
    double derivNormSq = derivative*derivative;
    if (std::abs(derivNormSq) < 1e-12) {
      break; // Derivative too small, cannot continue
    }

    // Newton-Raphson update: dt = -(error · derivative) / (derivative · derivative)
    double dt = -(error*derivative) / derivNormSq;
    t += dt;

    // Clamp t to valid range
    t = std::max(tMin, std::min(tMax, t));
  }

  return std::numeric_limits<double>::quiet_NaN(); // Failed to converge
}

double BSpline::invEval(double targetPt, double guess, bool debug) const {
  if(debug) {
    std::cerr << "targetPt " << targetPt << " guess " << guess << "\n";
  }

  // Find best initial guess by sampling
  const int numSamples = knots.size()*2;
  const double tolerance = 1e-10;
  const int maxIterations = 50;
  const double tMin = 0.0;
  const double tMax = 1.0;
  double bestT = tMin;
  double minDistance = std::numeric_limits<double>::max();

  std::vector<double> guesses;
  guesses.reserve(numSamples);
  for (int i = 0; i <= numSamples; ++i) {
    double t = tMin + (tMax - tMin) * i / numSamples;
    double point = eval(t);
    double distance = std::abs(point - targetPt);

    if (distance < minDistance) {
      minDistance = distance;
      bestT = t;
      guesses.push_back(t);
    }
  }
  if(guess != -1) {
    guesses.push_back(guess);
  }

  for(auto rit = guesses.rbegin(); rit != guesses.rend(); ++rit) {
    const auto guess = *rit;
    if(debug) {
      std::cerr << " guess " << guess << "\n";
    }
    const auto res = newtonRaphson(targetPt, guess, tolerance, maxIterations, tMin, tMax);
    if(! std::isnan(res) ) {
      return res;
    }
  }
  std::cerr << "ERROR: " << __func__ << " could not find parametric coordinate for input point "
            << targetPt << '\n';
  assert(false);
  exit(EXIT_FAILURE);
}

// H.Prautzsch Springer,2002
BSpline &BSpline::operator=(const PolyNomial &pn) {
  vector<double> coffs_p;
  pn.getcoeffs(coffs_p);
  order = coffs_p.size();
  ctrlPts.clear();
  knots.clear();
  weight.clear();
  ctrlPts_1st.clear();
  ctrlPts_2nd.clear();
  knots.resize(order + order);
  ctrlPts.resize(order);
  for (int i = 0; i < order; i++) {
    knots.at(i) = 0.;
    knots.at(order + i) = 1.0;
    double ctr = 0.0;
    for (int j = 0; j <= i; j++) {
      ctr += (double)calcuBinomial(order - 1 - j, i - j) *
             coffs_p.at(order - 1 - j);
    }
    ctr /= (double)calcuBinomial(order - 1, i);
    ctrlPts.at(i) = ctr;
  }
  double tol = 1e-6;
  assert(fabs(ctrlPts.at(0) - pn.eval(0.0)) < tol);
  assert(fabs(ctrlPts.at(order - 1) - pn.eval(1.0)) < tol);
  calcuDerivCoeff();
  return *this;
}
