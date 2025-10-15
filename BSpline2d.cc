#include <cmath>
#include <iostream>
#include "splineInterpolation.h"

namespace SplineInterp {

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
double BSpline2d::newtonRaphson(const Point2d& targetPt,
                                const double initialGuess,
                                const double tolerance,
                                const int maxIterations,
                                const double tMin,
                                const double tMax) const {
  double t = initialGuess;

  for (int i = 0; i < maxIterations; ++i) {
    // Evaluate spline and derivative at current t
    Point2d currentPoint = {x.eval(t), y.eval(t)};
    Point2d derivative = {x.evalFirstDeriv(t), y.evalFirstDeriv(t)};

    // Calculate error vector
    Point2d error = currentPoint - targetPt;

    // Check convergence
    if (error.norm() < tolerance) {
      return t;
    }

    // Check if derivative is too small (avoid division by zero)
    double derivNormSq = derivative.dot(derivative);
    if (derivNormSq < tolerance) {
      break; // Derivative too small, cannot continue
    }

    // Newton-Raphson update: dt = -(error · derivative) / (derivative · derivative)
    double dt = -error.dot(derivative) / derivNormSq;
    t += dt;

    // Clamp t to valid range
    t = std::max(tMin, std::min(tMax, t));
  }

  return std::numeric_limits<double>::quiet_NaN(); // Failed to converge
}

double BSpline2d::invEval(const Point2d& targetPt, double guess, bool debug) const {
  if(debug) {
    std::cerr << "targetPt " << targetPt.x << ", " << targetPt.y << " guess " << guess << "\n";
  }

  // Find best initial guess by sampling
  const int numSamples = x.getNumKnots()*2;
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
    Point2d point = {x.eval(t), y.eval(t)};
    double distance = (point - targetPt).norm();

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
            << targetPt.x << ", " << targetPt.y << '\n';
  return std::numeric_limits<double>::quiet_NaN(); // Failed to converge
}

} //end namespace
