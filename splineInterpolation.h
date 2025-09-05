#include "BSpline.h"
#include <vector>
#include <string>

namespace SplineInterp {

class BSpline2d {
public:
  Spline::BSpline x;
  Spline::BSpline y;
};

struct SplineInfo {
  std::vector<SplineInterp::BSpline2d> splines;
  SplineInfo(int numSplines) {
    splines.reserve(numSplines);
  }
  void addSpline(SplineInterp::BSpline2d spline) {
    splines.emplace_back(spline);
  }
  void writeToOsh(std::string filename);
  void writeSamplesToCsv(std::string fileName);
};

/**
 * Interpolate a spline through the given points using a clamped
 * b-spline as underlying representation.
 */
BSpline2d fitCubicSplineToPoints(std::vector<double> pts);
BSpline2d fitCubicSplineToPoints(std::vector<double> xpts,
                                 std::vector<double> ypts);
/**
 * create a PWL curve through the given points
 */
BSpline2d attach_piecewise_linear_curve(std::vector<double> points);

/**
 * Check the orientation of a curve. The method is applicable
 * to non-convex polygons too.
 * Returns true if the curve is clockwise, false if counter-clockwise.
 */
bool curveOrientation(std::vector<double> &curvePts);

} // namespace SplineInterp
