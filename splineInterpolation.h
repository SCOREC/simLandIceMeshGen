#include "BSpline.h"
#include <vector>

namespace SplineInterp {

class BSpline2d {
public:
  M3DC1::BSpline x;
  M3DC1::BSpline y;
};

/**
 * Interpolate a spline through the given points using a clamped
 * b-spline as underlying representation.
 */
BSpline2d fitCubicSplineToPoints(std::vector<double> xpts,
                                 std::vector<double> ypts);

/**
 * Check the orientation of a curve. The method is applicable
 * to non-convex polygons too.
 * Returns true if the curve is clockwise, false if counter-clockwise.
 */
bool curveOrientation(std::vector<double> &curvePts);

} // namespace SplineInterp
