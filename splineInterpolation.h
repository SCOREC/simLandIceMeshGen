#include "BSpline.h"
#include <vector>
#include <string>
#include <cmath>

namespace SplineInterp {

struct Point2d {
    double x, y;

    Point2d(double x = 0, double y = 0) : x(x), y(y) {}

    Point2d operator-(const Point2d& other) const {
        return Point2d(x - other.x, y - other.y);
    }

    Point2d operator+(const Point2d& other) const {
        return Point2d(x + other.x, y + other.y);
    }

    Point2d operator*(double scalar) const {
        return Point2d(x * scalar, y * scalar);
    }

    double dot(const Point2d& other) const {
        return x * other.x + y * other.y;
    }

    double norm() const {
        return std::sqrt(x * x + y * y);
    }
};

class BSpline2d {
public:
  Spline::BSpline x;
  Spline::BSpline y;
  double invEval(const Point2d& targetPt, double guess, bool debug=false) const;
private:
  double newtonRaphson(const Point2d& targetPt,
                       const double initialGuess,
                       const double tolerance = 1e-10,
                       const int maxIterations = 50,
                       const double tMin = 0,
                       const double tMax = 1) const;
};

struct SplineInfo {
  std::vector<SplineInterp::BSpline2d> splines;
  SplineInfo(int numSplines) {
    splines.reserve(numSplines);
  }
  void addSpline(SplineInterp::BSpline2d spline) {
    splines.emplace_back(spline);
  }
  int size() {
    return splines.size();
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
BSpline2d attach_piecewise_linear_curve(std::vector<double> xpts, std::vector<double> ypts);

/**
 * Check the orientation of a curve. The method is applicable
 * to non-convex polygons too.
 * Returns true if the curve is clockwise, false if counter-clockwise.
 */
bool curveOrientation(std::vector<double> &curvePts);

} // namespace SplineInterp
