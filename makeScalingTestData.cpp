#include "makeScalingTestData.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalConstants.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecutionSpace::memory_space;

//This function creates points that are on a circle 
//with specified center and radius
Kokkos::View<double*[2], MemSpace> makeCircle(double centerX, double centerY, double radius, int numPoints) {
  Kokkos::View<double*[2], MemSpace> res("result", numPoints);
  double angle = (2*Kokkos::numbers::pi)/numPoints;
  Kokkos::parallel_for("parallel calculation of points on a circle", numPoints, KOKKOS_LAMBDA(const int i) {
    res(i, 0) = centerX+radius*Kokkos::cos(i*angle);
    res(i, 1) = centerY+radius*Kokkos::sin(i*angle);
  });
  return res;
}

Kokkos::View<double*[2], MemSpace> makeEllipse(double centerX, double centerY, double xRadius, double yRadius, int numPoints) {
  Kokkos::View<double*[2], MemSpace> res("result", numPoints);
  double angle = (2*Kokkos::numbers::pi)/numPoints;
  Kokkos::parallel_for("parallel calculation of points on a circle", numPoints, KOKKOS_LAMBDA(const int i) {
    res(i, 0) = centerX+xRadius*Kokkos::cos(i*angle);
    res(i, 1) = centerY+yRadius*Kokkos::sin(i*angle);
  });
  return res;
}
