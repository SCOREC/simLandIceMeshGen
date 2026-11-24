#ifndef MAKESCALINGTESTDATA_H
#define MAKESCALINGTESTDATA_H

#include <Kokkos_Core.hpp>
#include <vector>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecutionSpace::memory_space;

Kokkos::View<double *[2], MemSpace> makeCircle(double centerX, double centerY,
                                               double radius, int numPts);
Kokkos::View<double *[2], MemSpace> makeEllipse(double centerX, double centerY,
                                                double xRadius, double yRadius,
                                                int numPts);

#endif
