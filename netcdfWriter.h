#ifndef NETCDF_WRITER_H
#define NETCDF_WRITER_H

#include "MeshSim.h"
#include "modelGen2d.h"
#include <string>
#include <vector>

int writeMeshSimToNetCDF(pMesh mesh, pGModel model, std::string outputFileName,
                          bool convertKmToMeters,
                          const BoundaryPolygons& boundaryPolygons = BoundaryPolygons(),
                          const std::vector<int>& boundaryOrder = {});

#endif
