#ifndef NETCDF_WRITER_H
#define NETCDF_WRITER_H

#include "MeshSim.h"
#include <string>
int writeMeshSimToNetCDF(pMesh mesh, pGModel model, std::string outputFileName, bool convertKmToMeters);

#endif
