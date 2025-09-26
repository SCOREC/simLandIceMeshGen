#include "simModelGen2d.h"
#include <numeric> //std::accumulate
#include "Omega_h_file.hpp"

void messageHandler(int type, const char *msg);

std::string getFileExtension(const std::string &filename) {
  size_t dotPos = filename.rfind('.');
  if (dotPos != std::string::npos) {
    return filename.substr(dotPos);
  }
  return "";
}

void writePointParametricCoords(const GeomInfo& geom, const PointClassification& ptClass, const SplineInterp::SplineInfo& sinfo, std::string filename) {
  assert(ptClass.splineIdx.size() == geom.numVtx);
  Omega_h::HostWrite<Omega_h::Real> paraCoords(geom.numVtx*2);
  for(int i=0; i<geom.numVtx; i++) {
    const auto sIdx = ptClass.splineIdx[i];
    const auto bspline = sinfo.splines[sIdx];
    const auto x = geom.vtx_x.at(i); 
    const auto y = geom.vtx_y.at(i); 
    paraCoords[i*2] = bspline.x.invEval(x);
    paraCoords[i*2+1] = bspline.y.invEval(y);
  }
  auto paraCoords_d = Omega_h::read(paraCoords.write());

  std::ofstream file(filename);
  assert(file.is_open());
  const int compressed = 0;
  //the following is from src/Omega_h_file.cpp write(...)
  unsigned char const magic[2] = {0xa1, 0x1a};
  file.write(reinterpret_cast<const char*>(magic), sizeof(magic));
  bool needs_swapping = !Omega_h::is_little_endian_cpu();
  Omega_h::binary::write_value(file, compressed, needs_swapping);
  Omega_h::binary::write_array(file, paraCoords_d, compressed, needs_swapping);
  file.close();
}

int main(int argc, char **argv) {
  const int numExpectedArgs = 8;
  if (argc != numExpectedArgs) {
    std::cerr << "Usage: <jigsaw .msh or .vtk file> <output prefix> "
                 "<coincidentVtxTolerance> <angleTolerance> <createMesh> <units>\n";
    std::cerr << "coincidentVtxTolerance is the mininum allowed "
                 "distance between adjacent vertices in the "
                 "input.  Vertices within the specified distance will "
                 "be merged.\n";
    std::cerr << "angleTolerance defines the upper bound and "
                 "-angleTolerance defines lower bound for determining "
                 "if a vertex bounding two consecutative edges should be "
                 "treated as a model vertex.\n";
    std::cerr << "onCurveAngleTolerance defines the upper bound on the angle "
                 "between adjacent edges in a sequence of four consecutive edges "
                 "used to determine if they are part of the same curve.\n";
    std::cerr << "createMesh = 1:generate mesh, otherwise, "
                 "skip mesh generation.\n";
    std::cerr << "units = m:meters, km:kilometers\n";
    return 1;
  }
  assert(argc == numExpectedArgs);

  ModelFeatures mdl;

  std::string filename = argv[1];
  std::string ext = getFileExtension(filename);
  const auto prefix = std::string(argv[2]);
  const auto coincidentPtTol = std::stof(argv[3]);
  const auto angleTol = std::atof(argv[4]);
  const auto onCurveAngleTol = std::atof(argv[5]);
  const bool doCreateMesh = (std::stoi(argv[6]) == 1);
  const std::string units = argv[7];
  std::cout << "input points file: " << argv[1] << " "
            << "coincidentPtTol: " << coincidentPtTol << " "
            << "output prefix: " << prefix << " "
            << "angleTol: " << angleTol << " "
            << "onCurveAngleTol: " << onCurveAngleTol << " "
            << "createMesh: " << doCreateMesh << " "
            << "units: " << units << "\n";

  assert(units == "m" || units == "km");

  if (ext == ".vtk") {
    mdl = readVtkGeom(filename);
  } else if (ext == ".msh") {
    mdl = readJigGeom(filename);
  } else {
    std::cerr << "Unsupported file extension: " << ext << "\n";
    return 1;
  }
  if(units == "m") {
    convertMetersToKm(mdl);
  }
  const double coincidentPtTolSquared = coincidentPtTol*coincidentPtTol;
  //FIXME For a mesh input, all input points need to be classified on the resulting
  //model and have parametetric coordinates assigned to them.  Otherwise, those
  //points would have to be marked as on the interior which would not be valid.
  //Some of the 'cleanup' can be skipped given the more robust discoverTopology
  //procedure.  Specifically, consecutative points that are very close to each
  //other don't need to be removed.
  cleanGeom(mdl, coincidentPtTolSquared, false); 
  std::string modelFileName = prefix + ".smd";
  std::string meshFileName = prefix + ".sms";

  const auto debug = false;

  // You will want to place a try/catch around all SimModSuite calls,
  // as errors are thrown.
  try {
    Sim_logOn("simMeshGen.log");
    SimModel_start(); // Call before Sim_readLicenseFile
    // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
    // pass in the location of a file containing your keys.  For a release
    // product, use Sim_registerKey()
    Sim_readLicenseFile(0);
    // Tessellation of GeomSim geometry requires Meshing to have started
    MS_init();

    Sim_setMessageHandler(messageHandler);
    pProgress progress = NULL;
    if(debug) {
      pProgress progress = Progress_new();
      Progress_setDefaultCallback(progress);
    }

    ModelTopo mdlTopo;
    mdlTopo.model = GM_new(1);
    mdlTopo.part = GM_rootPart(mdlTopo.model);
    mdlTopo.region = GIP_outerRegion(mdlTopo.part);

    for(auto& geom : mdl.geom) {
      auto [isPointOnCurve, isMdlVtx] = discoverTopology(geom, coincidentPtTolSquared, angleTol, onCurveAngleTol, debug);

      const auto numMdlVerts = isMdlVtx.size() ? std::accumulate(isMdlVtx.begin()+geom.firstContourPt, isMdlVtx.end(), 0) : 0;
      auto splines = SplineInterp::SplineInfo(numMdlVerts+4); //+4 splines for the bounding box
      createBoundingBoxGeom(mdlTopo, geom, splines);

      PointClassification ptClass(geom.numVtx);
      createEdges(mdlTopo, geom, ptClass, splines, isPointOnCurve, isMdlVtx, debug);
    }

    //For now, assume that there is are no near tangencies or contact between
    //features in ModelFeatures.
    //TODO - identify contact and/or near tangencies and fail
    //     - process the geometry to close near tangencies
    //     - process the topology to make it valid - e.g., model edges that
    //       intersect have a model vertex at the point of intersection
    //     - ...
    //     - For the above processing the linkage between the ModelTopo and
    //       ModelFeatures info needs to be such that modification of the 
    //       geometry is easily supported
    //       

    writePointParametricCoords(geom, ptClass, splines, modelFileName + "_parametric.oshb");
    //write the point classification to an omegah binary file
    ptClass.writeToOsh(modelFileName + "_class.oshb");
    //write the bsplines to an omegah binary file
    splines.writeToOsh(modelFileName + "_splines.oshb");
    //write the sampled bsplines to a csv file
    splines.writeSamplesToCsv(modelFileName + "_splines.csv");

    createFaces(mdlTopo, geom);

    printModelInfo(mdlTopo.model);

    GM_write(mdlTopo.model, modelFileName.c_str(), 0, 0);

    auto isValid = GM_isValid(mdlTopo.model, 2, NULL);
    if (!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    }

    if(doCreateMesh) {
      createMesh(mdlTopo, meshFileName, progress);
    }

    // cleanup
    GM_release(mdlTopo.model);
    if(debug) Progress_delete(progress);
    MS_exit();
    Sim_unregisterAllKeys();
    SimModel_stop();
    Sim_logOff();

  } catch (pSimInfo err) {
    std::cerr << "SimModSuite error caught:" << std::endl;
    std::cerr << "  Error code: " << SimInfo_code(err) << std::endl;
    std::cerr << "  Error string: " << SimInfo_toString(err) << std::endl;
    SimInfo_delete(err);
    return 1;
  } catch (...) {
    std::cerr << "Unhandled exception caught" << std::endl;
    return 1;
  }
  return 0;
}

void messageHandler(int type, const char *msg) {
  switch (type) {
  case Sim_InfoMsg:
    std::cout << "Info: " << msg << std::endl;
    break;
  case Sim_DebugMsg:
    std::cout << "Debug: " << msg << std::endl;
    break;
  case Sim_WarningMsg:
    std::cout << "Warning: " << msg << std::endl;
    break;
  case Sim_ErrorMsg:
    std::cout << "Error: " << msg << std::endl;
    break;
  }
  return;
}
