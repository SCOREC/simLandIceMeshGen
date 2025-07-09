#include "landIceMeshGen.h"

void messageHandler(int type, const char *msg);

std::string getFileExtension(const std::string &filename) {
  size_t dotPos = filename.rfind('.');
  if (dotPos != std::string::npos) {
    return filename.substr(dotPos);
  }
  return "";
}

int main(int argc, char **argv) {
  const int numExpectedArgs = 7;
  if (argc != numExpectedArgs) {
    std::cerr << "Usage: <jigsaw .msh or .vtk file> <output prefix> "
                 "<coincidentVtxToleranceSquared> <angleTolerance> <createMesh>\n";
    std::cerr << "coincidentVtxToleranceSquared is the mininum allowed "
                 "distance (squared) between adjacent vertices in the "
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
    return 1;
  }
  assert(argc == numExpectedArgs);

  GeomInfo dirty;

  std::string filename = argv[1];
  std::string ext = getFileExtension(filename);
  const auto coincidentPtTol = std::stof(argv[3]);
  const auto prefix = std::string(argv[2]);
  const auto angleTol = std::atof(argv[4]);
  const auto onCurveAngleTol = std::atof(argv[5]);
  const bool doCreateMesh = (std::stoi(argv[6]) == 1);
  std::cout << "input points file: " << argv[1] << " "
            << "coincidentPtTol: " << coincidentPtTol << " "
            << "output prefix: " << prefix << " "
            << "angleTol: " << angleTol << " "
            << "onCurveAngleTol: " << onCurveAngleTol << " "
            << "createMesh: " << doCreateMesh << "\n";

  if (ext == ".vtk") {
    dirty = readVtkGeom(filename);
  } else if (ext == ".msh") {
    dirty = readJigGeom(filename);
  } else {
    std::cerr << "Unsupported file extension: " << ext << "\n";
    return 1;
  }
  convertMetersToKm(dirty);
  auto geom = cleanGeom(dirty, coincidentPtTol, false);
  std::string modelFileName = prefix + ".smd";
  std::string meshFileName = prefix + ".sms";

  const auto debug = true;

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
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    ModelTopo mdlTopo;
    mdlTopo.model = GM_new(1);
    mdlTopo.part = GM_rootPart(mdlTopo.model);
    mdlTopo.region = GIP_outerRegion(mdlTopo.part);

    createBoundingBoxGeom(mdlTopo,geom);

    auto [isPointOnCurve, isMdlVtx] = discoverTopology(geom, angleTol, onCurveAngleTol, debug);

    createEdges(mdlTopo, geom, isPointOnCurve, isMdlVtx, debug);

    createFaces(mdlTopo, geom);

    printModelInfo(mdlTopo.model);

    GM_write(mdlTopo.model, modelFileName.c_str(), 0, 0);

    auto isValid = GM_isValid(mdlTopo.model, 2, NULL);
    if (!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    } else {
      std::cout << "Model is valid.\n";
    }

    if(doCreateMesh) {
      createMesh(mdlTopo, meshFileName, progress);
    }

    // cleanup
    GM_release(mdlTopo.model);
    Progress_delete(progress);
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
