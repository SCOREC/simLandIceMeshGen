#include "simModelGen2d.h"
#include "curveReader.h"
#include <numeric> //std::accumulate

void messageHandler(int type, const char *msg);

int main(int argc, char **argv) {
  const int numExpectedArgs = 3;
  if (argc != numExpectedArgs) {
    std::cerr << "Usage: <input csv file> <output prefix>\n";
    std::cerr << "input csv file with the following columns: "
                 "x,y,z,isOnCurve,angle,isMdlVtx\n";
    return 1;
  }
  assert(argc == numExpectedArgs);

  std::string filename = argv[1];
  const auto prefix = std::string(argv[2]);
  std::cout << "input csv file: " << argv[1] << " "
            << "output prefix: " << prefix << "\n";

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

    auto curveInfo = CurveReader::readCurveInfo(filename);
    CurveReader::printCurveInfo(curveInfo);
    GeomInfo geom;
    geom.numVtx = curveInfo.isMdlVtx.size();
    geom.vtx_x = curveInfo.x;
    geom.vtx_y = curveInfo.y;
    geom.firstContourPt = 0;
    const auto numMdlVerts = std::accumulate(curveInfo.isMdlVtx.begin(), curveInfo.isMdlVtx.end(), 0);
    auto splines = SplineInterp::SplineInfo(numMdlVerts);
    createEdges(mdlTopo, geom, splines, curveInfo.isOnCurve, curveInfo.isMdlVtx,
                debug);

    auto isValid = GM_isValid(mdlTopo.model, 2, NULL);
    if (!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    } else {
      std::cout << "Model is valid.\n";
    }

    printModelInfo(mdlTopo.model);

    GM_write(mdlTopo.model, modelFileName.c_str(), 0, 0);

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
