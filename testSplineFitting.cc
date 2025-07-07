#include "MeshSim.h"
#include "SimAdvModel.h"
#include "SimDisplay.h"
#include "SimInfo.h"
#include "SimModel.h"
#include "SimUtil.h"
#include "splineInterpolation.h"
#include <cassert>
#include <iostream>
#include <math.h>
#include <string>

using namespace std;

void messageHandler(int type, const char *msg);

pGEdge fitCurveToContour(pGRegion region, pGVertex first, pGVertex last, const int numPts,
                         std::vector<double>& pts, bool debug=false) {
  assert(numPts > 1);
  int order;
  if( numPts == 2 || numPts == 3) {
    order = numPts;
  } else {
    order = 4;
  }
  pCurve curve =
    SCurve_createInterpolatedBSpline(order, numPts, &pts[0], NULL);
  pGEdge edge = GR_createEdge(region, first, last, curve, 1);
  return edge;
}

int main(int argc, char **argv) {
  pGVertex *vertices; // array to store the returned model vertices
  pGEdge *edges;      // array to store the returned model edges
  pGFace *faces;      // array to store the returned model faces
  pGRegion outerRegion;
  pGIPart part;
  pGModel model; // pointer to the complete model

  const auto prefix = std::string("foo");
  std::string modelFileName = prefix + ".smd";
  std::string meshFileName = prefix + ".sms";

  try {
    Sim_logOn("simMeshGen.log");
    SimModel_start(); // Call before Sim_readLicenseFile
    Sim_readLicenseFile(0);
    MS_init();

    Sim_setMessageHandler(messageHandler);
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    model = GM_new(1);
    part = GM_rootPart(model);
    outerRegion = GIP_outerRegion(part);

    const auto numPts = 4;
    std::vector<double> pts = {0,  0, 0,
                               4,  1, 0,
                               8,  1, 0,
                               12, 0, 0};

    double startPt[3] = {pts[0], pts[1], pts[2]};
    auto startVtx = GR_createVertex(outerRegion, startPt);

    double endPt[3] = {pts[(numPts-1)*3],
                       pts[(numPts-1)*3+1],
                       pts[(numPts-1)*3+2]};
    auto endVtx = GR_createVertex(outerRegion, endPt);

    auto edge = fitCurveToContour(outerRegion, startVtx, endVtx, numPts, pts, true);

    auto isValid = GM_isValid(model, 2, NULL);
    if (!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    } else {
      cout << "Model is valid.\n";
    }
    assert(GM_numVertices(model) == 2);
    assert(GM_numEdges(model) == 1);
    assert(GM_numFaces(model) == 0);
    assert(GM_numRegions(model) == 0);
    //assert(std::fabs(GE_length(edge) - 12.2685) <= 1e-5);
    GM_write(model, modelFileName.c_str(), 0, 0);

    // cleanup
    GM_release(model);
    Progress_delete(progress);
    MS_exit();
    Sim_unregisterAllKeys();
    SimModel_stop();
    Sim_logOff();

  } catch (pSimInfo err) {
    cerr << "SimModSuite error caught:" << endl;
    cerr << "  Error code: " << SimInfo_code(err) << endl;
    cerr << "  Error string: " << SimInfo_toString(err) << endl;
    SimInfo_delete(err);
    return 1;
  } catch (...) {
    cerr << "Unhandled exception caught" << endl;
    return 1;
  }
  return 0;
}

void messageHandler(int type, const char *msg) {
  switch (type) {
  case Sim_InfoMsg:
    cout << "Info: " << msg << endl;
    break;
  case Sim_DebugMsg:
    cout << "Debug: " << msg << endl;
    break;
  case Sim_WarningMsg:
    cout << "Warning: " << msg << endl;
    break;
  case Sim_ErrorMsg:
    cout << "Error: " << msg << endl;
    break;
  }
  return;
}
