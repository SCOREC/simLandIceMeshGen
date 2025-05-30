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
    std::vector<double> xpts = {0, 4, 8, 12};
    std::vector<double> ypts = {0, 1, 1, 0};

    auto bspline = SplineInterp::fitCubicSplineToPoints(xpts, ypts);

    vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
    int order;
    bspline.x.getpara(order, ctrlPtsX, knots, weight);
    bspline.y.getpara(order, ctrlPtsY, knots, weight);
    const int numCtrlPts = ctrlPtsX.size();
    vector<double> ctrlPts3D(3 * (numCtrlPts));
    for (int k = 0; k < numCtrlPts; k++) {
      ctrlPts3D.at(3 * k) = ctrlPtsX.at(k);
      ctrlPts3D.at(3 * k + 1) = ctrlPtsY.at(k);
      ctrlPts3D[3 * k + 2] = 0.0;
    }
    // To make it consistent, we will define every edge in counter-clockwise
    // direction. If curve is clockwise, set edge dir to 0, otherwise 1 to
    // follow the above convention.
    int edgeDir = 1;
    bool clockwise = SplineInterp::curveOrientation(ctrlPts3D);
    if (clockwise)
      edgeDir = 0;

    pCurve curve =
        SCurve_createBSpline(order, numCtrlPts, &ctrlPts3D[0], &knots[0], NULL);

    vertices = new pGVertex[2];
    double startVtx[3] = {0, 0, 0};
    vertices[0] = GR_createVertex(outerRegion, startVtx);
    double endVtx[3] = {12, 0, 0};
    vertices[1] = GR_createVertex(outerRegion, endVtx);

    pGEdge edge =
        GR_createEdge(outerRegion, vertices[0], vertices[1], curve, edgeDir);

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
    assert(std::fabs(GE_length(edge) - 12.2685) <= 1e-5);
    GM_write(model, modelFileName.c_str(), 0, 0);

    // cleanup
    GM_release(model);
    delete[] vertices;
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
