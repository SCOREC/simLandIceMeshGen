#include "MeshSim.h"
#include "SimAdvModel.h"
#include "SimDisplay.h"
#include "SimInfo.h"
#include "SimModel.h"
#include "SimUtil.h"
#include "splineInterpolation.h"
#include "curveReader.h"
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
  const int numExpectedArgs = 3;
  if (argc != numExpectedArgs) {
    std::cerr << "Usage: <input csv file> <expected curve length>\n";
    std::cerr << "input csv file with the following columns: "
                 "x,y,z,isOnCurve,angle,isMdlVtx\n";
    return 1;
  }
  assert(argc == numExpectedArgs);

  std::string curveFilename = argv[1];
  double expectedCurveLength = std::stod(argv[2]);

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

    auto curve = CurveReader::readCurveInfo(curveFilename);
    const auto numPts = curve.x.size();

    double startPt[3] = {curve.x[0], curve.y[0], 0.0};
    auto startVtx = GR_createVertex(outerRegion, startPt);

    double endPt[3] = {curve.x[(numPts-1)], curve.y[(numPts-1)], 0.0};
    auto endVtx = GR_createVertex(outerRegion, endPt);

    vector<double> pts;
    for(int i = 0; i < numPts; ++i) {
      pts.push_back(curve.x[i]);
      pts.push_back(curve.y[i]);
      pts.push_back(0);
    }

    auto edge = fitCurveToContour(outerRegion, startVtx, endVtx, numPts, pts, true);

    //Fit curve using Spline2D Implementation
    auto bspline = SplineInterp::fitCubicSplineToPoints(curve.x, curve.y);
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

    std::cout << bspline.x.eval(0) << ", " << bspline.y.eval(0) << std::endl;
    std::cout << bspline.x.eval(1) << ", " << bspline.y.eval(1) << std::endl;

    pCurve spline2DCurve =
        SCurve_createBSpline(order, numCtrlPts, &ctrlPts3D[0], &knots[0], NULL);
    pGEdge spline2DEdge =
        GR_createEdge(outerRegion, startVtx, endVtx, spline2DCurve, edgeDir);

    auto isValid = GM_isValid(model, 2, NULL);
    if (!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    } else {
      cout << "Model is valid.\n";
    }
    assert(GM_numVertices(model) == 2);
    assert(GM_numEdges(model) == 2);
    assert(GM_numFaces(model) == 0);
    assert(GM_numRegions(model) == 0);
    std::cout << GE_length(edge) << " " << GE_length(spline2DEdge) << std::endl;
    assert(std::fabs(GE_length(edge) - GE_length(spline2DEdge)) <= 1e-5);
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
