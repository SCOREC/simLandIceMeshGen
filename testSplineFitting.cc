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
#include <fstream>

using namespace std;

void messageHandler(int type, const char *msg);

void writeDefinition(std::string fileNameNoExt, SplineInterp::BSpline2d& bspline) {
    vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
    int order;
    bspline.x.getpara(order, ctrlPtsX, knots, weight);
    bspline.y.getpara(order, ctrlPtsY, knots, weight);

    std::cout << fileNameNoExt << std::endl;
    std::ofstream splineInterpDefinitionFile(fileNameNoExt + "_Spline_Interp_Definition.csv");
      if (!splineInterpDefinitionFile) {
      std::cerr << "Failed to open output file for splineinterp definition.\n";
      exit(EXIT_FAILURE);
    }
    splineInterpDefinitionFile << "# Control Points" << std::endl << "X,Y,Z,Weight" << std::endl;

    const int numCtrlPts = ctrlPtsX.size();
    for(int i = 0; i < numCtrlPts; ++i) {
      splineInterpDefinitionFile << ctrlPtsX[i] << "," << ctrlPtsY[i] << ",0,";
      if(weight.size() > i)
        splineInterpDefinitionFile << weight[i] << std::endl;
      else
        splineInterpDefinitionFile << "1" << std::endl;
    }
    splineInterpDefinitionFile << std::endl << "# Knot Vector" << std::endl;
    for(int i = 0; i < knots.size(); ++i)
      splineInterpDefinitionFile << knots[i] << std::endl;
    splineInterpDefinitionFile.close();
}

void writeSamples(std::string fileNameNoExt, const CurveReader::CurveInfo& curve, SplineInterp::BSpline2d& bspline) {
    vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
    int order;
    bspline.x.getpara(order, ctrlPtsX, knots, weight);
    bspline.y.getpara(order, ctrlPtsY, knots, weight);

    std::ofstream splineInterpSampleFile(fileNameNoExt + "_Spline_Fitting_Samples.csv");
    if (!splineInterpSampleFile) {
      std::cerr << "Failed to open output file for splineinterp sampled points.\n";
      exit(EXIT_FAILURE);
    }
    auto numSamples = curve.x.size() * 25;
    splineInterpSampleFile << "x, y, isVertex\n";
    splineInterpSampleFile << curve.x[0] << ", " << curve.y[0] << ", 1\n";
    for(int i = 0; i < numSamples; ++i) {
      auto t = 1.0 * i / numSamples;
      splineInterpSampleFile << bspline.x.eval(t) << ", " << bspline.y.eval(t) << ", 0\n";
    }
    splineInterpSampleFile << curve.x.back() << ", " << curve.y.back() << ", 1\n";
    splineInterpSampleFile.close();
}

double createSimModelEdge(const CurveReader::CurveInfo& curve, SplineInterp::BSpline2d& bspline) {
  pGRegion outerRegion;
  pGIPart part;
  pGModel model; // pointer to the complete model

  const auto prefix = std::string("foo");
  std::string modelFileName = prefix + ".smd";
  std::string meshFileName = prefix + ".sms";

  double length = 0;

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

    double startPt[3] = {curve.x[0], curve.y[0], 0.0};
    auto startVtx = GR_createVertex(outerRegion, startPt);

    double endPt[3] = {curve.x.back(), curve.y.back(), 0.0};
    auto endVtx = GR_createVertex(outerRegion, endPt);

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
    assert(GM_numEdges(model) == 1);
    assert(GM_numFaces(model) == 0);
    assert(GM_numRegions(model) == 0);

    length = GE_length(spline2DEdge);
    std::cout << "Simmetrix Edge Length: " << length << "\n";
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
  return length;
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
  int extensionPos = curveFilename.rfind(".");
  int slashPos = curveFilename.rfind("/");
  std::string fileNameNoExt = curveFilename.substr(slashPos + 1, extensionPos);
  double expectedCurveLength = std::stod(argv[2]);

    auto curve = CurveReader::readCurveInfo(curveFilename);

    //Fit curve using Spline2D Implementation
    auto bspline = SplineInterp::fitCubicSplineToPoints(curve.x, curve.y);
    writeDefinition(fileNameNoExt, bspline);
    writeSamples(fileNameNoExt, curve, bspline);

    //Fit curve using the simmetrix APIs
    auto length = createSimModelEdge(curve, bspline);

    //compare length
    std::cout << " Expected Length: " << expectedCurveLength << std::endl;
    assert(std::fabs(length - expectedCurveLength) <= 1e-5);
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
