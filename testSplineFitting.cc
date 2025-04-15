#include "SimUtil.h"
#include "SimModel.h"
#include "SimInfo.h"
#include "SimCreateModel.h"
#include "SimDisplay.h"
#include "MeshSim.h"
#include <iostream>
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <tuple>
#include <array>
#include <vector>
#include <cassert>
#include <algorithm> //[min|max]element
#include <map>
#include <limits> //std::numeric_limits

#include <unsupported/Eigen/Splines>

using namespace std;

void messageHandler(int type, const char *msg);

namespace LandiceMeshGen {
class Spline {
  public:
    Spline(std::vector<double> pts_in, std::vector<double> knots_in, int dim_in, int order_in, int numPts_in) :
      pts(pts_in), knots(knots_in), dim(dim_in), order(order_in), numPts(numPts_in) {}
    std::vector<double>& getKnots() { return knots; }
    std::vector<double>& getPts() { return pts; }
    size_t getOrder() { return order; }
    size_t getDim() { return dim; }
    size_t getNumPts() { return numPts; }
  private:
    std::vector<double> pts;
    std::vector<double> knots;
    size_t order;
    size_t dim;
    size_t numPts;
};
}

pCurve createPCurveFromSpline(LandiceMeshGen::Spline& spline) {
  bool debug = true;
  //transpose the points and knots vector for simmetrix
  const auto& sPts = spline.getPts();
  const auto dim = spline.getDim();
  const auto numPts = spline.getNumPts();
  std::vector<double> pts(dim*numPts);
  for( int row = 0; row < dim; row++) {
    for( int col = 0; col < numPts; col++) {
      pts[col*dim + row] = sPts.at(row*numPts + col);
    }
  }
  if (debug) {
    std::cout << "sPts: "; for (auto val : sPts) std::cout << val << " "; std::cout << "\n";
    std::cout << "pts: "; for (auto val : pts) std::cout << val << " "; std::cout << "\n";
  }
  return SCurve_createBSpline(spline.getOrder(), spline.getNumPts(),
      spline.getPts().data(),spline.getKnots().data(),NULL);

}

//assumption here is that all the x coordinates are listed first, then y, and
//then z
LandiceMeshGen::Spline fit3dCubicSplineToPoints(std::vector<double> pts)
{
  bool debug = true;
  const auto dim = 3;
  const auto degree = 3;
  assert(pts.size() % dim == 0);
  const auto numPts = pts.size() / dim;
  Eigen::Array<double,dim,Eigen::Dynamic> points(dim,numPts);
  int ptsIdx = 0;
  for( int row = 0; row < dim; row++) {
    for( int col = 0; col < numPts; col++) {
      points(row,col) = pts.at(ptsIdx++);
    }
  }
  const Eigen::Spline3d spline = Eigen::SplineFitting<Eigen::Spline3d>::Interpolate(points, degree);
  const Eigen::Spline3d::KnotVectorType eKnots = spline.knots(); //Array<scalar, 1, dynamic>
  std::vector<double> knots {eKnots.begin(), eKnots.end()};
  if(debug)
    std::cout << "knots: "; for (auto val : knots) std::cout << val << " "; std::cout << "\n";
  return LandiceMeshGen::Spline(pts, knots, dim, degree, numPts);
}

int main(int argc, char **argv)
{
  pGImporter importer;  // the importer object used to create the geometry
  pGVertex* vertices;    // array to store the returned model vertices
  pGEdge* edges;     // array to store the returned model edges
  pGFace* faces;      // array to store the returned model faces
  pGModel model;        // pointer to the complete model

  const auto prefix = std::string("foo");
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
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    // Create the pGImporter object
    importer = GImporter_new();

    const auto numPts = 4;
    std::vector<double> pts = {0, 4, 8, 12,  //x
                               0, 1, 1 ,0,   //y
                               0, 0, 0 ,0};  //z

    auto spline = fit3dCubicSplineToPoints(pts);

    auto curve = createPCurveFromSpline(spline);

    vertices = new pGVertex[2];
    double startVtx[3] = {0,0,0};
    vertices[0] = GImporter_createVertex(importer, startVtx);
    double endVtx[3] = {12,0,0};
    vertices[1] = GImporter_createVertex(importer, endVtx);

    edges = new pGEdge[1];
    edges[0] = GImporter_createEdge(importer, vertices[0], vertices[1],
        curve, 0, 1, 1);

    // Now complete the model and delete the importer
    model = GImporter_complete(importer);
    auto isValid = GM_isValid(model,2,NULL);
    if(!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    } else {
      cout << "Model is valid.\n";
    }
    GImporter_delete(importer);
    
    cout<<"Number of vertices in model: "<<GM_numVertices(model)<<endl;
    cout<<"Number of edges in model: "<<GM_numEdges(model)<<endl;
    cout<<"Number of faces in model: "<<GM_numFaces(model)<<endl;
    cout<<"Number of regions in model: "<<GM_numRegions(model)<<endl;
    GM_write(model,modelFileName.c_str(),0,0);

    // cleanup
    GM_release(model);
    Progress_delete(progress);
    MS_exit();
    Sim_unregisterAllKeys();
    SimModel_stop();
    Sim_logOff();

  } catch (pSimInfo err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimInfo_code(err)<<endl;
    cerr<<"  Error string: "<<SimInfo_toString(err)<<endl;
    SimInfo_delete(err);
    return 1;
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
    return 1;
  }
  return 0;
}

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    cout<<"Info: "<<msg<<endl;
    break;
  case Sim_DebugMsg:
    cout<<"Debug: "<<msg<<endl;
    break;
  case Sim_WarningMsg:
    cout<<"Warning: "<<msg<<endl;
    break;
  case Sim_ErrorMsg:
    cout<<"Error: "<<msg<<endl;
    break;
  }
  return;
}
