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

using namespace std;

void messageHandler(int type, const char *msg);

/* DGESV prototype */
extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                 double* b, int* ldb, int* info );

// **********************************************
void interpolateCubicBSpline( vector<double>& points,vector<double>& knots, vector<double> &ctrlPoints, int bc)
// **********************************************
{
  int numPts=points.size();
  int order_p=4;
  assert(numPts>1);
  ctrlPoints.resize(numPts+2);
  vector< double> coeffs((numPts+2)*(numPts+2),0.); // numPts+2 * numPts+2 linear system
  // first find the constraint to coninside with points
  // 2 natural cubic spline at the boundary
  for( int i=0; i<numPts+2; i++)
  {
    vector<double> points_tmp(numPts+2,0.);
    points_tmp.at(i)=1.0;
    vector<double>  weight_p;
    M3DC1::BSpline  basis(order_p, points_tmp, knots, weight_p);
    for( int j=0; j<numPts; j++)
    {
      double para=knots.at(order_p+j-1);
      double res= basis.eval(para);
      double secondDeriv0=basis.evalSecondDeriv(0);
      double secondDeriv1=basis.evalSecondDeriv(1);
      coeffs.at(i*(numPts+2)+j)=res;
    } 
    double secondDeriv0=basis.evalSecondDeriv(0);
    double secondDeriv1=basis.evalSecondDeriv(1);
    if(bc==0) // natural
    {
      coeffs.at(i*(numPts+2)+numPts)=secondDeriv0;
      coeffs.at(i*(numPts+2)+numPts+1)=secondDeriv1;
    }
    else // periodic
    {
      double firstDeriv0=basis.evalFirstDeriv(0);
      double firstDeriv1=basis.evalFirstDeriv(1);
      coeffs.at(i*(numPts+2)+numPts)=firstDeriv0-firstDeriv1;
      coeffs.at(i*(numPts+2)+numPts+1)=secondDeriv0-secondDeriv1; 
    }
  }

  // set up the linear system and solve
  vector<double> rhs(numPts+2,0.0);
  for( int i=0; i<numPts; i++)
    rhs.at(i)=points.at(i);

  int info,one=1, dim=numPts+2;
  vector<int> ipiv(dim,0);
  dgesv_( &dim, &one,& (coeffs.at(0)), &dim, &(ipiv.at(0)), &(rhs.at(0)), &dim, &info );
  assert( info==0);
  for ( int i=0; i<numPts+2; i++)
    ctrlPoints.at(i)=rhs.at(i);  
}

// From Usman's m3dc1_model.cc code
// use clamped b-spline as underlying representation
// **********************************************
void fitCubicSplineToPoints(int * numPts, double * points)
// **********************************************
{
  int order_p=4; 
  int knotsize=2*order_p+*numPts-2;
  vector<double> knots(knotsize,0.);
  vector<double> ctrlPointsX(*numPts+2),ctrlPointsY(*numPts+2),weight;
  for( int i=0; i<order_p; i++)
  {
    knots.at(knotsize-i-1)=1.0;
  }
  double increment=1.0/(*numPts-1);
  for (int i=0; i<*numPts-2; i++)
  {
    //double increment=inter_len.at(i)/len;
    knots.at(order_p+i)=knots.at(order_p+i-1)+increment;
  }
  vector<double> pointsX(*numPts),pointsY(*numPts);
  for( int i=0; i<*numPts; i++)
  {
    pointsX.at(i)=points[2*i];
    pointsY.at(i)=points[2*i+1];
  }
  interpolateCubicBSpline(pointsX,knots,ctrlPointsX,0);
  interpolateCubicBSpline(pointsY,knots,ctrlPointsY,0);
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
    std::vector<double> pts = {0.,0.,0., 4.,1.,0., 8.,1.,0., 12.,0.,0.};

    auto curve = SCurve_createPiecewiseLinear(numPts, pts.data());

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
