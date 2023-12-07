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

using namespace std;

void messageHandler(int type, const char *msg);

struct PlaneBounds {
  double minX;
  double maxX;
  double minY;
  double maxY;
};

struct JigGeom {
  int numVtx;
  int numEdges;
  std::vector<double> vtx_x;
  std::vector<double> vtx_y;
  std::vector<std::array<int,2> > edges;
};

std::tuple<string, int> readKeyValue(std::ifstream& in, bool debug=true) {
  std::string key, value;
  std::getline(in, key, '=');
  std::getline(in, value);
  if(debug) std::cout << "line: " << key << " " << value << std::endl;
  const int val = std::stoi(value);
  return {key, val};
}

void skipLine(std::ifstream& in, bool debug=true) {
  std::string line;
  std::getline(in, line);
  if(debug) std::cout << "skip line: " << line << std::endl;
}

PlaneBounds getBoundingPlane(JigGeom& geom) {
  auto minX = std::min_element(geom.vtx_x.begin(), geom.vtx_x.end());
  auto maxX = std::max_element(geom.vtx_x.begin(), geom.vtx_x.end());
  auto minY = std::min_element(geom.vtx_y.begin(), geom.vtx_y.end());
  auto maxY = std::max_element(geom.vtx_y.begin(), geom.vtx_y.end());
  return {*minX, *maxX, *minY, *maxY};
}


std::array<double, 3> readPoint(std::ifstream& in, bool debug=true) {
  std::array<double, 3> pt;
  std::string value;
  std::getline(in, value, ';');
  pt[0] = std::stoi(value);
  std::getline(in, value, ';');
  pt[1] = std::stoi(value);
  std::getline(in, value);
  pt[2] = std::stoi(value);
  return pt;
}

std::array<int, 2> readEdge(std::ifstream& in, bool debug=true) {
  std::array<int, 2> edge;
  std::string value;
  std::getline(in, value, ';');
  edge[0] = std::stoi(value);
  std::getline(in, value, ';');
  edge[1] = std::stoi(value);
  //not using the id of the edge
  std::getline(in, value);
  return edge;
}

JigGeom readJigGeom(std::string fname, bool debug=false) {
  std::ifstream mshFile(fname);
  if ( ! mshFile.is_open() ) {
    fprintf(stderr, "failed to open jigsaw geom file %s\n", fname.c_str());
    exit(EXIT_FAILURE);
  }

  JigGeom geom;

  //header - skip
  skipLine(mshFile,debug);
  //MSHID - skip
  skipLine(mshFile,debug);
  //NDIMS
  {
    auto [key,value] = readKeyValue(mshFile,debug);
    if(debug) std::cout << "key: " << key << " val: " << value << std::endl;
    assert(value == 2);
  }
  //POINT
  {
    auto [key,value] = readKeyValue(mshFile,debug);
    if(debug) std::cout << "key: " << key << " val: " << value << std::endl;
    geom.numVtx = value;
  }
  geom.vtx_x.reserve(geom.numVtx);
  geom.vtx_y.reserve(geom.numVtx);
  //point coordinates
  for(int i=0; i<geom.numVtx; i++) {
    auto pt = readPoint(mshFile,debug);
    geom.vtx_x.push_back(pt[0]);
    geom.vtx_y.push_back(pt[1]);
    if(debug) std::cout << "pt " << geom.vtx_x[i] << ", " << geom.vtx_y[i] << std::endl;
  }
  //EDGE
  {
    auto [key,value] = readKeyValue(mshFile,debug);
    if(debug) std::cout << "key: " << key << " val: " << value << std::endl;
    geom.numEdges = value;
  }
  geom.edges.reserve(geom.numEdges);
  //edge indices
  for(int i=0; i<geom.numEdges; i++) {
    geom.edges.push_back(readEdge(mshFile,debug));
    if(debug) std::cout << "edge " << geom.edges[i][0] << ", " << geom.edges[i][1] << std::endl;
  }

  return geom;
}

int main(int argc, char **argv)
{

  pGImporter importer;  // the importer object used to create the geometry
  pGVertex* vertices;    // array to store the returned model vertices
  pGEdge* edges;     // array to store the returned model edges
  pGFace faces[2];      // array to store the returned model faces
  pGRegion region;      // pointer to returned model region
  pGModel model;        // pointer to the complete model

  auto geom = readJigGeom(argv[1]);
  const auto prefix = std::string(argv[2]);
  std::string modelFileName = prefix + ".smd";
  std::string meshFileName = prefix + ".sms";

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

    vertices = new pGVertex[geom.numVtx];
    edges = new pGEdge[geom.numEdges];

    // First we'll add the vertices
    int i;
    for(i=0; i<4; /*geom.numVtx;*/ i++) {
      double vtx[3] = {geom.vtx_x[i], geom.vtx_y[i], 0};
      vertices[i] = GImporter_createVertex(importer, vtx);
    }

    // Now we'll add the edges
    double point0[3],point1[3];  // xyz locations of the two vertices
    pCurve linearCurve;

    // First, the bottom edges at z=0, connecting the first four vertices in the array
    // 0->1, 1->2, 2->3, 3->0  (indices of the vertices)
    for(i=0; i<4; /*geom.numEdges;*/ i++) {
      const auto startVertIdx = geom.edges[i][0];
      const auto endVertIdx = geom.edges[i][1];
      auto startVert = vertices[startVertIdx];
      auto endVert = vertices[endVertIdx];
      GV_point(startVert, point0);
      GV_point(endVert, point1);
      linearCurve = SCurve_createLine(point0, point1);
      edges[i] = GImporter_createEdge(importer, startVert, endVert, linearCurve, 0, 1, 1);
    }

    auto planeBounds = getBoundingPlane(geom);

    // Now add the faces
    double corner[3], xPt[3], yPt[3];  // the points defining the surface of the face
    pGEdge* faceEdges;                 // the array of edges connected to the face
    int* faceDirs;                     // the direction of the edge with respect to the face

    // When defining the loop, will always start with the first edge in the faceEdges array
    int loopDef[1] = {0}; //FIXME needs to be length two for the outer face, and length one for the inner
    pSurface planarSurface;

    // First the face between the bounding rectangle and the grounding line
    // Define the surface
    corner[0] = planeBounds.minX;
    corner[1] = planeBounds.minY;
    corner[2] = 0;
    xPt[0] = planeBounds.maxX;
    xPt[1] = planeBounds.minY;
    xPt[2] = 0;
    yPt[0] = planeBounds.minX;
    yPt[1] = planeBounds.maxY;
    yPt[2] = 0;
    planarSurface = SSurface_createPlane(corner,xPt,yPt);

    // Create the face
    faceEdges = new pGEdge[4]; //geom.numEdges];
    faceDirs = new int[4]; //geom.numEdges];
    // the first four edges define the outer bounding rectangle
    // - the jigsaw geometry file order is clockwise, reverse the order as
    //   required by GImporter_createFace
    for(i=0; i<4; i++) {
      faceDirs[i] = 0;
      faceEdges[i] = edges[3-i];
    }
    faces[0] = GImporter_createFace(importer,4,faceEdges,faceDirs,1,loopDef,planarSurface,0);

    // Now complete the model and delete the importer
    model = GImporter_complete(importer);
    auto isValid = GM_isValid(model,2,NULL);
    if(!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    }
    GImporter_delete(importer);
    
    cout<<"Number of vertices in model: "<<GM_numVertices(model)<<endl;
    cout<<"Number of edges in model: "<<GM_numEdges(model)<<endl;
    cout<<"Number of faces in model: "<<GM_numFaces(model)<<endl;
    cout<<"Number of regions in model: "<<GM_numRegions(model)<<endl;
    GM_write(model,modelFileName.c_str(),0,0);

    // This next section creates a surface mesh from the model.  You can comment out this section
    // if you don't want to mesh
    pMesh mesh = M_new(0,model);
    pACase meshCase = MS_newMeshCase(model);
    
    pModelItem domain = GM_domain(model);
    MS_setMeshSize(meshCase,domain,2,0.05,NULL);
    
    pSurfaceMesher surfMesh = SurfaceMesher_new(meshCase,mesh);
    SurfaceMesher_execute(surfMesh,progress);
    SurfaceMesher_delete(surfMesh);
    cout<<"Number of mesh faces in surface: "<<M_numFaces(mesh)<<endl;
    
    M_write(mesh,meshFileName.c_str(),0,progress);
    cout<<"Number of mesh regions in volume: "<<M_numRegions(mesh)<<endl;
    MS_deleteMeshCase(meshCase);
    M_release(mesh);
    // end of meshing section

    delete [] vertices;
    delete [] edges;
    delete [] faceEdges;
    delete [] faceDirs;
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
