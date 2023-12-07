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
    geom.vtx_x[i] = pt[0];
    geom.vtx_y[i] = pt[1];
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
    geom.edges[i] = readEdge(mshFile,debug);
    if(debug) std::cout << "edge " << geom.edges[i][0] << ", " << geom.edges[i][1] << std::endl;
  }

  return geom;
}

int main(int argc, char **argv)
{

  pGImporter importer;  // the importer object used to create the geometry
  pGVertex* vertices;    // array to store the returned model vertices
  pGEdge* edges;     // array to store the returned model edges
  pGFace faces[6];      // array to store the returned model faces
  pGRegion region;      // pointer to returned model region
  pGModel model;        // pointer to the complete model

  auto geom = readJigGeom(argv[1]);

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
    for(i=0; i<geom.numVtx; i++) {
      double vtx[3] = {geom.vtx_x[i], geom.vtx_y[i], 0};
      vertices[i] = GImporter_createVertex(importer, vtx);
    }

    // Now we'll add the edges
    double point0[3],point1[3];  // xyz locations of the two vertices
    pCurve linearCurve;

    // First, the bottom edges at z=0, connecting the first four vertices in the array
    // 0->1, 1->2, 2->3, 3->0  (indices of the vertices)
    for(i=0; i<geom.numEdges; i++) {
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
    // There is only one loop per surface for the faces of the box
    // When defining the loop, will always start with the first edge in the faceEdges array
    int loopDef[1] = {0};        
    pSurface planarSurface;

    // First the face between the bounding box and the grounding line
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
    return 0;
/*
    // Create the face 
    faceEdges = new pGEdge[geom.numEdges];
    faceDirs = new int[geom.numEdges];
    for(i=0; i<geom.numEdges; i++) {
      faceDirs[i] = 0;
      faceEdges[i] = edges[3-i]; // edge order 3->2->1->0
    }
    faces[0] = GImporter_createFace(importer,4,faceEdges,faceDirs,1,loopDef,planarSurface,1);

    // Now the side faces of the box - each side face has the edges defined in the same way
    // for the first side face, the edge order is 0->5->8->4
    for(i=0; i<4; i++) {
      //Define surface such that normals all point out of the box
      for(int j=0; j<3; j++) {
        corner[j] = vert_xyz[i][j];      // the corner is the lower left vertex location
        xPt[j] = vert_xyz[(i+1)%4][j];   // the xPt the lower right vertex location
        yPt[j] = vert_xyz[i+4][j];       // the yPt is the upper left vertex location
      }
      planarSurface = SSurface_createPlane(corner,xPt,yPt);

      faceEdges[0] = edges[i];
      faceDirs[0] = 1;
      faceEdges[1] = edges[(i+1)%4+4];
      faceDirs[1] = 1;
      faceEdges[2] = edges[i+8];
      faceDirs[2] = 0;
      faceEdges[3] = edges[i+4];
      faceDirs[3] = 0;

      faces[i+1] = GImporter_createFace(importer,4,faceEdges,faceDirs,1,loopDef,planarSurface,1);
    }

    // Finally the top face of the box
    // Define the surface - we want the normal to point out of the box
    for(i=0; i<3; i++) {
      corner[i] = vert_xyz[4][i];  // the corner is at {0,0,10}
      xPt[i] = vert_xyz[5][i];     // the xPt is at {10,0,10}
      yPt[i] = vert_xyz[7][i];     // the yPt is at {0,10,10}
    }
    planarSurface = SSurface_createPlane(corner,xPt,yPt);
    // Define and insert the face
    for(i=0; i<4; i++) {
      faceDirs[i] = 1;
      faceEdges[i] = edges[i+8]; // edge order 8->9->10->11
    }
    faces[5] = GImporter_createFace(importer,4,faceEdges,faceDirs,1,loopDef,planarSurface,1);
    
    // Create the region of the box 
  
    // the side of the face used by the region is stored in regionDirs.  The normal of the 
    // face points out of side 1, so since we created the normals to point away from the box, 
    // the region uses side 0 of each face.
    int regionDirs[6] = {0,0,0,0,0,0};
    int shellIndex[1] = {0};
    
    region = GImporter_createRegion(importer,6,faces,regionDirs,1,shellIndex);
    
    // Now complete the model and delete the importer
    model = GImporter_complete(importer);
    GImporter_delete(importer);
    
    cout<<"Number of vertices in model: "<<GM_numVertices(model)<<endl;
    cout<<"Number of edges in model: "<<GM_numEdges(model)<<endl;
    cout<<"Number of faces in model: "<<GM_numFaces(model)<<endl;
    cout<<"Number of regions in model: "<<GM_numRegions(model)<<endl;
    GM_write(model,"model.smd",0,0);

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
    
    pVolumeMesher volMesh = VolumeMesher_new(meshCase,mesh);
    VolumeMesher_execute(volMesh,progress);
    VolumeMesher_delete(volMesh);
    M_write(mesh,"volume.sms",0,progress);
    cout<<"Number of mesh regions in volume: "<<M_numRegions(mesh)<<endl;
    MS_deleteMeshCase(meshCase);
    M_release(mesh);
    // end of meshing section
*/

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
