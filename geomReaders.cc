#include "modelGen2d.h"
#include <set>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

struct HostGraph {
  Omega_h::HostRead<Omega_h::LO> values;
  Omega_h::HostRead<Omega_h::LO> offsets;
};

void markAsVisited(std::set<Omega_h::LO>& visitedEdges,
    Omega_h::LO edge) {
  assert(visitedEdges.count(edge) == 0);
  visitedEdges.insert(edge);
}

// get the edges adjacent to vtx (b) that have not been
//   visited and have class_dim == 1
// returns the index of the found edge or -1 otherwise
Omega_h::LO getNextEdge(Omega_h::LO edge, Omega_h::LO vtx,
    HostGraph& vtxToEdge, 
    Omega_h::HostRead<Omega_h::LO>& edgeToVtx,
    Omega_h::HostRead<Omega_h::I8>& edgeClassDim,
    std::set<Omega_h::LO>& visitedEdges) {
  Omega_h::LO nextEdge = -1;
  for (auto edgeIdx = vtxToEdge.offsets[vtx]; 
            edgeIdx < vtxToEdge.offsets[vtx + 1]; 
            ++edgeIdx) {
    const auto otherEdge = vtxToEdge.values[edgeIdx];
    const auto isDifferentEdge = (otherEdge != edge);
    const auto isNotVisited = (visitedEdges.count(otherEdge) == 0);
    const auto isClassifiedOnEdge = (edgeClassDim[otherEdge] == 1);
    if( isDifferentEdge && isNotVisited && isClassifiedOnEdge ) {
      nextEdge = otherEdge;
      break;
    }
  }
  return nextEdge;
}

Omega_h::LO findFirstEdge(Omega_h::HostRead<Omega_h::I8> classDim) {
  for(int i=0; i<classDim.size(); i++) {
    if(classDim[i] == 1) {
      return i;
    }
  }
  assert(false);
  return -1;
}

Omega_h::LO getDownVtx(Omega_h::HostRead<Omega_h::LO> edgeToVtx, Omega_h::LO edge, int which) {
  assert(which == 0 || which == 1);
  return edgeToVtx[edge*2+which];
}

Omega_h::LO getOtherVtx(Omega_h::HostRead<Omega_h::LO> edgeToVtx,
    Omega_h::LO edge, Omega_h::LO vtx) {
  auto vtxA = getDownVtx(edgeToVtx, edge, 0);
  auto vtxB = getDownVtx(edgeToVtx, edge, 1);
  assert( vtx == vtxA || vtx == vtxB );
  if( vtx != vtxA ) {
    return vtxA; 
  } else {
    return vtxB;
  }
}

void addVtx(GeomInfo& geom, Omega_h::HostRead<Omega_h::Real> coords, Omega_h::LO vtx) {
  assert(vtx >= 0);
  assert(vtx*2 < coords.size());
  assert(vtx*2+1 < coords.size());
  const auto x = coords[vtx*2];
  const auto y = coords[vtx*2+1];
  geom.addVtx(x,y);
}

GeomInfo readOmegahGeom(std::string fname, bool debug) {
  auto lib = Omega_h::Library();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(fname, lib.world(), &mesh);

  std::set<Omega_h::LO> visitedEdges;

  auto edgeClassDim_d = mesh.get_array<Omega_h::I8>(1, "class_dim");
  auto edgeToVtx_d = mesh.ask_down(1,0); 
  auto vtxToEdge_d = mesh.ask_up(0,1); 
  auto coords_d = mesh.coords();

  auto edgeClassDim = Omega_h::HostRead<Omega_h::I8>(edgeClassDim_d);
  auto edgeToVtx = Omega_h::HostRead<Omega_h::LO>(edgeToVtx_d.ab2b);
  HostGraph vtxToEdge{
    Omega_h::HostRead<Omega_h::LO>(vtxToEdge_d.ab2b),
    Omega_h::HostRead<Omega_h::LO>(vtxToEdge_d.a2ab)
  };
  auto coords = Omega_h::HostRead<Omega_h::Real>(coords_d);

  //HACK - the following won't work for all meshes
  GeomInfo geom; 
  geom.numVtx = geom.numEdges = 0;
  auto edge = findFirstEdge(edgeClassDim);
  auto vtx = getDownVtx(edgeToVtx, edge, 0);
  geom.firstContourPt = 0; //using geom indexing, not omegah's 
  const auto firstVtx = vtx;
  addVtx(geom, coords, vtx);
  while( -1 != (edge = getNextEdge(edge, vtx, vtxToEdge, edgeToVtx, edgeClassDim, visitedEdges)) ) {
    markAsVisited(visitedEdges, edge);
    auto nextVtx = getOtherVtx(edgeToVtx, edge, vtx);
    if(nextVtx != firstVtx) { //don't add the first vtx twice
      addVtx(geom, coords, nextVtx);
    }
    geom.addEdge(vtx, nextVtx);
    vtx = nextVtx;
  }
  return geom;
}
