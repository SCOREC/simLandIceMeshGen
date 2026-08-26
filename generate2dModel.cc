#include "simModelGen2d.h"
#include "netcdfWriter.h"
#include <algorithm> //std::sort
#include <numeric> //std::accumulate
#include "Omega_h_file.hpp"
#include "Omega_h_library.hpp"

void messageHandler(int type, const char *msg);

std::string getFileExtension(const std::string &filename) {
  size_t dotPos = filename.rfind('.');
  if (dotPos != std::string::npos) {
    return filename.substr(dotPos);
  }
  return "";
}

// A single contour supplied via a repeated --contour command line flag, e.g.:
//   --contour file=inner.msh,order=0,units=km,fail-if-cleaned
struct ContourSpec {
  std::string file;
  int order = -1;
  std::string units;
  bool failIfCleaned = false;
  bool boundaryTriangles = false;
};

std::vector<std::string> splitOn(const std::string &s, char delim) {
  std::vector<std::string> parts;
  size_t start = 0;
  while (start <= s.size()) {
    size_t pos = s.find(delim, start);
    if (pos == std::string::npos) {
      parts.push_back(s.substr(start));
      break;
    }
    parts.push_back(s.substr(start, pos - start));
    start = pos + 1;
  }
  return parts;
}

ContourSpec parseContourSpec(const std::string &arg) {
  ContourSpec spec;
  for (const auto &kv : splitOn(arg, ',')) {
    if (kv == "fail-if-cleaned") {
      spec.failIfCleaned = true;
      continue;
    }
    if (kv == "boundary-triangles") {
      spec.boundaryTriangles = true;
      continue;
    }
    auto eq = kv.find('=');
    if (eq == std::string::npos) {
      std::cerr << "ERROR: malformed --contour option '" << kv << "'\n";
      exit(EXIT_FAILURE);
    }
    const auto key = kv.substr(0, eq);
    const auto value = kv.substr(eq + 1);
    if (key == "file") {
      spec.file = value;
    } else if (key == "order") {
      spec.order = std::stoi(value);
    } else if (key == "units") {
      spec.units = value;
    } else {
      std::cerr << "ERROR: unknown --contour option '" << key << "'\n";
      exit(EXIT_FAILURE);
    }
  }
  if (spec.file.empty()) {
    std::cerr << "ERROR: --contour is missing required 'file=' option\n";
    exit(EXIT_FAILURE);
  }
  if (spec.order < 0) {
    std::cerr << "ERROR: --contour '" << spec.file
               << "' is missing required 'order=' option (0=innermost)\n";
    exit(EXIT_FAILURE);
  }
  if (spec.units != "m" && spec.units != "km") {
    std::cerr << "ERROR: --contour '" << spec.file
               << "' has missing/invalid 'units=' option (must be 'm' or 'km')\n";
    exit(EXIT_FAILURE);
  }
  return spec;
}

// Reads a single contour from a file. The file must contain exactly one
// contour; if the legacy bbox+inner splitting heuristic in
// splitIntoInnerAndOuter fires (i.e. the file looks like it defines two
// contours) this is treated as an error since it's ambiguous which contour
// the caller wants.
GeomInfo readSingleContour(const std::string &filename, bool boundaryTriangles, bool debug) {
  const auto ext = getFileExtension(filename);
  ModelFeatures features;
  if (ext == ".vtk") {
    features = readVtkGeom(filename, boundaryTriangles, debug);
  } else if (ext == ".msh") {
    if (boundaryTriangles) {
      std::cerr << "ERROR: '" << filename << "' the 'boundary-triangles' "
                   "contour flag is only supported for .vtk files\n";
      exit(EXIT_FAILURE);
    }
    features = readJigGeom(filename, debug);
  } else {
    std::cerr << "Unsupported file extension: " << ext << "\n";
    exit(EXIT_FAILURE);
  }
  const bool hasInner = features.inner.numVtx > 0;
  const bool hasOuter = features.outer.numVtx > 0;
  if (hasInner && hasOuter) {
    std::cerr << "ERROR: '" << filename << "' contains two contours "
                 "(bounding box + inner contour) but --contour expects "
                 "exactly one contour per file\n";
    exit(EXIT_FAILURE);
  }
  return hasInner ? features.inner : features.outer;
}

Omega_h::Reals setParametricCoords(GeomInfo& geom, const PointClassification& ptClass, const SplineInterp::SplineInfo& sinfo, bool debug=false) {
  const int numVtx = geom.numVtx;
  assert(ptClass.splineIdx.size() == numVtx);
  Omega_h::HostWrite<Omega_h::Real> paraCoords(geom.numVtx*2);
  if(debug)
    std::cerr << __func__ << " sIdx, x, y, paraX, paraY\n";
  for(int i=0; i<geom.numVtx; i++) {
    const auto sIdx = ptClass.splineIdx.at(i);
    const auto bspline = sinfo.splines.at(sIdx);
    const auto x = geom.vtx_x.at(i); 
    const auto y = geom.vtx_y.at(i); 
    const auto noGuess = -1;
    if(debug)
      std::cerr << __func__ << " " << sIdx << " " << x << ", " << y << '\n';
    const double t = bspline.invEval({x,y}, noGuess, debug);
    paraCoords[i*2] = t;
    paraCoords[i*2+1] = t;
    if(debug)
      std::cerr << __func__ << " " << paraCoords[i*2] << ", " << paraCoords[i*2+1] << "\n";
    assert(!std::isnan(paraCoords[i*2]));
    assert(!std::isnan(paraCoords[i*2+1]));
  }
  auto paraCoords_d = Omega_h::read(paraCoords.write());
  return paraCoords_d;
}

void writePointParametricCoords(Omega_h::Reals paraCoords_d, std::string filename) {
  std::ofstream file(filename);
  assert(file.is_open());
  const int compressed = 0;
  //the following is from src/Omega_h_file.cpp write(...)
  unsigned char const magic[2] = {0xa1, 0x1a};
  file.write(reinterpret_cast<const char*>(magic), sizeof(magic));
  bool needs_swapping = !Omega_h::is_little_endian_cpu();
  Omega_h::binary::write_value(file, compressed, needs_swapping);
  Omega_h::binary::write_array(file, paraCoords_d, compressed, needs_swapping);
  file.close();
}

int main(int argc, char **argv) {
  std::vector<ContourSpec> contourSpecs;
  std::vector<std::string> positional;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--contour") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --contour requires an argument\n";
        return 1;
      }
      contourSpecs.push_back(parseContourSpec(argv[++i]));
    } else {
      positional.push_back(arg);
    }
  }

  const size_t numExpectedPositionalArgs = 5;
  if (positional.size() != numExpectedPositionalArgs || contourSpecs.empty()) {
    std::cerr << "Usage: --contour file=<jigsaw .msh or .vtk file>,order=<N>,units=<m|km>"
                 "[,fail-if-cleaned][,boundary-triangles] [--contour ...] "
                 "<output prefix> <coincidentVtxTolerance> <angleTolerance> "
                 "<onCurveAngleTolerance> <createMesh>\n";
    std::cerr << "--contour specifies one nested contour and may be repeated. "
                 "order is the nesting order of the contour, 0=innermost, "
                 "N-1=outermost, where N is the total number of contours. "
                 "Exactly one or two --contour flags are currently supported.\n";
    std::cerr << "  file: path to a jigsaw .msh or .vtk file containing exactly "
                 "one contour\n";
    std::cerr << "  order: nesting order, 0=innermost, N-1=outermost\n";
    std::cerr << "  units: m=meters, km=kilometers\n";
    std::cerr << "  fail-if-cleaned: exit with error if cleaning removes any "
                 "input points from this contour\n";
    std::cerr << "  boundary-triangles: the .vtk file's cells are boundary "
                 "triangles (VTK POLYGONS) rather than boundary edges (VTK "
                 "LINES), and POINT_DATA gives each vertex's 0-based position "
                 "in the CW/CCW boundary traversal (-1 for interior points). "
                 "Only supported for .vtk files.\n";
    std::cerr << "coincidentVtxTolerance is the mininum allowed "
                 "distance between adjacent vertices in the "
                 "input.  Vertices within the specified distance will "
                 "be merged.\n";
    std::cerr << "angleTolerance defines the upper bound and "
                 "-angleTolerance defines lower bound for determining "
                 "if a vertex bounding two consecutative edges should be "
                 "treated as a model vertex.\n";
    std::cerr << "onCurveAngleTolerance defines the upper bound on the angle "
                 "between adjacent edges in a sequence of four consecutive edges "
                 "used to determine if they are part of the same curve.\n";
    std::cerr << "createMesh = 1:generate mesh, otherwise, "
                 "skip mesh generation.\n";
    return 1;
  }

  const auto prefix = positional[0];
  const auto coincidentPtTol = std::stof(positional[1]);
  const auto angleTol = std::atof(positional[2].c_str());
  const auto onCurveAngleTol = std::atof(positional[3].c_str());
  const bool doCreateMesh = (std::stoi(positional[4]) == 1);

  const size_t numContours = contourSpecs.size();
  if (numContours > 2) {
    std::cerr << "ERROR: " << numContours << " --contour flags given; only "
                 "one or two nested contours are currently supported\n";
    return 1;
  }
  {
    std::vector<bool> orderSeen(numContours, false);
    for (const auto &spec : contourSpecs) {
      if (spec.order < 0 || static_cast<size_t>(spec.order) >= numContours ||
          orderSeen[spec.order]) {
        std::cerr << "ERROR: --contour order values must be exactly "
                     "{0.."
                  << (numContours - 1) << "} with no gaps or duplicates\n";
        return 1;
      }
      orderSeen[spec.order] = true;
    }
  }
  std::sort(contourSpecs.begin(), contourSpecs.end(),
      [](const ContourSpec &a, const ContourSpec &b) { return a.order < b.order; });

  std::cout << "output prefix: " << prefix << " "
            << "coincidentPtTol: " << coincidentPtTol << " "
            << "angleTol: " << angleTol << " "
            << "onCurveAngleTol: " << onCurveAngleTol << " "
            << "createMesh: " << doCreateMesh << "\n";
  for (const auto &spec : contourSpecs) {
    std::cout << "contour order " << spec.order << ": file=" << spec.file
              << " units=" << spec.units
              << " fail-if-cleaned=" << spec.failIfCleaned
              << " boundary-triangles=" << spec.boundaryTriangles << "\n";
  }

  const auto debug = true;
  const double coincidentPtTolSquared = coincidentPtTol*coincidentPtTol;

  //load, clean, and orient each contour, in nesting order (0=innermost)
  std::vector<GeomInfo> contours;
  contours.reserve(numContours);
  for (const auto &spec : contourSpecs) {
    auto geom = readSingleContour(spec.file, spec.boundaryTriangles, debug);
    //simmetrix operations are done in km to avoid problems with floating
    //point operations
    if (spec.units == "m") {
      convertMetersToKm(geom);
    }
    const int preCleanNumVtx = geom.numVtx;
    //force the contour to be positive (CCW)
    if( !spec.boundaryTriangles )
      geom = cleanGeom(geom, coincidentPtTolSquared, debug);
    makeOrientationPositive(geom);
    if (spec.failIfCleaned && geom.numVtx < preCleanNumVtx) {
      std::cerr << "ERROR: cleaning removed " << (preCleanNumVtx - geom.numVtx)
                << " point(s) from contour '" << spec.file << "'\n";
      return 1;
    }
    contours.push_back(std::move(geom));
  }

  //TODO: generalize createEdges/createFaces/createMesh to support an
  //arbitrary number of nested contours. For now only one or two are
  //supported: a single contour is the domain boundary (features.outer);
  //two contours are the domain boundary (outermost, features.outer) and
  //a single interior contour (innermost, features.inner).
  ModelFeatures features;
  if (numContours == 1) {
    features.outer = contours[0];
  } else {
    features.inner = contours[0];
    features.outer = contours[1];
  }

  std::string modelFileName = prefix + ".smd";
  std::string meshFileName = prefix + ".sms";


  // You will want to place a try/catch around all SimModSuite calls,
  // as errors are thrown.
  try {
    auto ohLib = Omega_h::Library(&argc, &argv);
    Sim_logOn("simMeshGen.log");
    SimModel_start(); // Call before Sim_readLicenseFile
    // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
    // pass in the location of a file containing your keys.  For a release
    // product, use Sim_registerKey()
    Sim_readLicenseFile(0);
    // Tessellation of GeomSim geometry requires Meshing to have started
    MS_init();

    Sim_setMessageHandler(messageHandler);
    pProgress progress = NULL;
    if(debug) {
      progress = Progress_new();
      Progress_setDefaultCallback(progress);
    }

    ModelTopo mdlTopo;
    mdlTopo.model = GM_new(1);
    mdlTopo.part = GM_rootPart(mdlTopo.model);
    mdlTopo.region = GIP_outerRegion(mdlTopo.part);

    auto [isPointOnCurveInner, isMdlVtxInner] = discoverTopology(features.inner, coincidentPtTolSquared, angleTol, onCurveAngleTol, debug);
    auto [isPointOnCurveOuter, isMdlVtxOuter] = discoverTopology(features.outer, coincidentPtTolSquared, angleTol, onCurveAngleTol, debug);

    const auto numInnerMdlVerts = isMdlVtxInner.size() ? std::accumulate(isMdlVtxInner.begin(), isMdlVtxInner.end(), 0) : 0;
    const auto numOuterMdlVerts = isMdlVtxOuter.size() ? std::accumulate(isMdlVtxOuter.begin(), isMdlVtxOuter.end(), 0) : 0;
    auto splinesInner = SplineInterp::SplineInfo(numInnerMdlVerts);
    auto splinesOuter = SplineInterp::SplineInfo(numOuterMdlVerts);
    PointClassification ptClassInner(features.inner.numVtx);
    PointClassification ptClassOuter(features.outer.numVtx);
    createEdges(mdlTopo, features.outer, ptClassOuter, splinesOuter, isPointOnCurveOuter, isMdlVtxOuter, debug);
    createEdges(mdlTopo, features.inner, ptClassInner, splinesInner, isPointOnCurveInner, isMdlVtxInner, debug);

    const auto paraCoordsOuter = setParametricCoords(features.outer, ptClassOuter, splinesOuter);
    const auto paraCoordsInner = setParametricCoords(features.inner, ptClassInner, splinesInner);

    writePointParametricCoords(paraCoordsInner, prefix + "_parametricInner.oshb");
    writePointParametricCoords(paraCoordsOuter, prefix + "_parametricOuter.oshb");
    //write the point classification to an omegah binary file
    ptClassInner.writeToOsh(prefix + "_classInner.oshb");
    ptClassOuter.writeToOsh(prefix + "_classOuter.oshb");
    //write the bsplines to an omegah binary file
    splinesInner.writeToOsh(prefix + "_splinesInner.oshb");
    splinesOuter.writeToOsh(prefix + "_splinesOuter.oshb");
    //write the sampled bsplines to a csv file
    splinesInner.writeSamplesToCsv(prefix + "_splinesInner.csv");
    splinesOuter.writeSamplesToCsv(prefix + "_splinesOuter.csv");

    auto planeBounds = getBoundingPlane(features.outer);
    const bool hasSingleContour = (numContours == 1);
    createFaces(mdlTopo, planeBounds, hasSingleContour, debug);

    printModelInfo(mdlTopo.model);

    GM_write(mdlTopo.model, modelFileName.c_str(), 0, 0);

    auto isValid = GM_isValid(mdlTopo.model, 2, NULL);
    if (!isValid) {
      fprintf(stderr, "ERROR: model is not valid... exiting\n");
      exit(EXIT_FAILURE);
    }

    if(doCreateMesh) {
      auto mesh = createMesh(mdlTopo, meshFileName, progress);
      std::string netcdfFileName = prefix + ".nc";
      //compass assumes units of meters, need to convert back to meters,
      //simmetrix operations are done in units of km
      const auto convertBackToMeters = true;
      writeMeshSimToNetCDF(mesh, mdlTopo.model, netcdfFileName, convertBackToMeters);
      M_release(mesh);
    }

    // cleanup
    GM_release(mdlTopo.model);
    if(debug) Progress_delete(progress);
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
