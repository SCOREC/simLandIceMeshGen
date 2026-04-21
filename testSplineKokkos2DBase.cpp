//Constructor test for new representation of Kokkos spline
//We will be comparing this against the serial version
#include <Kokkos_Core.hpp>
#include "BSplineKokkos2D.h"
#include "BSpline.h"
#include "splineInterpolation.h"
#include "curveReader.h"

#include <vector>
#include <iostream>
#include <string>
#include <cassert>
#include <math.h>
#include <fstream>
#include <numeric>

int main(int argc, char* argv[]) {
    //We check how many arguments are given
    int retVal = 0;
    if (argc != 3) {
        std::cerr<< "Input arguments: <input csv file> <expected curve length>" << std::endl;
	std::cerr << "input csv need these columns: ";
	std::cerr << "coordinate x, coordinate y, coordinate z,isOnCurve,angle,isMdlVtx" << std::endl;
	return 1;
    }
    Kokkos::initialize(argc, argv); 
    {
        
        #ifdef KOKKOS_ENABLE_CUDA
        #define MemSpace Kokkos::CudaSpace
        #endif
        #ifndef MemSpace
        #define MemSpace Kokkos::HostSpace
        #endif
	
	using ExecutionSpace = MemSpace::execution_space;
	std::string inputCSV = argv[1];
	int extensionPos = inputCSV.rfind(".");
        int slashPos = inputCSV.rfind("/");
        std::string fileNameNoExt = inputCSV.substr(slashPos+1, extensionPos);
        double expectedCurveLength = std::stod(argv[2]);
        auto curve = CurveReader::readCurveInfo(inputCSV);

        //Construct BSpline2d object
        SplineInterp::BSpline2d serialBSP;
        if (curve.x.size() == 2) {
            serialBSP = SplineInterp::attach_piecewise_linear_curve(curve.x, curve.y);
        } else {
            serialBSP = SplineInterp::fitCubicSplineToPoints(curve.x, curve.y);
        }

	std::vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
	int order;

	serialBSP.x.getpara(order, ctrlPtsX, knots, weight);
	serialBSP.y.getpara(order, ctrlPtsY, knots, weight);

	BSplineKokkos2D<ExecutionSpace> kokkosBSP(order, ctrlPtsX, ctrlPtsY, knots);

    }
    Kokkos::finalize();
    return retVal;
}

