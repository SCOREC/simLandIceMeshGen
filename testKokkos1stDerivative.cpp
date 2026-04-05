#include <Kokkos_Core.hpp>
#include "BSplineKokkos.h"
#include "BSpline.h"
#include <vector>

//These imports are from testSplineFitting.cc
#include "splineInterpolation.h"
#include "curveReader.h"
#include <cassert>
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <numeric>


template<typename MemSpace>
void printDView(Kokkos::View<double*, MemSpace>& dView) {
    auto host_dView = Kokkos::create_mirror_view(dView);
    Kokkos::deep_copy(host_dView, dView);
    for (int i = 0; i < host_dView.extent(0); i++) {
        std::cout << host_dView(i) << std::endl;
    }
}

template<typename MemSpace>
void printIView(Kokkos::View<int*, MemSpace>& iView) {
    auto host_iView = Kokkos::create_mirror_view(iView);
    Kokkos::deep_copy(host_iView, iView);
    for (int i = 0; i < host_iView.extent(0); i++) {
        std::cout << host_iView(i) << std::endl;
    }
}

template <typename MemSpace>
void printSpline(BSplineKokkos<MemSpace> spline) {
    std::cout << "Spline content" << std::endl;
    Kokkos::View<double*, MemSpace> dView = spline.getCtrlPts();
    std::cout << "CtrlPts" << std::endl;
    printDView(dView);
    Kokkos::View<int*, MemSpace> iView = spline.getCPOffset();
    std::cout << "CtrlPts Offset" << std::endl;
    printIView(iView);
    dView = spline.getKnots();
    std::cout << "Knots" << std::endl;
    printDView(dView);
    iView = spline.getKnotsOffset();
    std::cout << "Knots Offset" << std::endl;
    printIView(iView);
    iView = spline.getOrder();
    std::cout << "Order" << std::endl;
    printIView(iView);
}

void printVector(const std::string& name, std::vector<double>& vector) {
    std::cout << name << ": " << vector.size() << std::endl;
    for (int i = 0; i < vector.size(); i++) {
        std::cout << vector[i] << std::endl;
    }
}

double epsilon = 0.005;

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        
	//Check for what memory space we should use
        //Sticking to just cuda space or host space for now

        #ifdef KOKKOS_ENABLE_CUDA
        #define MemSpace Kokkos::CudaSpace
        #endif
        #ifndef MemSpace
        #define MemSpace Kokkos::HostSpace
        #endif

        using ExecutionSpace = MemSpace::execution_space;


        if (argc != 3) {                                                    //Check for the arguments needed                                std::cerr<<"Input arguments: <input csv file> <expected curve length>" << std::endl;                                            std::cerr << "input csv need these columns: ";
            std::cerr << "x,y,z,isOnCurve,angle,isMdlVtx" << std::endl;
            return 1;                                                   }
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
        }                                                               else {
            serialBSP = SplineInterp::fitCubicSplineToPoints(curve.x, curve.y);
        }
	
	std::vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
	int order;

	//Get the info from serial spline, we will feed this to kokkos spline
	serialBSP.x.getpara(order, ctrlPtsX, knots, weight);
	serialBSP.y.getpara(order, ctrlPtsY, knots, weight);

	BSplineKokkos<ExecutionSpace> kokkosBSP(order, ctrlPtsX, ctrlPtsY, knots);
        //printSpline(kokkosBSP);
	

	/*INITALIZATION COMPLETE, TEST STARTS BELOW*/
	//TEST 1 HERE
	 std::cout << "TEST 1, EVAL AT 0" << std::endl;
	 double derivX = serialBSP.x.evalFirstDeriv(0);
	 double derivY = serialBSP.y.evalFirstDeriv(0);
	 std::cout << "SERIAL: 1st derivative, input x = 0: " << derivX << ", " << derivY << std::endl;

	 std::vector<double> kokkos1stDeriv = kokkosBSP.eval1stDeriv(0, 0);
	 double kokkosDerivX = kokkos1stDeriv[0];
	 double kokkosDerivY = kokkos1stDeriv[1];

	 double derivXDiff = derivX - kokkosDerivX;
	 double derivYDiff = derivY - kokkosDerivY;
	 std::cout << "KOKKOS: 1st derivative, input x = 0: " << kokkosDerivX <<", " << kokkosDerivY  << std::endl;
	
	 std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;

	
	//TEST 2 HERE
	std::cout << "TEST 2, EVAL AT 0.2" << std::endl;
	derivX = serialBSP.x.evalFirstDeriv(0.2);
	derivY = serialBSP.y.evalFirstDeriv(0.2);
        std::cout << "SERIAL: 1st derivative, input x = 0.2: " << derivX << ", " << derivY << std::endl;
	kokkos1stDeriv = kokkosBSP.eval1stDeriv(0.2, 0);
	kokkosDerivX = kokkos1stDeriv[0];
	kokkosDerivY = kokkos1stDeriv[1];

	derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 0.2: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;


	//TEST 3
	std::cout << "TEST 3: EVAL AT 0.41" << std::endl;
	derivX = serialBSP.x.evalFirstDeriv(0.41);
        derivY = serialBSP.y.evalFirstDeriv(0.41);
	std::cout << "SERIAL: 1st derivative, input x = 0.41: " << derivX << ", " << derivY << std::endl;

        kokkos1stDeriv = kokkosBSP.eval1stDeriv(0.41, 0);
        kokkosDerivX = kokkos1stDeriv[0];
        kokkosDerivY = kokkos1stDeriv[1];
        
	derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 0.41: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;


	//TEST 4
	std::cout << "TEST 4: EVAL AT 0.5" << std::endl;
        derivX = serialBSP.x.evalFirstDeriv(0.5);
        derivY = serialBSP.y.evalFirstDeriv(0.5);
        std::cout << "SERIAL: 1st derivative, input x = 0.5: " << derivX << ", " << derivY << std::endl;

        kokkos1stDeriv = kokkosBSP.eval1stDeriv(0.5, 0);
        kokkosDerivX = kokkos1stDeriv[0];
        kokkosDerivY = kokkos1stDeriv[1];

        derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 0.5: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;

	//TEST 5
	std::cout << "TEST 5: EVAL AT 0.66" << std::endl;
        derivX = serialBSP.x.evalFirstDeriv(0.66);
        derivY = serialBSP.y.evalFirstDeriv(0.66);
        std::cout << "SERIAL: 1st derivative, input x = 0.66: " << derivX << ", " << derivY << std::endl;
                                                
	kokkos1stDeriv = kokkosBSP.eval1stDeriv(0.66, 0);
        kokkosDerivX = kokkos1stDeriv[0];
        kokkosDerivY = kokkos1stDeriv[1];

        derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 0.66: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;


	//TEST 6
	std::cout << "Test 6: EVAL AT 0.73" << std::endl;
        derivX = serialBSP.x.evalFirstDeriv(0.73);
        derivY = serialBSP.y.evalFirstDeriv(0.73);
        std::cout << "SERIAL: 1st derivative, input x = 0.73: " << derivX << ", " << derivY << std::endl;
                                                                        kokkos1stDeriv = kokkosBSP.eval1stDeriv(0.73, 0);
        kokkosDerivX = kokkos1stDeriv[0];
        kokkosDerivY = kokkos1stDeriv[1];

        derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 0.73: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;


	//TEST 7 
	std::cout << "TEST 7: EVAL AT 0.75" << std::endl;
        derivX = serialBSP.x.evalFirstDeriv(0.75);
        derivY = serialBSP.y.evalFirstDeriv(0.75);
        std::cout << "SERIAL: 1st derivative, input x = 0.75: " << derivX << ", " << derivY << std::endl;
                                                                        kokkos1stDeriv = kokkosBSP.eval1stDeriv(0.75, 0);
        kokkosDerivX = kokkos1stDeriv[0];
        kokkosDerivY = kokkos1stDeriv[1];

        derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 0.75: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;


	//TEST 8
	std::cout << "TEST 8: EVAL AT 0.89" << std::endl;
        derivX = serialBSP.x.evalFirstDeriv(0.89);
        derivY = serialBSP.y.evalFirstDeriv(0.89);
        std::cout << "SERIAL: 1st derivative, input x = 0.89: " << derivX << ", " << derivY << std::endl;
                          
	kokkos1stDeriv = kokkosBSP.eval1stDeriv(0.89, 0);
        kokkosDerivX = kokkos1stDeriv[0];
        kokkosDerivY = kokkos1stDeriv[1];

        derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 0.89: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;

	
	//TEST 9 
	std::cout << "TEST 9: EVAL AT 0.94" << std::endl;
        derivX = serialBSP.x.evalFirstDeriv(0.94);
        derivY = serialBSP.y.evalFirstDeriv(0.94);
        std::cout << "SERIAL: 1st derivative, input x = 0.94: " << derivX << ", " << derivY << std::endl;
                                                                        kokkos1stDeriv = kokkosBSP.eval1stDeriv(0.94, 0);
        kokkosDerivX = kokkos1stDeriv[0];
        kokkosDerivY = kokkos1stDeriv[1];

        derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 0.94: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;


	//TEST 10
	std::cout << "TEST 10: EVAL AT 1" << std::endl;
        derivX = serialBSP.x.evalFirstDeriv(1);
        derivY = serialBSP.y.evalFirstDeriv(1);
        std::cout << "SERIAL: 1st derivative, input x = 1: " << derivX << ", " << derivY << std::endl;

	kokkos1stDeriv = kokkosBSP.eval1stDeriv(1, 0);
        kokkosDerivX = kokkos1stDeriv[0];
        kokkosDerivY = kokkos1stDeriv[1];

        derivXDiff = derivX - kokkosDerivX;
        derivYDiff = derivY - kokkosDerivY;

        std::cout << "KOKKOS: 1st derivative, input x = 1: " << kokkosDerivX << ", " << kokkosDerivY << std::endl;

        std::cout << "Serial and Kokkos Difference: " << derivXDiff << ", " << derivYDiff << std::endl;

    }
    Kokkos::finalize();
}
