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
	//We will try to test with dummy data for BSplineKokkos creation
	std::vector<double> dummyCtrlPtsX{2.14, 3.98, 1.56, 9.10, 11.87};

	std::vector<double> dummyCtrlPtsY{7.928, 12.18, 9.17, 6.67, 12.78};
	std::vector<double> dummyKnotsX{5.99, 8.98, 17.71, 21.1, 2.154, 9.17, 6.32};
	//std::vector<double> dummyKnotsY{2.85, 1.16, 0.82, 1.89, 9.82, 12.0, 7.19};
	int dOrder = 2;
	//std::vector<double> dummyWeights{0.87, 0.125, 0.03};
	//Initialize 1 Kokkos Spline
	BSplineKokkos<ExecutionSpace> ex1(dOrder, dummyCtrlPtsX, dummyCtrlPtsY, dummyKnotsX);

	//Check if we actually initialized with correct values
	std::cout << "ex1" << std::endl;
	printSpline(ex1);
	std::vector<BSplineKokkos<ExecutionSpace>> splineList;
	splineList.push_back(ex1);

	//Test constructor that takes in multiple BSplineKokkos objects
	dummyCtrlPtsX.assign({1.20, 12.91, 7.11});
	dummyCtrlPtsY.assign({3.31, 7.89, 9.12});
	dummyKnotsX.assign({2.18, 5.56, 3.18, 16.15});
	//dummyKnotsY.assign({1.971,4.12, 8.192, 12.11});
	dOrder = 3;
	//dummyWeights.assign({0.19, 0.27});
	BSplineKokkos<ExecutionSpace> ex2(dOrder, dummyCtrlPtsX, dummyCtrlPtsY, dummyKnotsX);
	splineList.push_back(ex2);
	std::cout << "ex2" << std::endl;
	printSpline(ex2);

	dummyCtrlPtsX.assign({18.29, 12.12, 5.92, 7.182, 1.16});
	dummyCtrlPtsY.assign({4.72, 5.78, 9.12, 1.29, 10.2});
	dummyKnotsX.assign({4.34, 29.3, 18.7, 15.14, 16.87, 9.99});
	dOrder = 5;
	std::cout << "\t\tSTART OF MULTISPLINE CREATION" << std::endl;
	BSplineKokkos<ExecutionSpace> ex3(dOrder, dummyCtrlPtsX, dummyCtrlPtsY, dummyKnotsX);
	splineList.push_back(ex3);
	std::cout << "ex3" << std::endl;
	printSpline(ex3);

	BSplineKokkos<ExecutionSpace> fromList(splineList);

	std::cout << "Spline created from list of splines" << std::endl;
	printSpline(fromList);

	std::cout << "\n\nChecking derivative calculation against serial version" << std::endl;
	std::cout << "Number of commandline arguments: " << argc << std::endl;
	for (int i = 0; i < argc; i++) {
	    std::cout << argv[i] << std::endl;
	}
	if (argc != 3) {
	    //Check for the arguments needed
	    std::cerr<<"Input arguments: <input csv file> <expected curve length>" << std::endl;
	    std::cerr << "input csv need these columns: ";
	    std::cerr << "x,y,z,isOnCurve,angle,isMdlVtx" << std::endl;
	    return 1;
	}

	//Set up before using curve reader
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
	}
	else {
	    serialBSP = SplineInterp::fitCubicSplineToPoints(curve.x, curve.y);
	}
	
	//Get the parameters and print the spline we initialized
	std::vector<double> ctrlPtsX, ctrlPtsY, knots, weight;
	int order;
	serialBSP.x.getpara(order, ctrlPtsX, knots, weight);
	std::cout << "x component of the spline" << std::endl;
	std::cout << "order: " << order << std::endl;
	printVector("ctrlPtsX", ctrlPtsX);
	printVector("knots", knots);
	printVector("weight", weight);
	serialBSP.y.getpara(order, ctrlPtsY, knots, weight);
        std::cout << "y component of the spline" << std::endl;
        std::cout << "order: " << order << std::endl;
        printVector("ctrlPtsY", ctrlPtsY);
        printVector("knots", knots);
        printVector("weight", weight);

	//Now that we ensured that serial BSpline is created correctly we will now feed these data to our BSpline Kokkos Object
	BSplineKokkos<ExecutionSpace> kokkosBSP(order, ctrlPtsX, ctrlPtsY, knots);
	printSpline(kokkosBSP);
	

	//We can now test taking the first derivative
	std::cout << "----Serial 1st derivative----" << std::endl;
	for (int i = 0; i < ctrlPtsX.size(); i++) {
	    std::cout <<"Point: " << ctrlPtsX[i] << ", ";
	    std::cout << ctrlPtsY[i] << std::endl;
	}
	std::cout << "1st derivative, input x = 0: " << serialBSP.x.evalFirstDeriv(0) << ", " << serialBSP.y.evalFirstDeriv(0) << std::endl;
        std::cout << "1st derivative, input x = 1: " << serialBSP.x.evalFirstDeriv(1) << ", " << serialBSP.y.evalFirstDeriv(1) << std::endl;

	std::cout << "----Kokkos 1st derivative----" << std::endl;
	//Get and copy the control points here so we could print them
	 Kokkos::View<double*, MemSpace> dView = kokkosBSP.getCtrlPts();
	auto host_dView = Kokkos::create_mirror_view(dView);
    	Kokkos::deep_copy(host_dView, dView);                                                      
        for (int i = 1; i < host_dView.extent(0); i+=2) {
            std::cout <<"Point: " <<  host_dView(i-1) << ", " << host_dView(i) << std::endl;
        }
        std::cout << "1st derivative, input x = 0: " << kokkosBSP.eval1stDeriv(0, 0) << std::endl;

        std::cout << "1st derivative, input x = 1: " << kokkosBSP.eval1stDeriv(1, 0) << std::endl;
        
	std::cout << "----Serial 2nd derivative----" << std::endl;
        std::cout << "2nd derivative, input x = 0: " << serialBSP.x.evalSecondDeriv(0) << ", " << serialBSP.y.evalSecondDeriv(0) << std::endl;
        std::cout << "2nd derivative, input x = 1: " << serialBSP.x.evalSecondDeriv(1) << ", " << serialBSP.y.evalSecondDeriv(1) << std::endl;

	std::cout << "----Kokkos 2nd derivative----" << std::endl;
        //Get and copy the control points here so we could print them
        dView = kokkosBSP.getCtrlPts();
        host_dView = Kokkos::create_mirror_view(dView);
        Kokkos::deep_copy(host_dView, dView);
        for (int i = 1; i < host_dView.extent(0); i+=2) {
            std::cout <<"Point: " <<  host_dView(i-1) << ", " << host_dView(i) << std::endl;
        }
	std::cout << "2nd derivative, input x = 0: " << kokkosBSP.eval2ndDeriv(0, 0) << std::endl;

        std::cout << "2nd derivative, input x = 1: " << kokkosBSP.eval2ndDeriv(1, 0) << std::endl;

    }
    Kokkos::finalize();
}
