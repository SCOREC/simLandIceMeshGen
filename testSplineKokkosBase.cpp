#include <Kokkos_Core.hpp>
#include "BSplineKokkos.h"
#include <vector>
#include <iostream>

template<typename MemSpace>
void printDView(Kokkos::View<double*, MemSpace>& dView) {
	for (int i = 0; i < dView.extent(0); i++) {
		std::cout << dView(i) << std::endl;
	}
}

template<typename MemSpace>
void printIView(Kokkos::View<int*, MemSpace>& iView) {
	for (int i = 0; i < iView.extent(0); i++) {
		std::cout << iView(i) << std::endl;
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
		std::vector<double> dummyKnotsY{2.85, 1.16, 0.82, 1.89, 9.82, 12.0, 7.19};
		int dOrder = 2;
		std::vector<double> dummyWeights{0.87, 0.125, 0.03};
		//Initialize 1 Kokkos Spline
		BSplineKokkos<ExecutionSpace> ex1(dOrder, dummyCtrlPtsX, dummyCtrlPtsY, dummyKnotsX, dummyKnotsY, dummyWeights);

		//Check if we actually initialized with correct values
		std::cout << "ex1" << std::endl;
		printSpline(ex1);
		std::vector<BSplineKokkos<ExecutionSpace>> splineList;
		splineList.push_back(ex1);

		//Test constructor that takes in multiple BSplineKokkos objects
		dummyCtrlPtsX.assign({1.20, 12.91, 7.11});
		dummyCtrlPtsY.assign({3.31, 7.89, 9.12});
		dummyKnotsX.assign({2.18, 5.56, 3.18, 16.15});
		dummyKnotsY.assign({1.971,4.12, 8.192, 12.11});
		dOrder = 3;
		dummyWeights.assign({0.19, 0.27});
		BSplineKokkos<ExecutionSpace> ex2(dOrder, dummyCtrlPtsX, dummyCtrlPtsY, dummyKnotsX, dummyKnotsY, dummyWeights);
		splineList.push_back(ex2);
		std::cout << "ex2" << std::endl;
		printSpline(ex2);

		dummyCtrlPtsX.assign({18.29, 12.12, 5.92, 7.182, 1.16});
		dummyCtrlPtsY.assign({4.72, 5.78, 9.12, 1.29, 10.2});
		dummyKnotsX.assign({4.34, 29.3, 18.7});
		dummyKnotsY.assign({19.12, 10.21, 12.78});
		dOrder = 5;
                dummyWeights.assign({0.12, 0.182, 0.98, 0.05});
		std::cout << "\t\tSTART OF MULTISPLINE CREATION" << std::endl;
		BSplineKokkos<ExecutionSpace> ex3(dOrder, dummyCtrlPtsX, dummyCtrlPtsY, dummyKnotsX, dummyKnotsY, dummyWeights);
		splineList.push_back(ex3);
		std::cout << "ex3" << std::endl;
		printSpline(ex3);

		BSplineKokkos<ExecutionSpace> fromList(splineList);

		std::cout << "Spline created from list of splines" << std::endl;
		printSpline(fromList);
		 		 
	}
	Kokkos::finalize();
}
