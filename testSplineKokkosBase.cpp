#include <Kokkos_Core.hpp>
#include "BSplineKokkos.h"
#include <vector>

int main(int argc, char* argv[]) {
	Kokkos::initialize(argc, argv);

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
	BSplineKokkos<ExecutionSpace> ex1(dOrder, dummyCtrlPtsX, dummyCtrlPtsY, dummyKnotsX, dummyKnotsY, dummyWeights);
	Kokkos::finalize();
}
