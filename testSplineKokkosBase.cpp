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
	std::vector<double> dummyCtrlPts{2.14, 3.98, 1.56, 9.10, 11.87};
	std::vector<double> dummyKnots{5.99, 8.98, 17.71, 21.1, 2.154, 9.17, 6.32};
	int dOrder = 2;
	std::vector<double> dummyWeights{0.87, 0.125, 0.03};
	BSplineKokkos<ExecutionSpace> ex1(dOrder, dummyCtrlPts, dummyKnots, dummyWeights);
	Kokkos::finalize();
}
