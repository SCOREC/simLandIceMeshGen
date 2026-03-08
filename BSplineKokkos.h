#ifndef BSPLINEKOKKOS_H
#define BSPLINEKOKKOS_H

#include <Kokkos_Core.hpp>
#include <vector>
#include <string>

class BSplineKokkos : public Expression {
public:
	//Genuinely I don't think these should be placed here. Please ask.
	#ifdef KOKKOS_ENABLE_CUDA
	#define MemSpace Kokkos::CudaSpace
	#endif
	#ifdef KOKKOS_ENABLE_HIP
	#define MemSpace Kokkos::Experimental::HIPSpace
	#endif
	#ifndef MemSpace
  	#define MemSpace Kokkos::HostSpace
	#endif

	using ViewVectorType = Kokkos::View<double*, MemSpace>;
	//Constructors
	BSplineKokkos(int order_p, std::vector<double>& ctrlPts, std::vector<double>& knots, std::vector<double>& weight);

	//Accessors
	int getOrder() const {return order;}
	int getNumCtrlPts() const {return ctrlPts.extent(0);}
	int getNumKnots() const {return knots.extent(0);}
	double getCtrlPt(std::size_t i) const {return ctrlPts(i);}
	double getKnot(std::size_t i) const {return knot(i);}

	
	private:
		int order;
		//ADD WHAT MEMORY SPACE
		ViewVectorType ctrlPts;
		ViewVectorType knots;
		ViewVectorType weights;
		ViewVectorType ctrlPts_1stD;
		ViewVectorType ctrlPts_2ndD;

}
