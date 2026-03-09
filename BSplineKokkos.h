#ifndef BSPLINEKOKKOS_H
#define BSPLINEKOKKOS_H

#include <Kokkos_Core.hpp>
#include <vector>
#include <string>

template<typename ExecutionSpace>
class BSplineKokkos : public Expression {
public:
	using MemSpace = typename ExecutionSpace::memory_space;
	
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

		//CtrlPts, knots, weights and their offsets
		Kokkos::View<double*, MemSpace> ctrlPts;
		Kokkos::View<double*, MemSpace> knots;
		Kokkos::View<double*, MemSpace> weights;

		//The 1st and 2nd derivatives
		Kokkos::View<double*, MemSpace> ctrlPts_1stD;
		Kokkos::View<double*, MemSpace> ctrlPts_2ndD;

}
