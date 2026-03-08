#include "BSplineKokkos.h"
#include <Kokkos_Core.hpp>
#include <iostream>
#include <vector>

BSplineKokkos::BSplineKokkos(int orderC, std::vector<double>& ctrlPtsC, std::vector<double>& knotsC, std::vector<double>& weightsC) {
	order = orderC;
	//Allocate appropriate view space based on the number of control points, copy the data over to view
	ctrlPts("ctrlPts", ctrlPtsC.size());
	for (int i = 0; i < ctrlPtsC.size(); i++) {
		ctrlPts(i) = ctrlPtsC[i];
	}
	//Same for knots and weights
	knots("knots", knotsC.size());
	for (int i = 0; i < knotsC.size(); i++) {
		knots(i) = knotsC[i];
	}

	weights("weights", weightsC.size());
	for (int i = 0; i < weightsC.size(); i++) {
		weights(i) = weightsC[i];
	}
	
}
