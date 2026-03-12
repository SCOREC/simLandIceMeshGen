#include "BSplineKokkos.h"
#include <iostream>
/*template<typename ExecutionSpace>
BSplineKokkos<ExecutionSpace>::BSplineKokkos(int orderC, std::vector<double>& ctrlPtsC, std::vector<double>& knotsC, std::vector<double>& weightsC) {
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
	//Call the calculateDerivCoeff() to populate
	//1st and 2nd deriavtive views
	
	//NOTE: UNCOMMENT THIS ONCE ALL IN FRONT ARE RESOLVED
	//calculateDerivCoeff();
	
}*/


/*template<typename ExecutionSpace>
void BSplineKokkos<ExecutionSpace>::calculateDerivCoeff() {
	//Calculate first order derivative
	//Allocate space for ctrlPts_1stD
	ctrlPts_1stD("ctrlPts1Derivative", ctrlPts.extent(0)-1);	
	for (int i = 1; i < ctrlPts.extent(0); i++) {	
		double delta = double(order - 1)/(knots(i+order-1)-knots(i));
		ctrlPts_1stD(i) = ((ctrlPts(i) - ctrlPts(i-1)*delta));

	}

	//Calculate second order derivative
	//Allocate space for ctrlPts_2ndD
	ctrlPts_2ndD("ctrlPts2Derivative", ctrlPts.extent(0)-2);
	for (int i = 0; i < ctrlPts_1stD.extent(0); i++) {
		double delta = double(order - 2) / (knots(i+order-1)-knots(i+1));
		ctrlPts_2ndD(i) = ((ctrlPts_1stD(i) - ctrlPts_1stD(i-1))*delta);
	}
	//TODO: find another way to verify the size of the second derivative view

}*/

//Explicit instantiation of the templated class for Kokkos::serial
//template class BSplineKokkos<Kokkos::Serial>;
