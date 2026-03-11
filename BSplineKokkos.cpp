#include "BSplineKokkos.h"
#include <iostream>
template<typename ExecutionSpace>
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
	
	/*NOTE: UNCOMMENT THIS ONCE ALL IN FRONT ARE RESOLVED*/
	//calculateDerivCoeff();
	
}
template<typename ExecutionSpace>
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

}
template<typename ExecutionSpace>
double BSplineKokkos<ExecutionSpace>::eval(double x) const {
	//Implemented based on the serial BSpline
	//Find the interval of the 1D coordinate given
	int leftKnot = order - 1;
	int leftPt = 0;
	while (knots(leftKnot) < x) {
		leftKnot++;
		leftPt++;
		if (leftKnot == knots.extent(0)-1) {
		break;
		} 
	}
	//Find the points and local knots
	Kokkos::View<double*, MemSpace> pts("PtsView", order);
	for (int i = leftKnot; i < leftKnot+order; i++) {
		pts(i) = ctrlPts(i);
	}
	//TODO: find the size of the local knot and copy them to a view
	Kokkos::View<double*, MemSpace> localKnots;
	for (int r = 1; r <= order; r++) {
		for (int i = order-1; i>= r; i--) {
			double a_left = localKnots(i-1);
			double a_right = localKnots(i-1);
			double alpha;
			if (a_right == a_left) {
				alpha = 0.;
			}
			else {
				alpha = (x-a_left)/(a_right-a_left);
			}
			pts(i) = (1. - alpha) * pts(i-1) + alpha*pts(i);
		}
	}
	return pts(order-1);

}

//Explicit instantiation of the templated class for Kokkos::serial
//template class BSplineKokkos<Kokkos::Serial>;
