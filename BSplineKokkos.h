#ifndef BSPLINEKOKKOS_H
#define BSPLINEKOKKOS_H


#include <Kokkos_Core.hpp>
#include <vector>
#include <string>

//For the Kokkos version of the implementation, we will not be separating BSplineKokkos and BSplineKokkos2D, their respective functions will be in this file.

template<typename ExecutionSpace>
class BSplineKokkos {
public:
	using MemSpace = typename ExecutionSpace::memory_space;
	
	//Constructor for just 1 spline(Would this be used often?)
	BSplineKokkos(int order_p, std::vector<double>& ctrlPts_x, std::vector<double>& ctrlPts_y, std::vector<double>& knots_x, std::vector<double>& knots_y, std::vector<double>& weight_p) {
		order("Orders", sizeof(int));
		order(0) = order_p;
		
		//Initialize control points and their offset
		ctrlPts("CtrlPoints", 2*ctrlPts_x.size());
		for (int i = 0; i < ctrlPts_x.size(); i++) {
			ctrlPts(i) = ctrlPts_x[i];
		}

		for (int i = 0; i < ctrlPts_y.size(); i++) {
			ctrlPts(ctrlPts_x.size()+i) = ctrlPts_y[i];
		}

		cPOffset("CtrlPointsOffset", sizeof(int)*2);
		cPOffset(0) = 0;		//x coor start index
		cPOffset(1) = ctrlPts_x.size();	//y coor start index

		//Initialize knots and their offset
		knots("Knots", knots_x.size()*2);

		knotsOffset("KnotsOffset", 2*sizeof(int));
                knotsOffset(0) = 0;
		knotsOffset(1) = knots_x.size();

		for (int i = 0; i < knots_x.size(); i++) {
			knots(i) = knots_x[i];
		}
		for (int i = 0; i < knots_y.size(); i++) {
			knots(knotsOffset(1)+i) = knots_y[i];
		}

		//Copy the weights over
		weights("Weights", weights_p.size());
		for (int i = 0; i < weights_p.size(); i++) {
			weights(i) = weights[i];
		}

		//Derivative coefficients
		//calculateDerivCoeff();

	}

	//Create BSpline for a vector of BSplineKokkos objects
	BSplineKokkos(std::vector<BSplineKokkos>& splines) {
		//Preallocate space for Kokkos::View
		//They are fixed in size and resize is expensive
		int orderSize = 0;
		int ctrlPSize = 0;
		int knotsSize = 0;
		int cPOffSize = 0;
		int kOffSize = 0;
		int weightSize = 0;

		//Input the offsets
		//Track space needed to allocate for ctrlPts & knots
		for (int i = 0; i < splines.size(); i++) {
			cPOffSize += 2*(splines[i].getNumCtrlPts());
			knotsSize += 2*(splines[i].getNumKnots());
			ctrlPSize += splines[i].getCtrlPtsSize();
			knotsSize += splines[i].getKnotsSize();
			orderSize += splines[i].getOrderSize();
			weightSize += splines[i].getWeightSize();
		}
		//Pre-Allocate space for views 
		order("Order", orderSize);
		ctrlPts("CtrlPts", ctrlPSize);
		cPOffset("CtrlPtsOffset", cPOffSize);
		knotsSize("Knots", knotsSize);
		kOffSize("KnotsOffset", kOffSize);
		weights("Weights", weightSize);
		
		//Populate the views with data

	}



	//Accessors
	int getOrder(int idx) const {return order(idx);}
	int getOrderSize() const {return order.extent(0);}
	int getNumCtrlPts() const {return ctrlPtsOffset.extent(0)/2;}
	int getNumKnots() const {return knotsOffset.extent(0)/2;}
	int getCtrlPtsSize() const {return ctrlPts.extent(0);}
	int getKnotsSize() const {return knots.extent(0);}
	int getWeightSize() const {return weights.extent(0);}

	double getCtrlPtCoor(int BspIdx, int cPIdx, char coor) const {
		if (coor == 'x') {
			//Find the correct offset
			return ctrlPts(ctrlPtsOffset(2*BspIdx)+cPIdx);
		}
		else {
			return ctrlPts(ctrlPtsOffset(2*BspIdx+1)+cPIdx);	
		}
	}
	double getKnotCoor(int BspIdx, int kIdx, char coor) const {
		if (coor == 'x') {
			return knots(knotsOffset(2*BspIdx)+kIdx);
		}
		else {
			return knots(knotsOffset(2*BspIdx+1)+kIdx);
		}
	}

	void calculateDerivCoeff();
	
	private:
		Kokkos::View<int*, MemSpace> order;

		//CtrlPts, knots, weights and their offsets
		Kokkos::View<double*, MemSpace> ctrlPts;
		Kokkos::View<int*, MemSpace> cPOffset;
		Kokkos::View<double*, MemSpace> knots;
		Kokkos::View<int*, MemSpace> knotsOffset;
		Kokkos::View<double*, MemSpace> weights;
		//weight offset may be added later, not entirely sure for now

		//The 1st and 2nd derivatives
		Kokkos::View<double*, MemSpace> ctrlPts_1stD;
		Kokkos::View<int*, MemSpace> cP1stDOffset;
		Kokkos::View<double*, MemSpace> ctrlPts_2ndD;
		Kokkos::View<int*, MemSpace> cP2ndDOffset;

};
#endif
