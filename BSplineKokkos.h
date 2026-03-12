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
	BSplineKokkos(int order_p, std::vector<double>& ctrlPts_x, std::vector<double>& ctrlPts_y, std::vector<double>& knots_x, std::vector<double>& knots_y, std::vector<double>& weights_p) {
		Kokkos::View<int*, MemSpace> orderV("Orders", 1);
		orderV(0) = order_p;
		order = orderV;
		
		//Initialize control points and their offset
		Kokkos::View<double*, MemSpace> ctrlPtsV("CtrlPoints", 2*ctrlPts_x.size());
		for (int i = 0; i < ctrlPts_x.size(); i++) {
			ctrlPtsV(i) = ctrlPts_x[i];
		}

		for (int i = 0; i < ctrlPts_y.size(); i++) {
			ctrlPtsV(ctrlPts_x.size()+i) = ctrlPts_y[i];
		}

		ctrlPts = ctrlPtsV;

		Kokkos::View<int*, MemSpace> cPOffsetV("CtrlPointsOffset", sizeof(int)*2);
		cPOffsetV(0) = 0;		 //x coor start index
		cPOffsetV(1) = ctrlPts_x.size(); //y coor start index

		cPOffset = cPOffsetV;
/*
		//Initialize knots and their offset
		Kokkos::View<double*, MemSpace> knotsV("Knots", knots_x.size()*2);

		Kokkos::View<int*, MemSpace> knotsOffsetV("KnotsOffset", 2*sizeof(int));
                knotsOffsetV(0) = 0;
		knotsOffsetV(1) = knots_x.size();

		for (int i = 0; i < knots_x.size(); i++) {
			knotsV(i) = knots_x[i];
		}
		for (int i = 0; i < knots_y.size(); i++) {
			knotsV(knotsOffset(1)+i) = knots_y[i];
		}

		knots = knotsV;
		knotsOffset = knotsOffsetV;

		//Copy the weights over
		weights("Weights", weights_p.size());
		for (int i = 0; i < weights_p.size(); i++) {
			weights(i) = weights[i];
		}

		//Derivative coefficients
		//calculateDerivCoeff();
*/
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
		knots("Knots", knotsSize);
		knotsOffset("KnotsOffset", kOffSize);
		weights("Weights", weightSize);
		
		//Populate the views with data

	}



	//Accessors
	int getOrder(int idx) const {return order(idx);}
	int getOrderSize() const {return order.extent(0);}
	int getNumCtrlPts() const {return cPOffset.extent(0)/2;}
	int getNumKnots() const {return knotsOffset.extent(0)/2;}
	int getCtrlPtsSize() const {return ctrlPts.extent(0);}
	int getKnotsSize() const {return knots.extent(0);}
	int getWeightSize() const {return weights.extent(0);}

	double getCtrlPtCoor(int BspIdx, int cPIdx, char coor) const {
		if (coor == 'x') {
			//Find the correct offset
			return ctrlPts(cPOffset(2*BspIdx)+cPIdx);
		}
		else {
			return ctrlPts(cPOffset(2*BspIdx+1)+cPIdx);	
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

	void calculateDerivCoeff() {
		//Moved here from the .cpp file since it uses <ExecutionSpace>
		//Calculate first order derivative
		//Allocate the space
		ctrlPts_1stD("ctrlPts1Derivative", ctrlPts.extent(0)-cPOffset.extent(0));
		cP1stDOffset("ctrlPts1DerivativeOffset", cPOffset.extent(0));
		ctrlPts_2ndD("ctrlPts1Derivative", ctrlPts.extent(0)-2*(cPOffset.extent(0)));
                cP2ndDOffset("ctrlPts2DerivativeOffset", cPOffset.extent(0));
                cP2ndDOffset(0) = 0;

		//Adjust the offset to include be 1 & 2 less than the ctrlPtsOffset
		cP1stDOffset(0) = 0;
		for (int i = 1; i < cPOffset.extent(0); i++) {
			cP1stDOffset(i) = cPOffset(i)-1;
			cP2ndDOffset(i) = cPOffset(i)-2;
		}

		//We need to partition the x and y while we calculate the coefficient
		int idx = 1;
		int oidx = 0;
		for (int i = 1; i < cPOffset(cPOffset.extend(0)-2); i++) {
			if (i == cPOffset(idx)) {
				//Do not calculate, delta will be based on both x and y
				idx++;
				if (idx % 2 == 0) {
					oidx++;
				}
				continue;
			}
			double delta = double(order(oidx) - 1)/(knots(i+order(oidx)-1)-knots(i));
			ctrlPts_1stD(i-1) = ((ctrlPts(i) - ctrlPts(i-1)*delta));
		}

		//Calculate second order derivative
		ctrlPts_2ndD("ctrlPts1Derivative", ctrlPts.extent(0)-2*(cPOffset.extent(0)));
		cP2ndDOffset("ctrlPts2DerivativeOffset", cPOffset.extent(0));
		cP2ndDOffset(0) = 0;

		idx = 1;
		oidx = 0;

		for (int i = 1; i < ctrlPts_1stD.extent(0); i++) {
			if (i == cP1stDOffset(idx)) {
				idx++;
				if (idx %2 == 0) {
					oidx++;
				}
				continue;
			}
			double delta = double((order(oidx)-2))/(knots(i+order(oidx)-1)-knots(i+1));
			ctrlPts_2ndD(i-1) = ctrlPts_1stD(i) - ctrlPts_1stD(i-1)*delta;
		}
	}
	
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
