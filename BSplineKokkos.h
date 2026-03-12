#ifndef BSPLINEKOKKOS_H
#define BSPLINEKOKKOS_H


#include <Kokkos_Core.hpp>
#include <vector>
#include <string>
#include <iostream>

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
		std::cout << "Order view added" << std::endl; 
		//Initialize control points and their offset
		Kokkos::View<double*, MemSpace> ctrlPtsV("CtrlPoints", 2*ctrlPts_x.size());
		for (int i = 0; i < ctrlPts_x.size(); i++) {
			ctrlPtsV(i) = ctrlPts_x[i];
		}

		for (int i = 0; i < ctrlPts_y.size(); i++) {
			ctrlPtsV(ctrlPts_x.size()+i) = ctrlPts_y[i];
		}
		std::cout << "ctrlPts initialized" << std::endl; 

		ctrlPts = ctrlPtsV;

		Kokkos::View<int*, MemSpace> cPOffsetV("CtrlPointsOffset", 2);
		cPOffsetV(0) = 0;		 //x coor start index
		cPOffsetV(1) = ctrlPts_x.size(); //y coor start index

		cPOffset = cPOffsetV;

		std::cout << "cp offset initialized" << std::endl;

		//Initialize knots and their offset
		Kokkos::View<double*, MemSpace> knotsV("Knots", knots_x.size()*2);

		Kokkos::View<int*, MemSpace> knotsOffsetV("KnotsOffset", 2);
                knotsOffsetV(0) = 0;
		knotsOffsetV(1) = knots_x.size();

		for (int i = 0; i < knots_x.size(); i++) {
			knotsV(i) = knots_x[i];
		}
		for (int i = 0; i < knots_y.size(); i++) {
			knotsV(knotsOffsetV(1)+i) = knots_y[i];
		}

		knots = knotsV;
		knotsOffset = knotsOffsetV;

		std::cout << "knots initialized" << std::endl;
		std::cout << "knots offset initialized" << std::endl;

		//Copy the weights over
		Kokkos::View<double*, MemSpace> weightsV("Weights", weights_p.size());
		for (int i = 0; i < weights_p.size(); i++) {
			weightsV(i) = weights_p[i];
		}
		weights = weightsV;

		std::cout << "weights initialized" << std::endl;

		//Derivative coefficients
		calculateDerivCoeff();

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
		Kokkos::View<int*, MemSpace> orderV ("Order", orderSize);
		Kokkos::View<double*, MemSpace> ctrlPtsV ("CtrlPts", ctrlPSize);
		Kokkos::View<int*, MemSpace> cPOffsetV ("CtrlPtsOffset", cPOffSize);
		Kokkos::View<double*, MemSpace> knotsV("Knots", knotsSize);
		Kokkos::View<int*, MemSpace> knotsOffsetV ("KnotsOffset", kOffSize);
		Kokkos::View<double*, MemSpace> weightsV("Weights", weightSize);
		
		int oidx = 0;
		int cPidx = 0;
		int cOidx = 0;
		int kidx = 0;
		int kOidx = 0;
		int widx = 0;
		//Populate the views with data
		int lastOffset = 0;
		
		for (int i = 0; i < splines.size(); i++) {
			Kokkos::View<int*, MemSpace> intView = splines[i].getOrder();
			//Populate order
			for (int j = 0; j < intView.extent(0); j++) {
				orderV(oidx+j)= intView(j);
			}
			oidx += intView.extent(0);
			//Populate ctrlpts offset
			//Offset the offset by the last offset in the view
			intView = splines[i].getCPOffset();
			for (int j = 0; j < intView.extent(0); j++) {
				cPOffsetV(cOidx+j) = intView(j)+lastOffset; 
			}
			cOidx += intView.extent(0);
			lastOffset = intView(intView.extent(0)-1)+(intView(intView.extent(0)-1) - intView(intView.extent(0)-2));
			//Populate knots offset
			//Offset the offset by the last offset in the view
			intView = splines[i].getKnotsOffset();
			for (int j = 0; j < intView.extent(0); j++) {
				knotsOffsetV(kOidx+j) = intView(j)+lastOffset;
			}
			kOidx += intView.extent(0);
			lastOffset = intView(intView.extent(0)-1)+(intView(intView.extent(0)-1) - intView(intView.extent(0)-2));
			//Populate ctrlpts
			Kokkos::View<double*, MemSpace> doubleView = splines[i].getCtrlPts();
			for (int j = 0; j < doubleView.extent(0); j++) {
				ctrlPtsV(cPidx+j) = doubleView(j);
			}
			cPidx += doubleView.extent(0);
			//Populate knots
			doubleView = splines[i].getKnots();
			for (int j = 0; j < doubleView.extent(0); j++) {
				knotsV(kidx+j) = doubleView(j);
			}
			kidx += doubleView.extent(0);
			//Populate the weights
			doubleView = splines[i].getWeights();
			for (int j = 0; j < doubleView.extent(0); j++) {
				weightsV(widx+j) = doubleView(j);
			}
			widx += doubleView.extent(0);	
		}
		order = orderV;
		knots = knotsV;
		ctrlPts = ctrlPtsV;
		weights = weightsV;
		cPOffset = cPOffsetV;
		knotsOffset = knotsOffsetV;

	}



	//Accessors
	Kokkos::View<int*, MemSpace> getOrder() const {return order;}
	Kokkos::View<double*, MemSpace> getCtrlPts() const {return ctrlPts;}
	Kokkos::View<double*, MemSpace> getKnots() const {return knots;}
	Kokkos::View<double*, MemSpace> getWeights() const {return weights;}
	Kokkos::View<int*, MemSpace> getCPOffset() const {return cPOffset;}
	Kokkos::View<int*, MemSpace> getKnotsOffset() const {return knotsOffset;}

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
		Kokkos::View<double*, MemSpace> ctrlPts_1stDV("ctrlPts1Derivative", ctrlPts.extent(0)-cPOffset.extent(0));
		Kokkos::View<int*, MemSpace> cP1stDOffsetV("ctrlPts1DerivativeOffset", cPOffset.extent(0));
		Kokkos::View<int*, MemSpace> cP2ndDOffsetV("ctrlPts2DerivativeOffset", cPOffset.extent(0));
                cP2ndDOffsetV(0) = 0;

		//Adjust the offset to include be 1 & 2 less than the ctrlPtsOffset
		cP1stDOffsetV(0) = 0;
		for (int i = 1; i < cPOffset.extent(0); i++) {
			cP1stDOffsetV(i) = cPOffset(i)-1;
			cP2ndDOffsetV(i) = cPOffset(i)-2;
		}

		//We need to partition the x and y while we calculate the coefficient
		int idx = 1;
		int oidx = 0;
		for (int i = 1; i < cPOffset(cPOffset.extent(0)-2); i++) {
			if (i == cPOffset(idx)) {
				//Do not calculate, delta will be based on both x and y
				idx++;
				if (idx % 2 == 0) {
					oidx++;
				}
				continue;
			}
			double delta = double(order(oidx) - 1)/(knots(i+order(oidx)-1)-knots(i));
			ctrlPts_1stDV(i-1) = ((ctrlPts(i) - ctrlPts(i-1)*delta));
		}

		//Calculate second order derivative
		Kokkos::View<double*, MemSpace> ctrlPts_2ndDV("ctrlPts2Derivative", ctrlPts.extent(0)-2*(cPOffset.extent(0)));

		idx = 1;
		oidx = 0;

		for (int i = 1; i < ctrlPts_1stDV.extent(0); i++) {
			if (i == cP1stDOffsetV(idx)) {
				idx++;
				if (idx %2 == 0) {
					oidx++;
				}
				continue;
			}
			double delta = double((order(oidx)-2))/(knots(i+order(oidx)-1)-knots(i+1));
			ctrlPts_2ndDV(i-1) = ctrlPts_1stDV(i) - ctrlPts_1stDV(i-1)*delta;
		}
		ctrlPts_2ndD = ctrlPts_2ndDV;
		ctrlPts_1stD = ctrlPts_1stDV;
		cP1stDOffset = cP1stDOffsetV;
		cP2ndDOffset = cP2ndDOffsetV;
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
