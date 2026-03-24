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
	//In the original implementation, we had the x and y components separated and was storing redundant knots information
	//The implementation now change to interleave the x and y coordinates
	//We store just 1 set of knots rather than 2
	std::cout << "Start BSpline Construction" << std::endl;
	Kokkos::View<int*, MemSpace> orderV("Orders", 1);
	auto host_orderV = Kokkos::create_mirror_view(orderV);
	host_orderV(0) = order_p;
	Kokkos::deep_copy(orderV, host_orderV);
        order = orderV;
	std::cout << "Order inputted" << std::endl;
	//Initialize control points and their offset
	Kokkos::View<double*, MemSpace> ctrlPtsV("CtrlPoints", 2*ctrlPts_x.size());
	//Create the mirror on host so we can access the vector values
	auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);

	for (int i = 0; i < ctrlPts_x.size(); i+=2) {
	    host_ctrlPtsV(i) = ctrlPts_x[i];
	    host_ctrlPtsV(i+1) = ctrlPts_y[i];
	}

	Kokkos::deep_copy(ctrlPtsV, host_ctrlPtsV);
	ctrlPts = ctrlPtsV;

        Kokkos::View<int*, MemSpace> cPOffsetV("CtrlPointsOffset", 1);
	auto host_cPOffsetV = Kokkos::create_mirror_view(cPOffsetV);
        host_cPOffsetV(0) = 0;	//Start of ctrlPts

	Kokkos::deep_copy(cPOffsetV, host_cPOffsetV);
	cPOffset = cPOffsetV;

	//Initialize knots and their offset
	Kokkos::View<double*, MemSpace> knotsV("Knots", knots_x.size()*2);

	Kokkos::View<int*, MemSpace> knotsOffsetV("KnotsOffset", 2);
	auto host_knotsV = Kokkos::create_mirror_view(knotsV);
	auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
        host_knotsOffsetV(0) = 0;

	for (int i = 0; i < knots_x.size(); i+=2) {
	    host_knotsV(i) = knots_x[i];
	    host_knotsV(i+1) = knots_y[i];
	}

	Kokkos::deep_copy(knotsV, host_knotsV);
	Kokkos::deep_copy(knotsOffsetV, host_knotsOffsetV);
	knots = knotsV;
	knotsOffset = knotsOffsetV;

	//Copy the weights over
	Kokkos::View<double*, MemSpace> weightsV("Weights", weights_p.size());
	auto host_weightsV = Kokkos::create_mirror_view(weightsV);
	for (int i = 0; i < weights_p.size(); i++) {
	    host_weightsV(i) = weights_p[i];
	}
	Kokkos::deep_copy(weightsV, host_weightsV);
	weights = weightsV;
	
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
	    cPOffSize += splines[i].getNumCtrlPts();
	    kOffSize += splines[i].getNumKnots();
	    ctrlPSize += splines[i].getCtrlPtsSize();
	    knotsSize += splines[i].getKnotsSize();
	    orderSize += splines[i].getOrder().extent(0);
	    weightSize += splines[i].getWeightSize();
	}
	//Pre-Allocate space for views 
	Kokkos::View<int*, MemSpace> orderV ("Order", orderSize);
	Kokkos::View<double*, MemSpace> ctrlPtsV ("CtrlPts", ctrlPSize);
	Kokkos::View<int*, MemSpace> cPOffsetV ("CtrlPtsOffset", cPOffSize);
	Kokkos::View<double*, MemSpace> knotsV("Knots", knotsSize);
	Kokkos::View<int*, MemSpace> knotsOffsetV ("KnotsOffset", kOffSize);
	Kokkos::View<double*, MemSpace> weightsV("Weights", weightSize);
	auto host_orderV = Kokkos::create_mirror_view(orderV);
	auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
	auto host_cPOffsetV = Kokkos::create_mirror_view(cPOffsetV);
	auto host_knotsV = Kokkos::create_mirror_view(knotsV);
	auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
	auto host_weightsV = Kokkos::create_mirror_view(weightsV);

	int oidx = 0;
	int cPidx = 0;
	int cOidx = 0;
	int kidx = 0;
	int kOidx = 0;
	int widx = 0;

	//Populate the views with data
	int cpLast = 0;
	int kLast = 0;
		
	//NOTE TO SELF: INDEXING MAY BE WRONG
	for (int i = 0; i < splines.size(); i++) {

	    Kokkos::View<int*, MemSpace> intView = splines[i].getOrder();
	    auto host_intView = Kokkos::create_mirror_view(intView);
	    //Populate order
	    for (int j = 0; j < intView.extent(0); j++) {
		host_orderV(oidx+j)= host_intView(j);
	    }
	    oidx += host_intView.extent(0);
	    //Populate ctrlpts offset
	    //Offset the offset by the last offset in the view
	    intView = splines[i].getCPOffset();
	    host_intView = Kokkos::create_mirror_view(intView);
	    for (int j = 0; j < host_intView.extent(0); j++) {
	        host_cPOffsetV(cOidx+j) = host_intView(j)+cpLast;	
	    }
	    cOidx += host_intView.extent(0);
	    cpLast = host_cPOffsetV(cOidx-1)+(host_cPOffsetV(cOidx-1) - host_cPOffsetV(cOidx-2));
	    //Populate knots offset
	    //Offset the offset by the last offset in the view
	    
	    intView = splines[i].getKnotsOffset();
	    host_intView = Kokkos::create_mirror_view(intView);
	    for (int j = 0; j < host_intView.extent(0); j++) {
		host_knotsOffsetV(kOidx+j) = host_intView(j)+kLast;
	    }
	    kOidx += host_intView.extent(0);
	    kLast = host_knotsOffsetV(kOidx-1)+(host_knotsOffsetV(kOidx-1) - host_knotsOffsetV(kOidx-2));

	    //Populate ctrlpts
	    Kokkos::View<double*, MemSpace> doubleView = splines[i].getCtrlPts();
	    auto host_doubleView = Kokkos::create_mirror_view(doubleView);
	    for (int j = 0; j < host_doubleView.extent(0); j++) {
		host_ctrlPtsV(cPidx+j) = host_doubleView(j);
	    }
	    cPidx += host_doubleView.extent(0);
            //Populate knots
	    doubleView = splines[i].getKnots();
	    host_doubleView = Kokkos::create_mirror_view(doubleView);
	    for (int j = 0; j < host_doubleView.extent(0); j++) {
	        host_knotsV(kidx+j) = host_doubleView(j);
	    }
	    kidx += host_doubleView.extent(0);
	    //Populate the weights
	    doubleView = splines[i].getWeights();
	    host_doubleView = Kokkos::create_mirror_view(doubleView);
	    for (int j = 0; j < host_doubleView.extent(0); j++) {
		host_weightsV(widx+j) = host_doubleView(j);
	    }
	    widx += host_doubleView.extent(0);	
	}
	//Copy the modification made from the host data to the device view
	Kokkos::deep_copy(orderV, host_orderV);
	Kokkos::deep_copy(knotsV, host_knotsV);
	Kokkos::deep_copy(ctrlPtsV, host_ctrlPtsV);
	Kokkos::deep_copy(cPOffsetV, host_cPOffsetV);
	Kokkos::deep_copy(weightsV, host_weightsV);
	Kokkos::deep_copy(knotsOffsetV, host_knotsOffsetV);

	//Update the member variables
	order = orderV;
	knots = knotsV;
	ctrlPts = ctrlPtsV;
	weights = weightsV;
	cPOffset = cPOffsetV;
	knotsOffset = knotsOffsetV;
	calculateDerivCoeff();

    }

    double evalFirstDeriv(double x, int splineo, char coor) const {
	//Find the order based on the spline number given
	int lKnot = order(splineo)-1;

	//Find the range that contains this spline
	//Need to make sure it is within the range of this spline
	int leftPt;
	int bound;
	if (coor == 'x') {
	    //We are calculating the x component
	    leftPt = knotsOffset(2*splineo);
	    bound = knotsOffset(2*splineo+1);
	}
	else {
	    leftPt = knotsOffset(2*splineo+1);
	    bound = knotsOffset(2*(splineo+1));
	}

	//Search within this range
	while (knots(leftPt+1) < x && leftPt < bound) {
	    lKnot++;
	    leftPt++;
	}

	int order_t = order(splineo)-1;

	//Copy to view;
	int idx = 0;
	Kokkos::View<double*, MemSpace> pts("pts", order_t);
	for (int i = leftPt; i < leftPt+order_t; i++) {
	    std::cout << i << "|" << ctrlPts_1stD(i) << std::endl;
	    pts(idx) = ctrlPts_1stD(i);
	    idx++;
	}
	idx = 0;
	Kokkos::View<double*, MemSpace> localKnots("local knots", 2*order_t-2);
	for (int i = lKnot-order_t+2; i < lKnot+order_t+1; i++) {
	    localKnots(idx) = knots(i);
	}

	Kokkos::parallel_for("1st derivative loop", order_t, KOKKOS_LAMBDA(int r) {
	    for (int i = order_t-1; i >= r+1; i--) {
	        double aLeft = localKnots(i-1);
	        double aRight = localKnots(i+order_t-(r+1)-1);
	        double alpha;
	        if (aLeft == aRight) {
		    alpha = 0.;
	        }
	        else {
		    alpha = (x-aLeft)/(aRight-aLeft);
	        }
		pts(i) = (1. - alpha) * pts(i-1)+alpha*pts(i);
	    }
	});
	return pts(order_t-1);
    }
    
    double evalSecondDeriv(double x, int splineo, char coor) const {
	if (order == 2) {
	    return 0;
	}
		
	int lKnot = order(splineo)-1;
	int leftPt;
	int bound;
	if (coor == 'x') {
	    leftPt = (2*splineo);
	    bound = knotsOffset(2*splineo+1);
	} else {
	    leftPt = knotsOffset(2*splineo+1);
	    bound = knotsOffset(2*(splineo+1));
	}

	while (knots(lKnot+1) < x && leftPt < bound) {
	    lKnot++;
	    leftPt++;
	}

	//Populate pts and local knot views
	int order_t = order(splineo)-2;
	Kokkos::View<double*, MemSpace> pts("pts", order_t);
	for (int i = leftPt; i < leftPt+order_t; i++) {
	    pts(i-leftPt) = ctrlPts_2ndD(i);
	} 

	Kokkos::View<double*, MemSpace> localKnots("localKnots", 2*order_t-2);
	int idx = 0;
	for (int i = lKnot-order_t+2; i < lKnot+order; i++) {
	    localKnots(idx) = knots(i);
	    idx++;
	}

	Kokkos::parallel_for ("2nd derivative loop", order_t, KOKKOS_LAMBDA(int r){
	    for (int i = order_t-1; i >= r+1; i--) {
		double aLeft = localKnots(i-1);
		double aRight = localKnots(i+order_t-(r+1)-1);
		double alpha;
		if (aLeft == aRight) {
		    alpha = 0;
		}
		else {
		    alpha = (x-aLeft)/(aRight-aLeft);
		}

		pts(i) = (1. - alpha) * pts(i-1)+alpha*pts(i);
			
	    }
	});
	return pts(order_t-1);
	
    }





	//Accessors
    Kokkos::View<int*, MemSpace> getOrder() const {return order;}
    Kokkos::View<double*, MemSpace> getCtrlPts() const {return ctrlPts;}
    Kokkos::View<double*, MemSpace> getKnots() const {return knots;}
    Kokkos::View<double*, MemSpace> getWeights() const {return weights;}
    Kokkos::View<int*, MemSpace> getCPOffset() const {return cPOffset;}
    Kokkos::View<int*, MemSpace> getKnotsOffset() const {return knotsOffset;}
    Kokkos::View<double*, MemSpace> get1stD() const {return ctrlPts_1stD;}
    Kokkos::View<double*, MemSpace> get2ndD() const {return ctrlPts_2ndD;}


    int getNumCtrlPts() const {return cPOffset.extent(0)/2;}
    int getNumKnots() const {return knotsOffset.extent(0)/2;}
    int getCtrlPtsSize() const {return ctrlPts.extent(0);}
    int getKnotsSize() const {return knots.extent(0);}
    int getWeightSize() const {return weights.extent(0);}
    int getOrderSize() const {return order.extent(0);}

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
	
	std::cout << "calculating derivative coeff" << std::endl; 
	Kokkos::View<double*, MemSpace> ctrlPts_1stDV("ctrlPts1Derivative", ctrlPts.extent(0)-cPOffset.extent(0));
	Kokkos::View<int*, MemSpace> cP1stDOffsetV("ctrlPts1DerivativeOffset", cPOffset.extent(0));
	Kokkos::View<int*, MemSpace> cP2ndDOffsetV("ctrlPts2DerivativeOffset", cPOffset.extent(0));
        
	auto host_cP2ndDOffsetV = Kokkos::create_mirror_view(cP2ndDOffsetV);

	auto host_cP1stDOffsetV = Kokkos::create_mirror_view(cP1stDOffsetV);
	host_cP2ndDOffsetV(0) = 0;
	
	std::cout << "offset views and 1stctrlPts views construct" << std::endl;

	//Adjust the offset to include be 1 & 2 less than the ctrlPtsOffset
	host_cP1stDOffsetV(0) = 0;
	for (int i = 1; i < cPOffset.extent(0); i++) {
	    host_cP1stDOffsetV(i) = cPOffset(i)-1;
	    host_cP2ndDOffsetV(i) = cPOffset(i)-2;
	}

	//We need to partition the x and y while we calculate the coefficient
	int idx = 1;
	int oidx = 0;


	for (int i = 1; i < ctrlPts.extent(0); i++) {
	    if (i == cPOffset(idx)) {
		//Do not calculate, delta will be based on both x and y
		idx++;
		if (idx % 2 != 0) {
		    oidx++;
		    //This is a different spline with different order
		}
		continue;
	    }
	    //std::cout << "Order: " << order(oidx) << std::endl;
	    //std::cout << "knots(i): " << knots(i) << std::endl;
	    //std::cout << "i+order(oidx)-1: " << knots(i+order(oidx)-1) << std::endl;
	    double delta = double(order(oidx) - 1)/(knots(i+order(oidx)-1)-knots(i));

	    //std::cout << ctrlPts(i-1) << "|" << ctrlPts(i) << std::endl;


	    //std::cout << "delta: " << delta << std::endl;
	    //std::cout <<"current point: "<<  ctrlPts(i) << std::endl;
	    ctrlPts_1stDV(i-1-(idx-1)) = ((ctrlPts(i) - ctrlPts(i-1)*delta));
	    //std::cout <<"cP coefficient: "<< ctrlPts_1stDV(i-1) << std::endl;
	}

	//Calculate second order derivative
	Kokkos::View<double*, MemSpace> ctrlPts_2ndDV("ctrlPts2Derivative", ctrlPts.extent(0)-2*cPOffset.extent(0));

	idx = 1;
	oidx = 0;

	for (int i = 1; i < ctrlPts_1stDV.extent(0); i++) {
	    if (i == host_cP1stDOffsetV(idx)) {
	        idx++;
		if (idx %2 != 0) {
		    oidx++;
		}
		continue;
	    }
	    if (order(oidx) == 2) {
		//This will cause problem, it should be 0 anyway
		ctrlPts_2ndDV(i-1-(idx-1)) = 0;
		continue;
			
	    }
		double delta = double((order(oidx)-2))/(knots(i+order(oidx)-1)-knots(i+1));
		ctrlPts_2ndDV(i-1-(idx-1)) = ctrlPts_1stDV(i) - ctrlPts_1stDV(i-1)*delta;
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
