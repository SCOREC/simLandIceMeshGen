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
    BSplineKokkos(int order_p, std::vector<double>& ctrlPts_x, std::vector<double>& ctrlPts_y, std::vector<double>& knotsI) {
	//In the original implementation, we had the x and y components separated and was storing redundant knots information
	//The implementation now change to interleave the x and y coordinates
	//We store just 1 set of knots rather than 2
	std::cout << "Start BSpline Construction" << std::endl;
	Kokkos::View<int*, MemSpace> orderV("Orders", 1);
	auto host_orderV = Kokkos::create_mirror_view(orderV);
	host_orderV(0) = order_p;
	order = orderV;
	Kokkos::deep_copy(order, host_orderV);
	std::cout << "Order inputted" << std::endl;
	//Initialize control points and their offset
	Kokkos::View<double*, MemSpace> ctrlPtsV("CtrlPoints", 2*ctrlPts_x.size());
	//Create the mirror on host so we can access the vector values
	auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
	std::cout << "ctrlPts size: " << host_ctrlPtsV.extent(0) << std::endl;

	for (int i = 0; i < ctrlPts_x.size(); i++) {
	    host_ctrlPtsV(2*i) = ctrlPts_x[i];
	    host_ctrlPtsV((2*i)+1) = ctrlPts_y[i];
	}

	ctrlPts = ctrlPtsV;
	Kokkos::deep_copy(ctrlPts, host_ctrlPtsV);

        Kokkos::View<int*, MemSpace> cPOffsetV("CtrlPointsOffset", 1);
	auto host_cPOffsetV = Kokkos::create_mirror_view(cPOffsetV);
        host_cPOffsetV(0) = 2*ctrlPts_x.size();	//End of ctrlPts

	cPOffset = cPOffsetV;
	Kokkos::deep_copy(cPOffset, host_cPOffsetV);

	//Initialize knots and their offset
	Kokkos::View<double*, MemSpace> knotsV("Knots", knotsI.size());

	Kokkos::View<int*, MemSpace> knotsOffsetV("KnotsOffset", 1);
	auto host_knotsV = Kokkos::create_mirror_view(knotsV);
	auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
        host_knotsOffsetV(0) = knotsI.size();   //End of knots

	for (int i = 0; i < knotsI.size(); i++) {
	    host_knotsV(i) = knotsI[i];
	}

	knots = knotsV;
	Kokkos::deep_copy(knots, host_knotsV);

	knotsOffset = knotsOffsetV;
	Kokkos::deep_copy(knotsOffset, host_knotsOffsetV);

	std::cout << "All good, before derivative coe" << std::endl;	
	//Derivative coefficients
	calculateDerivCoeff();
	std::cout << "single spline initialization ok" << std::endl;
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

	//Input the offsets
	//Track space needed to allocate for ctrlPts & knots
	for (int i = 0; i < splines.size(); i++) {
	    cPOffSize += splines[i].getCPOffsetSize();
	    kOffSize += splines[i].getKnotsOffsetSize();
	    ctrlPSize += splines[i].getCtrlPtsSize();
	    knotsSize += splines[i].getKnotsSize();
	    orderSize += splines[i].getOrder().extent(0);
	}
	//Pre-Allocate space for views 
	Kokkos::View<int*, MemSpace> orderV ("Order", orderSize);
	Kokkos::View<double*, MemSpace> ctrlPtsV ("CtrlPts", ctrlPSize);
	Kokkos::View<int*, MemSpace> cPOffsetV ("CtrlPtsOffset", cPOffSize);
	Kokkos::View<double*, MemSpace> knotsV("Knots", knotsSize);
	Kokkos::View<int*, MemSpace> knotsOffsetV ("KnotsOffset", kOffSize);
	auto host_orderV = Kokkos::create_mirror_view(orderV);
	auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
	auto host_cPOffsetV = Kokkos::create_mirror_view(cPOffsetV);
	auto host_knotsV = Kokkos::create_mirror_view(knotsV);
	auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);


	std::cout << "ctrlPtsOffset size: " << cPOffsetV.extent(0)<< std::endl;
	std::cout << "knotsOffset size: " << knotsOffsetV.extent(0) << std::endl;
	int oidx = 0;
	int cPidx = 0;
	int cOidx = 0;
	int kidx = 0;
	int kOidx = 0;

	//Populate the views with data
	int cpLast = 0;
	int kLast = 0;
		
	//NOTE TO SELF: INDEXING MAY BE WRONG
	for (int i = 0; i < splines.size(); i++) {
	    Kokkos::View<int*, MemSpace> intView = splines[i].getOrder();
	    auto host_intView = Kokkos::create_mirror_view(intView);
	    Kokkos::deep_copy(host_intView, intView);
	    //Populate order
	    for (int j = 0; j < intView.extent(0); j++) {
		host_orderV(oidx+j)= host_intView(j);
		std::cout << "Spline Order: " << host_orderV(oidx+j) << std::endl;
	    }
	    oidx += host_intView.extent(0);
	    //Populate ctrlpts offset
	    //Offset the offset by the last offset in the view
	    intView = splines[i].getCPOffset();
	    host_intView = Kokkos::create_mirror_view(intView);
	    Kokkos::deep_copy(host_intView, intView);
	    for (int j = 0; j < host_intView.extent(0); j++) {
		std::cout << "\tctrlPtsOffset: " << host_intView(j) << std::endl;
		std::cout << "\tcpLast: " << cpLast << std::endl;
		std::cout << "cOidx: " << cOidx << std::endl;
		host_cPOffsetV(cOidx+j) = host_intView(j)+cpLast;	

	    }
	    cOidx += host_intView.extent(0);
	    cpLast += host_intView(host_intView.extent(0)-1);
	    //Populate knots offset
	    //Offset the offset by the last offset in the view
	    
	    intView = splines[i].getKnotsOffset();
	    host_intView = Kokkos::create_mirror_view(intView);
	    Kokkos::deep_copy(host_intView, intView);
	    for (int j = 0; j < host_intView.extent(0); j++) {
		host_knotsOffsetV(kOidx+j) = host_intView(j)+kLast;
	    }
	    kOidx += host_intView.extent(0);
	    kLast += host_intView(host_intView.extent(0)-1);
	    //kLast = host_knotsOffsetV(kOidx-1)+(host_knotsOffsetV(kOidx-1) - host_knotsOffsetV(kOidx-2));

	    //Populate ctrlpts
	    Kokkos::View<double*, MemSpace> doubleView = splines[i].getCtrlPts();
	    auto host_doubleView = Kokkos::create_mirror_view(doubleView);
	    Kokkos::deep_copy(host_doubleView, doubleView);
	    for (int j = 0; j < host_doubleView.extent(0); j++) {
		host_ctrlPtsV(cPidx+j) = host_doubleView(j);
	    }
	    cPidx += host_doubleView.extent(0);
            //Populate knots
	    doubleView = splines[i].getKnots();
	    host_doubleView = Kokkos::create_mirror_view(doubleView);
	    Kokkos::deep_copy(host_doubleView, doubleView);
	    for (int j = 0; j < host_doubleView.extent(0); j++) {
	        host_knotsV(kidx+j) = host_doubleView(j);
	    }
	    kidx += host_doubleView.extent(0);	
	}

	//Update the member variables
	order = orderV;
	knots = knotsV;
	ctrlPts = ctrlPtsV;
	cPOffset = cPOffsetV;
	knotsOffset = knotsOffsetV;
	//Copy the modification made from the host data to the device view
        Kokkos::deep_copy(order, host_orderV);
        Kokkos::deep_copy(knots, host_knotsV);
        Kokkos::deep_copy(ctrlPts, host_ctrlPtsV);
        Kokkos::deep_copy(cPOffset, host_cPOffsetV);
        Kokkos::deep_copy(knotsOffset, host_knotsOffsetV);
	
	//TO DO: Work on calculate derivative coeff
	calculateDerivCoeff();

    }

    double eval1stDeriv(double x, int splineo) const {
	//Find the order based on the spline number given
	auto order_mv = Kokkos::create_mirror_view(order);
	Kokkos::deep_copy(order_mv, order);
	int lKnot = order_mv(splineo)-1;

	//std::cout << lKnot << std::endl;
	//Find the range that contains this spline
	//Need to make sure it is within the range of this spline

	int leftPtX = 0;    //Even indicies are x coordinates
	int leftPtY = 1;    //Odd indicies are y coordinates

	//Find the boundary to the current spline control points

        auto mv_knots = Kokkos::create_mirror_view(knots);
        Kokkos::deep_copy(mv_knots, knots);

	//Search within this range
	while (mv_knots(lKnot+1) < x) {
	    lKnot+=2;
	    leftPtX+=2;    //Increment xpts
	    leftPtY+=2;    //Increment ypts
	}

	int order_t = order_mv(splineo)-1;

	//Copy our previously calculated coefficient
	auto mv_ctrlPts_1stD = Kokkos::create_mirror_view(ctrlPts_1stD);
	Kokkos::deep_copy(mv_ctrlPts_1stD, ctrlPts_1stD);
	
	int idx = 0;
	Kokkos::View<double*, MemSpace> ptsX("ptsX", order_t);
	auto mv_ptsX = Kokkos::create_mirror_view(ptsX);
	Kokkos::View<double*, MemSpace> ptsY("ptsY", order_t);
	auto mv_ptsY = Kokkos::create_mirror_view(ptsY);
	for (int i = leftPtX; i < leftPtX+order_t-1; i+=2) {
	    //std::cout << i << "|" << mv_ctrlPts_1stD(i) << std::endl;
	    mv_ptsX(idx) = mv_ctrlPts_1stD(i);
	    mv_ptsY(idx) = mv_ctrlPts_1stD(i+1);
	    idx++;
	}

	for (int i = 0; i < mv_ptsX.extent(0); i++) {
	    std::cout << mv_ptsX(i) << std::endl;
	}
	std::cout << "|" << std::endl;

	//We only need 1 copy of the local knots
	idx = 0;
	Kokkos::View<double*, MemSpace> localKnots("local knots", 2*order_t-2);
	auto mv_localKnots = Kokkos::create_mirror_view(localKnots);
	for (int i = lKnot-order_t+2; i < lKnot+order_t+1; i++) {
	    mv_localKnots(idx) = mv_knots(i);
	}
	//std::cout << "Local knots collection complete" << std::endl;

	//Copy allocated value back to device
        Kokkos::deep_copy(localKnots, mv_localKnots);
	Kokkos::deep_copy(ptsX, mv_ptsX);	

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
		ptsX(i) = (1. - alpha) * ptsX(i-1)+alpha*ptsX(i);
	    }
	});
        //std::cout << "Before deep copy of pts" << std::endl;
	Kokkos::deep_copy(mv_ptsX, ptsX);
	//std::cout << "After deep copy of pts" << std::endl;
	for (int i = 0; i < mv_ptsX.extent(0); i++) {
	    std::cout << mv_ptsX(i) << std::endl;
	}
	return mv_ptsX(order_t-1);
    }
    
    double eval2ndDeriv(double x, int splineo) const {
	auto mv_order = Kokkos::create_mirror_view(order);
	Kokkos::deep_copy(mv_order);
	if (mv_order(splineo) == 2) {
	    return 0;
	}
		
	int lKnot = mv_order(splineo)-1;
	int leftPt;
	int bound;

	auto mv_knotsOffset = Kokkos::create_mirror_view(knotsOffset);
	Kokkos::deep_copy(mv_knotsOffset, knotsOffset);
	if (splineo == 0) {
	    leftPt = 0;
	    bound = mv_knotsOffset(splineo);
	} else {
	    leftPt = mv_knotsOffset(splineo-1);
	    bound = mv_knotsOffset(splineo);
	}

	auto mv_knots = Kokkos::create_mirror_view(knots);
	Kokkos::deep_copy(mv_knots, knots);
	while (mv_knots(lKnot+1) < x && leftPt < bound) {
	    lKnot+=2;
	    leftPt+=2;
	}

	//Populate pts and local knot views
	int order_t = mv_order(splineo)-2;
	auto mv_ctrlPts_2ndD = Kokkos::create_mirror_view(ctrlPts_2ndD);
	Kokkos::deep_copy(mv_ctrlPts_2ndD, ctrlPts_2ndD);
	Kokkos::View<double*, MemSpace> pts("pts", order_t);
	auto mv_pts = Kokkos::create_mirror_view(pts);
	for (int i = leftPt; i < leftPt+order_t; i++) {
	    mv_pts(i-leftPt) = mv_ctrlPts_2ndD(i);
	} 

	Kokkos::View<double*, MemSpace> localKnots("localKnots", 2*order_t-2);
	int idx = 0;
	auto mv_localKnots = Kokkos::create_mirror_view(localKnots);
	for (int i = lKnot-order_t+2; i < lKnot+order; i++) {
	    mv_localKnots(idx) = mv_knots(i);
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

		mv_pts(i) = (1. - alpha) * mv_pts(i-1)+alpha*mv_pts(i);
			
	    }
	});
	Kokkos::deep_copy(pts, mv_pts);
	return pts(order_t-1);
	
    }





	//Accessors
    Kokkos::View<int*, MemSpace> getOrder() const {return order;}
    Kokkos::View<double*, MemSpace> getCtrlPts() const {return ctrlPts;}
    Kokkos::View<double*, MemSpace> getKnots() const {return knots;}
    Kokkos::View<int*, MemSpace> getCPOffset() const {return cPOffset;}
    Kokkos::View<int*, MemSpace> getKnotsOffset() const {return knotsOffset;}
    Kokkos::View<double*, MemSpace> get1stD() const {return ctrlPts_1stD;}
    Kokkos::View<double*, MemSpace> get2ndD() const {return ctrlPts_2ndD;}


    int getCPOffsetSize() const {return cPOffset.extent(0);}
    int getKnotsOffsetSize() const {return knotsOffset.extent(0);}
    int getCtrlPtsSize() const {return ctrlPts.extent(0);}
    int getKnotsSize() const {return knots.extent(0);}
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
	Kokkos::View<double*, MemSpace> ctrlPts_1stDV("ctrlPts1Derivative", ctrlPts.extent(0)-(2*cPOffset.extent(0)));
	Kokkos::View<int*, MemSpace> cP1stDOffsetV("ctrlPts1DerivativeOffset", cPOffset.extent(0));
	Kokkos::View<int*, MemSpace> cP2ndDOffsetV("ctrlPts2DerivativeOffset", cPOffset.extent(0));
        
	auto host_cP2ndDOffsetV = Kokkos::create_mirror_view(cP2ndDOffsetV);

	auto host_cP1stDOffsetV = Kokkos::create_mirror_view(cP1stDOffsetV);
	
	std::cout << "offset views and 1stctrlPts views construct" << std::endl;

	//Adjust the offset to include be 1 & 2 less than the ctrlPtsOffset

	auto host_cPOffset = Kokkos::create_mirror_view(cPOffset);
	Kokkos::deep_copy(host_cPOffset, cPOffset);
	if (host_cPOffset.extent(0) > 1) {
	    Kokkos::deep_copy(host_cPOffset, cPOffset);
	    for (int i = 1; i < cPOffset.extent(0); i++) {
	        host_cP1stDOffsetV(i) = host_cPOffset(i)-2;
	        host_cP2ndDOffsetV(i) = host_cPOffset(i)-4;
	        std::cout << "1st: " << host_cP1stDOffsetV(i);
	        std::cout << ", 2nd: " << host_cP2ndDOffsetV(i) << std::endl;
	    }
	}
	else {
	    host_cP1stDOffsetV(0) = host_cPOffset(0)-2;
	    host_cP2ndDOffsetV(0) = host_cPOffset(0)-4;
	}

	std::cout << "offset view populated" << std::endl;

	//We need to partition the x and y while we calculate the coefficient
	int idx = 0;
	int oidx = 0;

	auto host_knots = Kokkos::create_mirror_view(knots);
	auto host_order = Kokkos::create_mirror_view(order);
	auto host_ctrlPts = Kokkos::create_mirror_view(ctrlPts);
	auto host_ctrlPts_1stDV = Kokkos::create_mirror_view(ctrlPts_1stDV);
	
	//Copy these value over
	Kokkos::deep_copy(host_knots, knots);
	Kokkos::deep_copy(host_order, order);
	Kokkos::deep_copy(host_ctrlPts, ctrlPts);

	std::cout << "Copied ctrlPts, knots, order" << std::endl;
	//We need to keep track of the indicies of the splines
	for (int i = 2; i < host_ctrlPts.extent(0); i++) {
	    //Even indices stores all the x coordinates
	    //Odd indices stores all the y coordinates
	    if (i == host_cPOffset(idx)) {
	        //This is the start of the next spline
		oidx++;	//Advance to get the next spline order
		idx++;	//Advace offset to next spline cut off

	    }
	    double delta = double (host_order(oidx)-1)/(host_knots(i+host_order(oidx)-1)-host_knots(i));
	    host_ctrlPts_1stDV(i-2) = (host_ctrlPts(i)-host_ctrlPts(i-2))*delta;
	}

	for (int i = 0; i < host_ctrlPts_1stDV.extent(0); i++) {
            std::cout << "1stD Coeff: "<< host_ctrlPts_1stDV(i) << std::endl;
	}

	std::cout << "1stD coeff populated" << std::endl;

	Kokkos::View<double*, MemSpace> ctrlPts_2ndDV("ctrlPts2Derivative", ctrlPts.extent(0)-2*cPOffset.extent(0));
        auto host_ctrlPts_2ndDV = Kokkos::create_mirror_view(ctrlPts_2ndDV);
	idx = 0;
	oidx = 0;
	for (int i = 2; i < host_ctrlPts_1stDV.extent(0); i++) {
            if (i == host_cPOffset(idx)) {
		oidx++;
		idx++;
	    }
	    double delta = double(host_order(oidx)-2)/(host_knots(i+host_order(oidx)-1) - host_knots(i+1));
	    host_ctrlPts_2ndDV(i-2) = (host_ctrlPts(i) - host_ctrlPts(i-2))*delta;
	}

	Kokkos::deep_copy(cP1stDOffsetV, host_cP1stDOffsetV);
	Kokkos::deep_copy(cP2ndDOffsetV, host_cP2ndDOffsetV);
	Kokkos::deep_copy(ctrlPts_1stDV, host_ctrlPts_1stDV);
	Kokkos::deep_copy(ctrlPts_2ndDV, host_ctrlPts_2ndDV);

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

	//The 1st and 2nd derivatives
	Kokkos::View<double*, MemSpace> ctrlPts_1stD;
	Kokkos::View<int*, MemSpace> cP1stDOffset;
	Kokkos::View<double*, MemSpace> ctrlPts_2ndD;
	Kokkos::View<int*, MemSpace> cP2ndDOffset;

};
#endif
