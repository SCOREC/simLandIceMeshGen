#ifndef BSPLINEKOKKOS2D_H

#define BSPLINEKOKKOS2D_H

#include<Kokkos_Core.hpp>
#include<vector>
#include<iostream>
#include<string>

template<typename ExecutionSpace>

class BSplineKokkos2D {
public:
    using MemSpace = typename ExecutionSpace::memory_space;
    
    BSplineKokkos2D(int order_p, std::vector<double>& ctrlPts_x, std::vector<double>& ctrlPts_y, std::vector<double>& knotsI) {
	Kokkos::View<int*, MemSpace> orderV("Orders", 1);
	auto host_orderV = Kokkos::create_mirror_view(orderV);
	host_orderV(0) = order_p;
	order = orderV;
	Kokkos::deep_copy(order, host_orderV);

	Kokkos::View<double*[2], MemSpace> ctrlPtsV("ctrlPts", ctrlPts_x.size());
	auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
	for (int i = 0; i < ctrlPts_x.size(); i++) {
	    host_ctrlPtsV(i, 0) = ctrlPts_x[i];
	    host_ctrlPtsV(i, 1) = ctrlPts_y[i];
	}
	ctrlPts = ctrlPtsV;
	Kokkos::deep_copy(ctrlPts, host_ctrlPtsV);

	Kokkos::View<int*, MemSpace> cpOffsetV("cpOffset", 1);
	auto host_cpOffsetV = Kokkos::create_mirror_view(cpOffsetV);
	host_cpOffsetV(0) = ctrlPts_x.size();
	host_cpOffsetV(0) = ctrlPts_x.size();

	cpOffset = cpOffsetV;
	Kokkos::deep_copy(cpOffset, host_cpOffsetV);

	Kokkos::View<double*, MemSpace> knotsV("knots", knotsI.size());
	auto host_knotsV = Kokkos::create_mirror_view(knotsV);
	for (int i = 0; i < knotsI.size(); i++) {
	    host_knotsV(i) = knotsI[i];
	}
	knots = knotsV;
	Kokkos::deep_copy(knots, host_knotsV);
    
	Kokkos::View<int*, MemSpace> knotsOffsetV("knotsOffset", 1);
	auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
	host_knotsOffsetV(0) = knotsI.size();
	knotsOffset = knotsOffsetV;
	Kokkos::deep_copy(knotsOffset, host_knotsOffsetV);

	//TO DO: Implement 1st&2nd Coef function
	calculateDerivCoeff();

    }

    BSplineKokkos2D(std::vector<BSplineKokkos2D>& multiSplines) {
        //Loop over all the splines to know how much space to allocate
	int orderSize = 0;
	int ctrlPtsSize = 0;
	int cpOffsetSize = 0;
	int knotSize = 0;
	int knotsOffsetSize = 0;

	for (int i = 0; i < multiSplines.size(); i++) {
	    orderSize += multiSplines[i].getOrder().extent(0);
	    ctrlPtsSize += multiSplines[i].getCtrlPts().extent(0);
	    cpOffsetSize += multiSplines[i].getCPOffset().extent(0);
	    knotSize += multiSplines[i].getKnots().extent(0);
	    knotsOffsetSize += multiSplines[i].getKnotsOffset().extent(0);
	}

	//Populate the views
	Kokkos::View<int*, MemSpace> orderV ("order", orderSize);
	Kokkos::View<double*[2], MemSpace> ctrlPtsV("ctrlPts", ctrlPtsSize);
	Kokkos::View<int*, MemSpace> cPOffsetV("ctrlPtsOffset", cpOffsetSize);
	Kokkos::View<double*, MemSpace> knotsV("knots", knotSize);
	Kokkos::View<int*, MemSpace> knotsOffsetV("knotsOffset", knotsOffsetSize);

	//Create the mirror view on host to move the data over
	auto mvOrderV = Kokkos::create_mirror_view(orderV);
	auto mvCtrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
	auto mvCPOffsetV = Kokkos::create_mirror_view(cPOffsetV);
	auto mvKnotsV = Kokkos::create_mirror_view(knotsV);
	auto mvKnotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
	//Copy the data over
	int oidx = 0;
	int cPidx = 0;
	int cOidx = 0;
	int kidx = 0;
	int kOidx = 0;
	for (int i = 0; i < multiSplines.size(); i++) {
	    //Populate orders	
	    Kokkos::View<int*, MemSpace> intView = multiSplines[i].getOrder();
	    auto mvIntView = Kokkos::create_mirror_view(intView);
	    Kokkos::deep_copy(mvIntView, intView);
	    for (int j = 0; j < intView.extent(0); j++) {
		mvOrderV(oidx+j) = mvIntView(j);
	    }
	    oidx += mvIntView.extent(0);

	    //Populate ctrlPts
	    Kokkos::View<double*[2], MemSpace> double2DView = multiSplines[i].getCtrlPts();
	    auto mvDouble2DView = Kokkos::create_mirror_view(double2DView);
	    Kokkos::deep_copy(mvDouble2DView, double2DView);
	    for (int j = 0; j < mvDouble2DView.extent(0); j++) {
	        mvCtrlPtsV(cPidx+j, 0) = mvDouble2DView(j, 0);
		mvCtrlPtsV(cPidx+j, 1) = mvDouble2DView(j, 1);
	    }
	    cPidx += mvDouble2DView.extent(0);

	    //Populate ctrlPtsOffset
	    intView = multiSplines[i].getCPOffset();
	    mvIntView = Kokkos::create_mirror_view(intView);
	    Kokkos::deep_copy(mvIntView, intView);

	    for (int j = 0; j < intView.extent(0); j++) {
		mvCPOffsetV(cOidx + j) = intView(j);
	    }
	    cOidx += intView.extent(0);

	    //Populate knots
	    Kokkos::View<double*, MemSpace> doubleView = multiSplines[i].getKnots();
	    auto mvDoubleView = Kokkos::create_mirror_view(doubleView);
	    Kokkos::deep_copy(mvDoubleView, doubleView);

	    for (int j = 0; j < doubleView.extent(0); j++) {
		mvKnotsV(kidx+j) = doubleView(j);
	    }
	    kidx += doubleView.extent(0);

	    //Populate knots offset
	    intView = multiSplines[i].getKnotsOffset();
	    mvIntView = Kokkos::create_mirror_view(intView);
	    Kokkos::deep_copy(mvIntView, intView);

	    for (int j = 0; j < intView.extent(0); j++) {
		mvKnotsOffsetV(kOidx+j) = mvIntView(j);
	    }
	    kOidx += mvIntView.extent(0);
	}

	//Copy the data on host to device
	order = orderV;
	ctrlPts = ctrlPtsV;
	cpOffset = cPOffsetV;
	knots = knotsV;
	knotsOffset = knotsOffsetV;

	Kokkos::deep_copy(order, mvOrderV);
	Kokkos::deep_copy(ctrlPts, mvCtrlPtsV);
	Kokkos::deep_copy(cpOffset, mvCPOffsetV);
	Kokkos::deep_copy(knots, mvKnotsV);
	Kokkos::deep_copy(knotsOffset, mvKnotsOffsetV);
	
	//TO DO: Implement 1st&2nd Coef function
	calculateDerivCoeff();

    }

    void calculateDerivCoeff() {
        //Allocate the views we need
	Kokkos::View<double*[2], MemSpace> ctrlPts1stDV("ctrlPts1stDeriv", ctrlPts.extent(0)-cpOffset.extent(0));
	Kokkos::View<int*, MemSpace> cP1stDOffsetV("cP1stDOffset", cpOffset.extent(0));
	Kokkos::View<int*, MemSpace> cP2ndDOffsetV("cP2ndDOffset", cpOffset.extent(0));
	
	//Set up the offset views
	auto mvCP1stDOffsetV = Kokkos::create_mirror_view(cP1stDOffsetV);
	auto mvCP2ndDOffsetV = Kokkos::create_mirror_view(cP2ndDOffsetV);

	auto mvCPOffset = Kokkos::create_mirror_view(cpOffset);
	
	if (mvCPOffset.extent(0) > 1) {
	    Kokkos::deep_copy(mvCPOffset, cpOffset);
            for (int i = 1; i < mvCPOffset.extent(0); i++) {
	        mvCP1stDOffsetV(i) = mvCPOffset(i)-1;
	        mvCP2ndDOffsetV(i) = mvCPOffset(i)-2;
	    }
	} 
        else {
	    mvCP1stDOffsetV(0) = mvCPOffset(0)-1;
	    mvCP2ndDOffsetV(0) = mvCPOffset(0)-2;
	}

	auto mvOrder = Kokkos::create_mirror_view(order);
	auto mvKnots = Kokkos::create_mirror_view(knots);
	auto mvCtrlPts = Kokkos::create_mirror_view(ctrlPts);
	auto mvCtrlPts1stDV = Kokkos::create_mirror_view(ctrlPts1stDV);
	Kokkos::deep_copy(mvOrder, order);
	Kokkos::deep_copy(mvKnots, knots);
	Kokkos::deep_copy(mvCtrlPts, ctrlPts);

	//Calculate 1st derivative coef
	int offidx = 0;	//Offset index
	int oidx = 0; //Order index
	for (int i = 1; i < mvCtrlPts1stDV.extent(0); i++) {
	    //We need to check whether we are on the border for the next spline in our structure
	    if (i == mvCP1stDOffsetV(offidx)) {
		oidx++;
		offidx++;
		continue;	//Skip to the next
	    }
	    double delta = double(mvOrder(oidx-1)) / (mvKnots(i+mvOrder(oidx-1)) - mvKnots(i));
	    mvCtrlPts1stDV(i, 0) = (mvCtrlPts(i, 0) - mvCtrlPts(i-1, 0))*delta;
	    mvCtrlPts1stDV(i, 1) = (mvCtrlPts(i, 1) - mvCtrlPts(i-1, 1))*delta;
	}

	//Calculate 2nd derivative coef
	Kokkos::View<double*[2], MemSpace> ctrlPts2ndDV("ctrlPts2ndDeriv", ctrlPts.extent(0)-(2*cpOffset.extent(0)));
	auto mvCtrlPts2ndDV = Kokkos::create_mirror_view(ctrlPts2ndDV);

	for (int i = 1; i < ctrlPts1stDV.extent(0); i++) {
	    if (i == mvCP2ndDOffsetV(offidx)) {
		oidx++;
		offidx++;
		continue;
	    }
	    double delta = double(mvOrder(oidx)-2) / (mvKnots(i+mvOrder(oidx)-1) - mvKnots(i+1));
	    mvCtrlPts2ndDV(i, 0) = (mvCtrlPts(i, 0) - mvCtrlPts(i-1, 0)) * delta;
	    mvCtrlPts2ndDV(i, 1) = (mvCtrlPts(i, 1) - mvCtrlPts(i-1, 1)) * delta;
	}

	//Copy to device
	Kokkos::deep_copy(cP1stDOffsetV, mvCP1stDOffsetV);
	Kokkos::deep_copy(cP2ndDOffsetV, mvCP2ndDOffsetV);
	Kokkos::deep_copy(ctrlPts1stDV, mvCtrlPts1stDV);
	Kokkos::deep_copy(ctrlPts2ndDV, mvCtrlPts2ndDV);

	cp1stDOffset = cP1stDOffsetV;
	cp2ndDOffset = cP2ndDOffsetV;
	ctrlPts1stD = ctrlPts1stDV;
	ctrlPts2ndD = ctrlPts2ndDV;
    }

    //Accessors
    Kokkos::View<int*, MemSpace> getOrder() const {return order;}
    Kokkos::View<double*[2], MemSpace> getCtrlPts() const {return ctrlPts;}
    Kokkos::View<int*, MemSpace> getCPOffset() const {return cpOffset;}
    Kokkos::View<double*, MemSpace> getKnots() const {return knots;}
    Kokkos::View<int*, MemSpace> getKnotsOffset() const {return knotsOffset;}
    Kokkos::View<double*[2], MemSpace> getCP1stD() const {return ctrlPts1stD;}
    Kokkos::View<int*, MemSpace> getCP1stDOffset() const {return cp1stDOffset;}
    Kokkos::View<double*[2], MemSpace> getCP2ndD() const {return ctrlPts2ndD;}
    Kokkos::View<int*, MemSpace> getCP2ndDOffset() const {return cp2ndDOffset;}


private:
    Kokkos::View<int*, MemSpace> order;
    Kokkos::View<double*[2], MemSpace> ctrlPts;
    Kokkos::View<int*, MemSpace> cpOffset;
    Kokkos::View<double*, MemSpace> knots;
    Kokkos::View<int*, MemSpace> knotsOffset;
    Kokkos::View<double*[2], MemSpace> ctrlPts1stD;
    Kokkos::View<double*[2], MemSpace> ctrlPts2ndD;
    Kokkos::View<int*, MemSpace> cp1stDOffset;
    Kokkos::View<int*, MemSpace> cp2ndDOffset;

};

#endif
