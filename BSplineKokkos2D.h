#ifndef BSPLINEKOKKOS_H
#define BSPLINEKOKKOS_H

#include<Kokkos_Core.hpp>
#include<vector>
#include<iostream>
#include<string>

template<typename ExecutionSpace>

class BSPlineKokkos2DView {
public:
    using MemSpace = typename ExecutionSpace::memory_space;

    BSplineKokkos2DView(int order_p, std::vector<double>& ctrlPts_x, std::vector<double>& ctrlPts_y, std::vector<double>& knotsI) {
        Kokkos::View<int*, MemSpace> orderV("Orders", 1);
	auto host_orderV = Kokkos::create_mirror_view(orderV);
	host_orderV(0) = order_p;
	order = orderV;
	Kokkos::deep_copy(order, host_orderV);

	Kokkos::View<double*[2], MemSpace> ctrlPtsV("ctrlPts", ctrlPts.size());
	auto host_ctrlPtsV = Kokkos::create_mirror_view(ctrlPtsV);
	for (int i = 0; i < ctrlPts_x.size(); i++) {
	    host_ctrlPtsV(i, 0) = ctrlPts_x(i);
	    host_ctrlPtsV(i, 1) = ctrlPts_y(i);
	}
	ctrlPts = ctrlPtsV;
	Kokkos::deep_copy(ctrlPts, host_ctrlPtsV);

	Kokkos::View<int*[2], MemSpace> cpOffsetV("cpOffset", 1);
	auto host_cpOffsetV = Kokkos::create_mirror_view(cpOffsetV);
	host_cpOffsetV(0, 0) = ctrlPts_x.size();
	host_cpOffsetV(0, 1) = ctrlPts_x.size();

	cpOffset = cpOffsetV;
	Kokkos::deep_copy(cpOffset, host_cpOffsetV);

	Kokkos::View<double*, MemSpace> knotsV("knots", knotsI.size());
	auto host_knotsV = Kokkos::create_mirror_view(knotsV);
	for (int i = 0; i < knotsI.size(); i++) {
	    host_knotsV(i) = knotsI[i];
	}
	knots = knotsV;
	Kokkos::deep_copy(knots, host_knotsV(i));
    
	Kokkos::View<int*, MemSpace> knotsOffsetV("knotsOffset", 1);
	auto host_knotsOffsetV = Kokkos::create_mirror_view(knotsOffsetV);
	host_knotsOffsetV(0) = knotsI.size();
	knotsOffset = knotsOffsetV;
	Kokkos::deep_copy(knotsOffset, host_knotsOffsetV);

	//Later we need to address the coef calculation situation
    }

    BSplineKokkos2DView(std::vector<BSplineKokkos2DView>& multiSplines) {
        //Loop over all the splines to know how much space to allocate
	int orderSize = 0;
	int ctrlPtsSize = 0;
	int cpOffsetSize = 0;
	int knotSize = 0;
	int knotsOffsetSize = 0;

	for (int i = 0; i < multiSplines.size(); i++) {
	    orderSize += multiSplines[i];
	}
    }

    //Accessors
    Kokkos::View<int*, MemSpace> getOrder() const {return order;}
    Kokkos::View<double*[2], MemSpace> getCtrlPts() const {return ctrlPts;}
    Kokkos::View<int*, MemSpace> getCPOffset() const {return cPOffset;}
    Kokkos::View<double*, MemSpace> getKnots() const {return knots;}
    Kokkos::View<int*, MemSpace> getKnotsOffset() const {return knotsOffset;}
    Kokkos::View<double*[2], MemSpace> getCP1stD() const {return ctrlPts1stD;}
    Kokkos::View<int*, MemSpace> getCP1stDOffset() const {return cp1stDOffset;}
    Kokkos::View<double*[2], MemSpace> getCP2ndD() const {return ctrlPts2ndD;}
    Kokkos::View<int*, MemSpace> getCP2ndDOffset() const {return cp2ndDOffset;}


private:
    Kokkos::View<int*, MemSpace> order;
    Kokkos::View<double*[2], MemSpace> ctrlPts;
    kokkos::View<int*, MemSpace> cpOffset;
    Kokkos::View<double*, MemSpace> knots;
    Kokkos::View<int*, MemSpace> knotsOffset;
    Kokkos::View<double*[2], MemSpace> ctrlPts1stD;
    Kokkos::View<double*[2], MemSpace> ctrlPts2ndD;
    Kokkos::View<int*, MemSpace> cp1stDOffset;
    Kokkos::View<int*, MemSpace> cp2ndDOffset;

};

#endif
