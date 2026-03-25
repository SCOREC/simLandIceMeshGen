#include<Kokkos_Core.hpp>

template<typename MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct CoordSet {
    CoordSet() {}
    CoordSet(const size_t numPts) : pts("CoordSet", numPts) {}
    using view_type = Kokkos::View<double*[2], MemSpace>;
    view_type pts;

    KOKKOS_INLINE_FUNCTION double x(int i) const { return pts(i, 0); }
    KOKKOS_INLINE_FUNCTION double y(int i) const { return pts(i, 1); }
    KOKKOS_INLINE_FUNCTION int    size()   const { return pts.extent(0); }
};

template<typename MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
struct BSpline {
  BSpline(const size_t n) : 
    ctrlPts(n), 
    ctrlPtsDeriv1(n-1), 
    ctrlPtsDeriv2(n-2) {};
  CoordSet<MemSpace> ctrlPts;
  CoordSet<MemSpace> ctrlPtsDeriv1;
  CoordSet<MemSpace> ctrlPtsDeriv2;
};

void basic() {
  const auto numPts = 3;
  CoordSet c(numPts);
  auto pts_h = Kokkos::create_mirror_view(c.pts);
  for(int i=0; i<numPts; i++) {
    pts_h(i,0) = i;
    pts_h(i,1) = 42 + 0.1*i;
  }
  Kokkos::deep_copy(pts_h, c.pts);

  Kokkos::parallel_for("foo", numPts, KOKKOS_LAMBDA(int i) {
    const double x = c.x(i) + c.y(i);
  });
}

void doBSplineStuff() {
  const auto numPts = 5;
  BSpline b(numPts);
  Kokkos::parallel_for("setPts", numPts, KOKKOS_LAMBDA(int i) {
    b.ctrlPts.pts(i,0) = i;
    b.ctrlPts.pts(i,1) = 42 + 0.1*i;
  });
  Kokkos::parallel_for("foo", numPts, KOKKOS_LAMBDA(int i) {
    const double x = b.ctrlPts.x(i) + b.ctrlPts.y(i);
  });
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  { //or use kokkos scope guard
    basic();
    doBSplineStuff();
  }
  Kokkos::finalize();
  return 0;
}

