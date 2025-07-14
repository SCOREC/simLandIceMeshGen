#include "landIceMeshGen.h"

int main(int argc, char **argv) {
  const auto debug = true;

  {
    const auto firstPt=0;
    const auto lastPt=40;
    std::map<int,int> pairs;
    pairs.insert({0,24});
    pairs.insert({10,12});
    pairs.insert({20,35});
    auto unNestedPairs = removeNestedSegments(pairs, firstPt, lastPt);
    assert(unNestedPairs.size() == 2);
  }
  return 0;
}

