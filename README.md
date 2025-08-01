# simMeshGen
mesh generation with Simmetrix SimModSuite for ice sheets

The input '.msh' files are expected to list all edges in counter clockwise
order.  The first four vertices and edges define a bounding box for the domain.
The remaining vertices and edges should define a closed loop that is entirely
within the bounding box (i.e., no intersections with bounding box).

## build on SCOREC rhel9

env setup

```
module use /opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/
module load gcc/12.3.0-iil3lno mpich/4.1.1-xpoyz4t 
module load cmake
module load simmetrix-simmodsuite/2025.0-250108dev-llxq6sk
module load openblas/0.3.23-wqm7iud
```

clone

```
git clone git@github.com:scorec/simLandIceMeshGen
```

build

```
cmake -S simLandIceMeshGen -B buildSimLandIceMeshGen
cmake --build buildSimLandIceMeshGen
```

## run tests

```
ctest
```


## run clang formatting

```
#!/bin/bash 
module use /opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/
module load llvm/16.0.2-s2prjub 

# run clang-format on hpp and cpp files
for i in `ls $PWD/*.h`; do 
  echo $i
  clang-format -i $i
done

for i in `ls $PWD/*.cc`; do 
  echo $i
  clang-format -i $i
done
```
