# simMeshGen
mesh generation with Simmetrix SimModSuite for ice sheets

## build on SCOREC rhel9

env setup

```
module use /opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/
module load gcc/12.3.0-iil3lno mpich/4.1.1-xpoyz4t 
module load cmake
module load simmetrix-simmodsuite/2023.1-230907dev-jquex4z
```

build

```
cmake -S simLandIceMeshGen -B buildSimLandIceMeshGen
cmake --build buildSimLandIceMeshGen
```

## run

```
./landIceMeshGen /path/to/jigsaw/geom.msh outputPrefix
```

The resulting mesh will be extrememly coarse with a lot of spider webs near the
grounding line.  Using SimModeler and setting a relative size attribute of 0.5
on all the grounding line edges will give a reasonable mesh of about 1M
elements.

