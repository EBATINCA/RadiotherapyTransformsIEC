# RadiotherapyTransformsIEC
C++ library based on VTK for provinding affine transform matrices corresponding to the IEC 61217 standard. It is used in SlicerRT for e.g. the Room's Eye View, but can be used by any other external application.

[Check out our Wiki](https://github.com/EBATINCA/RadiotherapyTransformsIEC/wiki/IEC-coordinate-systems-summary) to visualize the translations and rotations defined by IEC 61217 standard.

## How to build
- First, be sure to have installed from source (or from binaries or package manager) a version of VTK >= 9.2
- `cd /opt/` (for example)
- `git clone git@github.com:EBATINCA/RadiotherapyTransformsIEC.git`
- `cd RadiotherapyTransformsIEC && mkdir build && cd build`
- `cmake -DVTK_DIR=/opt/VTK-9.3.1/install/lib/cmake/vtk-9.3/ ..` (VTK_DIR must be replaced with the path where you installed, or left away if system-wide install)
- `make`

## How to include library from external CMake projects

### Custom library
One way is to use the [External Project formalism](https://cmake.org/cmake/help/latest/module/ExternalProject.html) and then `target_link_libraries(... vtkIECTransformLogic)`. If you've already built by yourself the `vtkIEC` library in folder `/opt/RadiotherapyTransformsIEC/build/`, then add the following in your CMakeLists instead:

```
if(NOT DEFINED vtkIECTransformLogic_DIR)
set (vtkIECTransformLogic_DIR "/opt/RadiotherapyTransformsIEC/build/" CACHE STRING "IEC 61217 build directory containing vtkIECTransformLogicConfig.cmake")
endif()
find_package(vtkIECTransformLogic)
...
target_link_libraries(yourTarget PUBLIC vtkIECTransformLogic)
```

### From SlicerRt
There are two ways, assuming you have cloned SlicerRt repository in `/opt/SlicerRT_src` and have built Slicer in `/opt/Slicer-build`:
- `cmake ../SlicerRT_src/ -DSlicer_DIR=/opt/Slicer_bld/Slicer-build` (it will autodownload vtkIEC from git and build it in the process)
- `cmake ../SlicerRT_src/ -DSlicer_DIR=/opt/Slicer_bld/Slicer-build -DvtkIECTransformLogic_DIR=/opt/RadiotherapyTransformsIEC/build -DSlicerRT_SUPERBUILD:BOOL=OFF` (it will link instead to the custom-user build)

## Acknowledgment
This repository was partly supported by Conselleria de Educación, Investigación, Cultura y Deporte (Generalitat Valenciana), Spain under grant number CDEIGENT/2019/011.

## Contributors
- Csaba Pinter (Ebatinca)
- Márk Gerencsér (Ebatinca)
- Fernando Hueso-González (IFIC, CSIC - Universitat de València)
