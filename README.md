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

## Installing / packaging

By default only the runtime library is installed. To also install the
development files (headers, exported CMake targets and the package config
file) — e.g. when building a package or a standalone install tree — turn on
`vtkIECTransformLogic_INSTALL_DEVELOPMENT`:

    cmake -DvtkIECTransformLogic_INSTALL_DEVELOPMENT=ON -DCMAKE_INSTALL_PREFIX=/path/to/prefix ../vtkIECTransformLogic
    make install

The install tree is relocatable. Downstream projects consume it with:

    find_package(vtkIECTransformLogic REQUIRED)
    target_link_libraries(MyTarget PRIVATE vtkIECTransformLogic)

The following variables let a parent project (typically a superbuild) override
the install layout; each defaults to a `GNUInstallDirs`-based value when unset:

| Variable                                   | Default                                                  |
|--------------------------------------------|----------------------------------------------------------|
| `vtkIECTransformLogic_INSTALL_DEVELOPMENT` | `OFF`                                                    |
| `vtkIECTransformLogic_INSTALL_BIN_DIR`     | `${CMAKE_INSTALL_BINDIR}` (`bin`)                        |
| `vtkIECTransformLogic_INSTALL_LIB_DIR`     | `${CMAKE_INSTALL_LIBDIR}` (`lib` or `lib64`)             |
| `vtkIECTransformLogic_INSTALL_INCLUDE_DIR` | `${CMAKE_INSTALL_INCLUDEDIR}/vtkIECTransformLogic`       |
| `vtkIECTransformLogic_INSTALL_CMAKE_DIR`   | `${CMAKE_INSTALL_LIBDIR}/cmake/vtkIECTransformLogic`     |


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
Alternatively, if you installed the library with `-DvtkIECTransformLogic_INSTALL_DEVELOPMENT=ON` into `/opt/RadiotherapyTransformsIEC/install/`, point `vtkIECTransformLogic_DIR` to the installed CMake package directory (or pass `-DCMAKE_PREFIX_PATH=/opt/RadiotherapyTransformsIEC/install` when configuring):

```cmake
if(NOT DEFINED vtkIECTransformLogic_DIR)
  set(vtkIECTransformLogic_DIR "/opt/RadiotherapyTransformsIEC/install/lib/cmake/vtkIECTransformLogic/" CACHE STRING "IEC 61217 install directory containing vtkIECTransformLogicConfig.cmake")
endif()
find_package(vtkIECTransformLogic REQUIRED)
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
- Declan Garvey (IFIC, CSIC - Universitat de València)
