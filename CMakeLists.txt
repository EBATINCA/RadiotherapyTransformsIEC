cmake_minimum_required (VERSION 3.21)

project (vtkIECTransformLogic)

set(vtkIECTransformLogic_SRCS
  ../vtkIECTransformLogic/src/vtkIECTransformLogic.cxx
  ../vtkIECTransformLogic/src/vtkIECTransformLogic.h
  )

add_library(vtkIECTransformLogic ${vtkIECTransformLogic_SRCS})

target_include_directories(vtkIECTransformLogic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

find_package(VTK)

target_link_libraries(vtkIECTransformLogic ${VTK_LIBRARIES})
