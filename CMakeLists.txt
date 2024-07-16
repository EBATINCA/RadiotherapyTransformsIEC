project(vtkIECTransformLogic)

set(KIT ${PROJECT_NAME})

set(${KIT}_INCLUDE_DIRECTORIES
  #${SlicerRtCommon_INCLUDE_DIRS}
  #${vtkSlicerSubjectHierarchyModuleLogic_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

set(${KIT}_SRCS
  vtkIECTransformLogic.cxx
  vtkIECTransformLogic.h
  )

set(${KIT}_TARGET_LIBRARIES
  #vtkSlicerRtCommon
  #vtkSlicerSubjectHierarchyModuleLogic
  ${VTK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  #EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )
  
