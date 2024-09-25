#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)
#-----------------------------------------------------------------------------

project(vtkIECTransformLogic
  LANGUAGES CXX
  )

### Set release mode as default build type
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")# "MinSizeRel" "RelWithDebInfo"
endif()

set(vtkIECTransformLogic_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set(CMAKE_MODULE_PATH ${vtkIECTransformLogic_CMAKE_DIR} ${CMAKE_MODULE_PATH})

include(CMakeDependentOption)
include(vtkIECTransformLogicFunctionAddExecutable)

# --------------------------------------------------------------------------
# Dependencies
# --------------------------------------------------------------------------

#
# VTK
#
find_package(VTK 9.2 REQUIRED)
set(vtkIECTransformLogic_LIBS VTK::CommonTransforms)

# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------
if(NOT DEFINED BUILD_SHARED_LIBS)
  option(BUILD_SHARED_LIBS "Build with shared libraries." ON)
endif()

if(NOT DEFINED BUILD_TESTING)
  option(BUILD_TESTING "Test the project" OFF)
endif()

if(NOT DEFINED vtkIECTransformLogic_USE_UTF8)
  option(vtkIECTransformLogic_USE_UTF8 "Make applications use UTF-8 as code page." OFF)
endif()

if(NOT DEFINED vtkIECTransformLogic_WRAP_PYTHON)
  cmake_dependent_option(vtkIECTransformLogic_WRAP_PYTHON "Build python wrapping" OFF "VTK_WRAP_PYTHON" OFF)
endif()

if(NOT DEFINED vtkIECTransformLogic_LAUNCH_COMMAND)
  set(vtkIECTransformLogic_LAUNCH_COMMAND "" CACHE STRING "Command for setting up environment and running executables")
endif()

# --------------------------------------------------------------------------
# Configure headers
# --------------------------------------------------------------------------
set(configure_header_file vtkIECTransformLogicConfigure.h)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${configure_header_file}.in
  ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
  )

set(configure_export_header_file vtkIECTransformLogicExport.h)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${configure_export_header_file}
  ${CMAKE_CURRENT_BINARY_DIR}/${configure_export_header_file}
  )

set(MY_EXPORT_HEADER_PREFIX vtkIECTransformLogic)
set(MY_LIBRARY_EXPORT_DIRECTIVE VTK_IEC_TRANSFORM_LOGIC_EXPORT)
set(MY_LIBNAME vtkIECTransformLogic)

set(dynamicHeaders
  "${dynamicHeaders};${CMAKE_CURRENT_BINARY_DIR}/${MY_EXPORT_HEADER_PREFIX}Export.h")

#-----------------------------------------------------------------------------
# Source groups
# --------------------------------------------------------------------------
source_group("Generated" FILES
  ${dynamicHeaders}
  )

# --------------------------------------------------------------------------
# Install headers
# --------------------------------------------------------------------------

if(NOT DEFINED ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  set(${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT ON)
endif()

if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)

  if(NOT DEFINED ${PROJECT_NAME}_INSTALL_INCLUDE_DIR)
    set(${PROJECT_NAME}_INSTALL_INCLUDE_DIR include/${PROJECT_NAME})
  endif()

  file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  install(
    FILES ${headers} ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file} ${CMAKE_CURRENT_BINARY_DIR}/${configure_export_header_file}
    DESTINATION ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR} COMPONENT Development)
endif()

# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------
set(vtkIECTransformLogic_SRCS
  src/vtkIECTransformLogic.cxx
  src/vtkIECTransformLogic.h
)

# --------------------------------------------------------------------------
# Include dirs
# --------------------------------------------------------------------------
set(include_dirs
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )
#~ include_directories(BEFORE ${include_dirs})

# --------------------------------------------------------------------------
# Build library
# --------------------------------------------------------------------------
set(lib_name ${PROJECT_NAME})

set(srcs ${vtkIECTransformLogic_SRCS})
add_library(${lib_name} SHARED ${srcs}) #SHARED = THERE IS DLL FILE
target_include_directories(${lib_name} PUBLIC
  $<BUILD_INTERFACE:${include_dirs}>
  $<INSTALL_INTERFACE:include/${lib_name}>  # <prefix>/include/mylib
)
target_link_libraries(${lib_name} PUBLIC ${vtkIECTransformLogic_LIBS})

#~ message("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ")
#~ message(${CMAKE_BINARY_DIR})
#~ message(${CMAKE_BINARY_DIR})

# Optionally, specify the output directory for the DLL, LIB, PDB, etc.
#set_target_properties(${lib_name} PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY c/d/_Extensions/inner-build/lib/Slicer-5.7/qt-loadable-modules/   # For .dll
#    ARCHIVE_OUTPUT_DIRECTORY "c/d/_Extensions/inner-build/lib/Slicer-5.7/qt-loadable-modules/"   # For .lib
#    LIBRARY_OUTPUT_DIRECTORY "c/d/_Extensions/inner-build/lib/Slicer-5.7/qt-loadable-modules/"   # For .lib and .exp
#    PDB_OUTPUT_DIRECTORY c/d/_Extensions/inner-build/lib/Slicer-5.7/qt-loadable-modules/       # For .pdb
#)

# --------------------------------------------------------------------------
# Folder
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_FOLDER)
  set(${PROJECT_NAME}_FOLDER ${PROJECT_NAME})
endif()
if(NOT "${${PROJECT_NAME}_FOLDER}" STREQUAL "")
  set_target_properties(${lib_name} PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
endif()

# --------------------------------------------------------------------------
# Export target
# --------------------------------------------------------------------------

# Export targets for a build tree
if(NOT DEFINED ${PROJECT_NAME}_EXPORT_FILE)
  set(${PROJECT_NAME}_EXPORT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)
endif()
export(TARGETS ${lib_name} APPEND FILE ${${PROJECT_NAME}_EXPORT_FILE})

# --------------------------------------------------------------------------
# Install library
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_BIN_DIR)
  set(${PROJECT_NAME}_INSTALL_BIN_DIR bin)
endif()
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_LIB_DIR)
  set(${PROJECT_NAME}_INSTALL_LIB_DIR lib/${PROJECT_NAME})
endif()

install(TARGETS ${lib_name}
  RUNTIME DESTINATION ${CMAKE_BINARY_DIR} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${CMAKE_BINARY_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${CMAKE_BINARY_DIR} COMPONENT Development
  )

# --------------------------------------------------------------------------
# Python Wrapping
# --------------------------------------------------------------------------
if(vtkIECTransformLogic_WRAP_PYTHON)

  if(NOT DEFINED ${PROJECT_NAME}_INSTALL_PYTHON_MODULE_LIB_DIR)
    set(${PROJECT_NAME}_INSTALL_PYTHON_MODULE_LIB_DIR ${${PROJECT_NAME}_INSTALL_LIB_DIR})
  endif()
  if(NOT DEFINED ${PROJECT_NAME}_INSTALL_PYTHON_LIB_DIR)
    set(${PROJECT_NAME}_INSTALL_PYTHON_LIB_DIR ${${PROJECT_NAME}_INSTALL_LIB_DIR})
  endif()

  include(vtkMacroKitPythonWrap)

  vtkMacroKitPythonWrap(
    KIT_NAME ${lib_name}
    KIT_SRCS ${vtkIECTransformLogic_SRCS}
    KIT_INSTALL_BIN_DIR ${${PROJECT_NAME}_INSTALL_BIN_DIR}
    KIT_INSTALL_LIB_DIR ${${PROJECT_NAME}_INSTALL_PYTHON_LIB_DIR}
    KIT_MODULE_INSTALL_LIB_DIR ${${PROJECT_NAME}_INSTALL_PYTHON_MODULE_LIB_DIR}
    )
  # Export target
  export(TARGETS ${lib_name}Python APPEND FILE ${${PROJECT_NAME}_EXPORT_FILE})
  if(${VTK_VERSION} VERSION_LESS "8.90")
    export(TARGETS ${lib_name}PythonD APPEND FILE ${${PROJECT_NAME}_EXPORT_FILE})
  endif()
  # Folder
  if(NOT "${${PROJECT_NAME}_FOLDER}" STREQUAL "")
    set_target_properties(${lib_name}Python PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
    if(${VTK_VERSION} VERSION_LESS "8.90")
      set_target_properties(${lib_name}PythonD PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
    endif()
    if(TARGET ${lib_name}Hierarchy)
      set_target_properties(${lib_name}Hierarchy PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
    endif()
  endif()
endif()


# --------------------------------------------------------------------------
# Set INCLUDE_DIRS variable
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include dirs" FORCE)

# --------------------------------------------------------------------------
# Configure support file
# --------------------------------------------------------------------------

# Configure 'vtkIECTransformLogicConfig.cmake' for a build tree
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkIECTransformLogicConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/vtkIECTransformLogicConfig.cmake
  )

# Configuret vtkIECTransformLogicConfig.cmake for an install tree
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkIECTransformLogicInstallConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/install/vtkIECTransformLogicConfig.cmake
  )

# --------------------------------------------------------------------------
# Install configuration files
# --------------------------------------------------------------------------

# Set vtkIECTransformLogic_INSTALL_CMAKE_DIR
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_CMAKE_DIR)
  set(${PROJECT_NAME}_INSTALL_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
endif()

# Install vtkIECTransformLogicConfig.cmake
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/install/vtkIECTransformLogicConfig.cmake
  DESTINATION ${${PROJECT_NAME}_INSTALL_CMAKE_DIR}
  )

# Install vtkIECTransformLogic CMake files
install(
  FILES
  ${CMAKE_SOURCE_DIR}/CMake/vtkIECTransformLogicFunctionAddExecutable.cmake
  ${CMAKE_SOURCE_DIR}/CMake/vtkIECTransformLogicFunctionAddExecutable.cmake
  ${CMAKE_SOURCE_DIR}/CMake/vtkMacroKitPythonWrap.cmake
  ${CMAKE_SOURCE_DIR}/CMake/vtkWrapHierarchy.cmake
  ${CMAKE_SOURCE_DIR}/CMake/vtkWrapperInit.data.in
  ${CMAKE_SOURCE_DIR}/CMake/vtkWrapPython.cmake
  ${CMAKE_SOURCE_DIR}/CMake/WindowsApplicationUseUtf8.manifest
  DESTINATION ${${PROJECT_NAME}_INSTALL_CMAKE_DIR} COMPONENT Development
  )
