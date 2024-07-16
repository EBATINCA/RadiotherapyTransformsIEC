#! Usage:
#! \code
#! vtkiectransformlogic_add_executable(args)
#! \endcode
#!
#! This macro adds an executable that uses UTF-8 code page
#! if vtkIECTransformLogic_USE_UTF8 is enabled.
#!
#! Linux and macOS already uses this code page by default but on Windows
#! it has to be set explicitly, by using a manifest file.
#!
#! \ingroup CMakeUtilities

function(vtkiectransformlogic_add_executable)

  if(vtkIECTransformLogic_USE_UTF8 AND WIN32)
    add_executable(${ARGN} ${vtkIECTransformLogic_CMAKE_DIR}/WindowsApplicationUseUtf8.manifest)
  else()
    add_executable(${ARGN})
  endif()

endfunction()
