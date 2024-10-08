/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/ 

/// vtkIECTransformLogicExport
///
/// The vtkIECTransformLogicExport captures some system differences between Unix
/// and Windows operating systems.

#ifndef __vtkIECTransformLogicExport_h
#define __vtkIECTransformLogicExport_h

//#include <vtkIECTransformLogicConfigure.h>

#if defined(WIN32) && !defined(VTKIECTRANSFORMLOGIC_STATIC)
#if defined(vtkIECTransformLogic_EXPORTS)
#define VTK_IEC_TRANSFORM_LOGIC_EXPORT __declspec( dllexport )
#else
#define VTK_IEC_TRANSFORM_LOGIC_EXPORT __declspec( dllimport )
#endif
#else
#define VTK_IEC_TRANSFORM_LOGIC_EXPORT
#endif

#endif
