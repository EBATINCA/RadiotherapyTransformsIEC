/*==============================================================================

  Copyright (c) Laboratory for Percutaneous Surgery (PerkLab)
  Queen's University, Kingston, ON, Canada. All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Vinith Suriyakumar and Csaba Pinter,
  PerkLab, Queen's University and was supported through the Applied Cancer
  Research Unit program of Cancer Care Ontario with funds provided by the
  Ontario Ministry of Health and Long-Term Care

==============================================================================*/

#ifndef __vtkIECTransformLogic_h
#define __vtkIECTransformLogic_h

//#include "vtkSlicerBeamsModuleLogicExport.h"
#include "../vtkIECTransformLogicExport.h"

// STD includes
#include <map>
#include <vector>
#include <cstdint>
#include <array>
#include <list>

// VTK includes
#include <vtkNew.h>
#include <vtkTransform.h>

class vtkGeneralTransform;

/// @brief Logic representing the IEC standard coordinate systems and transforms.
///
/// The IEC standard describes coordinate systems and a transform hierarchy to
/// represent objects taking part in an external beam radiation therapy delivery in 3D space.
/// With this logic class it is possible to get a transform from any defined coordinate
/// system to another by simply inputting the coordinate systems. The logic can observe an
/// RT beam node to get the geometrical parameters defining the state of the objects involved.
/// Image describing these coordinate frames:
/// https://github.com/SlicerRt/SlicerRtDoc/blob/master/technical/IEC%2061217-2002_CoordinateSystemsDiagram_HiRes.png
/// The RAS coordinate system is not part of IEC but we leave it here as a helper function since it's used e.g. in Slicer https://slicer.readthedocs.io/en/latest/user_guide/coordinate_systems.html

/*
                          "IEC 61217:2011 Hierarchy"

                   -------------------("f")---------------------
                   |                    |                      |
        ---------("g")                ("i")                  ("s")
        |          |                    |                      |
      ("r")      ("b")                ("o")                  ("e")
                   |                                           |
                 ("w")                                       ("t")
                                                               |
                                                    ---------("p")
                                                    |          |
                                                *("ras")    *("dp")
                                                               |
                                                            *("pi")

Legend:
  ("f") - Fixed reference system
  ("g") - GANTRY coordinate system
  ("b") - BEAM LIMITING DEVICE or DELINEATOR coordinate system
  ("w") - WEDGE FILTER coordinate system
  ("r") - X-RAY IMAGE RECEPTOR coordinate system
  ("s") - PATIENT SUPPORT coordinate system
  ("e") - Table top eccentric rotation coordinate system
  ("t") - Table top coordinate system
  ("p") - PATIENT coordinate system (LSA)
 *("dp")- PATIENT coordinate system in LPS (DICOM)
 *("pi")- Patient image regular grid coordinate system
*("ras")- PATIENT coordinate system in RAS (3D Slicer)
  ("i") - Imager coordinate system
  ("o") - Focus coordinate system
   "*"  - Not part of standard IEC coordinate frames
*/
/*
 IEC Patient (LSA) to LPS (DICOM) Patient transformation:
     Counter clockwise rotation around X-axis, angle = -90

                            1 0  0 0
     TransformationMatrix = 0 0 -1 0
                            0 1  0 0
                            0 0  0 1

 IEC Patient (LSA) to RAS (3D Slicer) Patient transformation:
     Counter clockwise rotation around X-axis, angle = -90
     Clockwise rotation around Z-axis, angle = 180

                            -1 0 0 0
     TransformationMatrix =  0 0 1 0
                             0 1 0 0
                             0 0 0 1
*/

class VTK_IEC_TRANSFORM_LOGIC_EXPORT vtkIECTransformLogic : public vtkObject
{
public:
  enum CoordinateSystemIdentifier
  {
    RAS = 0,
    FixedReference,
    Gantry,
    Collimator,
    LeftImagingPanel,
    RightImagingPanel,
    PatientSupportRotation, // Not part of the standard, but useful for visualization
    PatientSupport,
    TableTopEccentricRotation,
    TableTop,
    FlatPanel,
    WedgeFilter,
    Patient,
    DICOM,
    PatientImageRegularGrid,
    Imager,
    Focus,
    LastIECCoordinateFrame // Last index used for adding more coordinate systems externally
  };
  typedef std::list< CoordinateSystemIdentifier > CoordinateSystemsList;

public:
  static vtkIECTransformLogic *New();
  vtkTypeMacro(vtkIECTransformLogic, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /// @brief Update GantryToFixedReference transform based on gantry rotation angle about the Y-axis with possibility of setting gantry pitch angle which is not part of IEC standard but a DICOM addition
  /// The order of rotations starting from the fixed reference frame is pitch rotation followed by gantry (roll) rotation
  /// @warning it is assumed that the same order of rotations than for "table top" (https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.8.8.14.12.html) applies here, although it is not explicitly said in the DICOM standard (https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.8.8.25.6.html#sect_C.8.8.25.6.5)
  /// @see section 3.4 of IEC61217:2011, p.12 and [URL to Slicer3D forum]
  /// @param gantryRotationAngleDeg the rotation in degrees of the gantry frame counter clockwise around the Y-axis
  /// @param gantryPitchAngleDeg the rotation in degrees of the gantry frame counter clockwise around the X-axis (not part of IEC standard but a DICOM addition)
  void UpdateGantryToFixedReferenceTransform(double gantryRotationAngleDeg, double gantryPitchAngleDeg=0);
  /// @brief Update CollimatorToGantry transform based on collimator angle parameter and Z displacement of the collimator
  /// @see section 3.5 of IEC61217:2011, p.13
  /// @param collimatorRotationAngleDeg the rotation in degrees of the collimator frame counter clockwise around the Z-axis starting from the gantry frame
  /// @param bz displacement of the collimator frame origin from the gantry frame origin along the Z-axis
  void UpdateCollimatorToGantryTransform(double collimatorRotationAngleDeg, double bz = 0);
  /// @brief Update WedgeFilterToCollimator transform based on wedge filter angle parameter and z displacement of the wedge filter
  /// @see section 3.6 of IEC61217:2011, p.13
  /// @param wedgefilterRotationAngleDeg the rotation in degrees of the wedge filter frame counter clockwise around the Z-axis starting from the collimator frame
  /// @param wz displacement of the wedge filter frame origin from the collimator frame origin along the Z-axis
  void UpdateWedgeFilterToCollimatorTransform(double wedgefilterRotationAngleDeg, double wz = 0);
  /// @brief Update PatientSupportRotationToFixedReference transform based on patient support rotation parameter
  /// @see section 3.8 of IEC61217:2011, p.14
  /// @param patientSupportRotationAngleDeg the rotation in degrees of the patient support rotation frame counter clockwise around the Z-axis starting from the fixed reference frame
  void UpdatePatientSupportRotationToFixedReferenceTransform(double patientSupportRotationAngleDeg);
  /// @brief Update TableTopEccentricRotationToPatientSupportRotation transform based on eccentric angle rotation parameter and z displacement of the eccentric device
  /// Starting from the patient support frame displacement of the origin along the Y-axis is performed followed by the rotation around the Z-axis
  /// @see section 3.9 of IEC61217:2011, p.15
  /// @param TableTopEccentricRotationAngleDeg the rotation in degrees of the table top eccentric frame counter clockwise around the Z-axis
  /// @param ey displacement of the table top eccentric frame origin along the Y-axis
  void UpdateTableTopEccentricRotationToPatientSupportRotationTransform(double TableTopEccentricRotationAngleDeg, double ey = 0);
  /// @brief Update TableTopToTableTopEccentricRotation transform based on Table Top displacement(in X,Y,Z), Table Top pitch rotation(around X axis), and Table Top roll rotation(around Y axis) parameters
  /// Starting from the table top eccentric frame displacement of the origin is performed followed by table top pitch rotation around the X-axis and then table top roll rotation around the Y-axis
  /// @see section 3.10 of IEC61217:2011, p.15
  /// @param tx displacement of the table top frame origin along the X-axis
  /// @param ty displacement of the table top frame origin along the Y-axis
  /// @param tz displacement of the table top frame origin along the Z-axis
  /// @param tableTopPitchAngleDeg the rotation in degrees of the table top eccentric frame counter clockwise around the X-axis
  /// @param tableTopRollAngleDeg the rotation in degrees of the table top eccentric frame counter clockwise around the Y-axis
  void UpdateTableTopToTableTopEccentricRotationTransform(double tx, double ty, double tz, double tableTopPitchAngleDeg = 0, double tableTopRollAngleDeg = 0);
  /// @brief Update PatientToTableTop transform based on patient displacement(in X,Y,Z) and patient rotation around X(Psi), Y(Phi) and Z(Theta) parameters
  /// Starting from the table top frame displacement of the origin is performed followed by Psi rotation around the X-axis, Phi rotation around the Y-axis, and then Theta rotation around the Z-axis respectively
  /// @see section 3.11 of IEC61217:2011, p.16
  /// @param px displacement of the patient frame origin along the X-axis
  /// @param py displacement of the patient frame origin along the Y-axis
  /// @param pz displacement of the patient frame origin along the Z-axis
  /// @param patientPsiAngleDeg the rotation in degrees of the table top eccentric frame counter clockwise around the X-axis
  /// @param patientPhiAngleDeg the rotation in degrees of the table top eccentric frame counter clockwise around the Y-axis
  /// @param patientThetaAngleDeg the rotation in degrees of the table top eccentric frame counter clockwise around the Z-axis
  void UpdatePatientToTableTopTransform(double px, double py, double pz, double patientPsiAngleDeg = 0, double patientPhiAngleDeg = 0, double patientThetaAngleDeg = 0);
  /// @brief Update PatientImageRegularGrid transform based on the parameters
  /// Using the input variables the corresponding transformation matrix from DICOM (patient frame in LPS format) to patient image regular grid frame is constructed and applied.
  /// @see https://nipy.org/nibabel/dicom/dicom_orientation.html
  /// @param columnPixelSpacing distance between each pixel column within the image
  /// @param rowPixelSpacing distance between each pixel row within the image
  /// @param sliceDistance the spacing between consecutive images
  /// @param sx the X-displacement of the zeroth pixel position in the zeroth image from patient frame origin in the patient coordinate frame (DICOM LPS, not IEC)
  /// @param sy the Y-displacement of the zeroth pixel position in the zeroth image from patient frame origin in the patient coordinate frame (DICOM LPS, not IEC)
  /// @param sz the Z-displacement of the zeroth pixel position in the zeroth image from patient frame origin in the patient coordinate frame (DICOM LPS, not IEC)
  /// @param directionCosineXx x element of the row direction cosine of the image orientation (see https://dicom.innolitics.com/ciods/ct-image/image-plane/00200037)
  /// @param directionCosineXy y element of the row direction cosine of the image orientation (see https://dicom.innolitics.com/ciods/ct-image/image-plane/00200037)
  /// @param directionCosineXz z element of the row direction cosine of the image orientation (see https://dicom.innolitics.com/ciods/ct-image/image-plane/00200037)
  /// @param directionCosineYx x element of the column direction cosine of the image orientation (see https://dicom.innolitics.com/ciods/ct-image/image-plane/00200037)
  /// @param directionCosineYy y element of the column direction cosine of the image orientation (see https://dicom.innolitics.com/ciods/ct-image/image-plane/00200037)
  /// @param directionCosineYz z element of the column direction cosine of the image orientation (see https://dicom.innolitics.com/ciods/ct-image/image-plane/00200037)
  /// @note default orientation (for BIPED) corresponds to the X-pixel number increasing from the right to the left of the patient, and Y-pixel number increasing from anterior to posterior, Image slice index increases from inferior to superior (corresponds to DICOM's LPS coordinate system)
  void UpdatePatientImageRegularGridToDICOMTransform(double columnPixelSpacing, double rowPixelSpacing, double sliceDistance, double sx, double sy, double sz,
                                                     double directionCosineXx = 1, double directionCosineXy = 0, double directionCosineXz = 0,
                                                     double directionCosineYx = 0, double directionCosineYy = 1, double directionCosineYz = 0);

  /// @brief Get transform from one coordinate frame to another
  /// @param fromFrame start transformation from frame
  /// @param toFrame proceed transformation to frame
  /// @param outputTransform General (linear) transform matrix fromFrame -> toFrame. Matrix is correct if return flag is true.
  /// @param transformForBeam calculate dynamic transformation for beam model or other models
  /// (e.g. transformation from Patient RAS frame to Collimation frame: RAS -> Patient -> TableTop -> Eccentric -> Patient Support -> Fixed reference -> Gantry -> Collimator)  //TODO: Deprecated
  /// @return Success flag (false on any error)
  bool GetTransformBetween(CoordinateSystemIdentifier fromFrame, CoordinateSystemIdentifier toFrame,
    vtkGeneralTransform* outputTransform, bool transformForBeam=false);
  //TODO: See this transformForBeam part if still needed

  /// @brief Get coordinate system identifiers from root system down to frame system
  vtkTransform* GetElementaryTransformBetween(CoordinateSystemIdentifier fromFrame, CoordinateSystemIdentifier toFrame);

public:
  //std::map<CoordinateSystemIdentifier, std::string> GetCoordinateSystemsMap()
  //{
  //  return CoordinateSystemsMap;
  //}

  /// @brief Get name of transform node between two coordinate systems
  /// @return Transform node name between the specified coordinate frames.
  /// @note If IEC does not specify a transform between the given coordinate frames, then there will be no node with the returned name.
  std::string GetTransformNameBetween(CoordinateSystemIdentifier fromFrame, CoordinateSystemIdentifier toFrame);

public:
  std::vector<std::pair<CoordinateSystemIdentifier, CoordinateSystemIdentifier>> GetIECTransforms()
  {
    return this->IECTransforms;
  }

  /// @brief Converts a 3D vector containing the indices (e0,e1,e2) in each axis of a regular grid to a linear index position when the 3D data are stored in a linear flat array
  /// @note In the case of DICOM images stacked by slice position as regular grid, dim 0: slice index, dim 1: row index, dim 2: column index, (all starting from zero), since PixelData are stored with row-major ordering
  /// @note C ordering is used, ie the last dimension is contiguous in memory, then the second dimension, with the first dimension being most distant.
  /// @see generalized row-major ordering https://en.wikipedia.org/wiki/Row-_and_column-major_order#Address_calculation_in_general
  /// @param vectorizedIndex 3-component array consisting of the indices in each dimension (e0,e1,e2)
  /// @param nElems 3D array containing the number of elements in each dimension
  /// @return The linearised pixel index starting at zero as a single int
  /// @note This algorithm uses row-major ordering to calculate indices, as is the case with DICOM images
  static inline uint64_t VectorizedToLinearizedIndex(const std::array<uint16_t, 3>& vectorizedIndex, const std::array<uint16_t, 3>& nElems)
  {
    const uint16_t n0 = nElems[0];
    const uint16_t n1 = nElems[1];
    const uint16_t n2 = nElems[2];
    const uint16_t e0 = vectorizedIndex[0];
    const uint16_t e1 = vectorizedIndex[1];
    const uint16_t e2 = vectorizedIndex[2];
    if(e0 >= n0 || e1 >= n1 || e2 >= n2)
    {
      throw std::runtime_error("Indices (" + std::to_string(e0) + "," + std::to_string(e1) + "," + std::to_string(e2) + ") out of range (" + std::to_string(n0) + "," + std::to_string(n1) + "," + std::to_string(n2) + ")" );
    }
    return static_cast<uint64_t>(e0)*n1*n2 + static_cast<uint64_t>(e1)*n2 + static_cast<uint64_t>(e2);
  }

  /// @brief Converts a linear index position of data (stored as linear flat array) in a 3D regular grid to a 3D index (e0,e1,e2) in each axis of the grid
  /// @note In the case of DICOM images stacked by slice position as regular grid, dim 0: slice index, dim 1: row index, dim 2: column index, (all starting from zero), since PixelData are stored with row-major ordering
  /// @note C ordering is used, ie the last dimension is contiguous in memory, then the second dimension, with the first dimension being most distant.
  /// @see generalized row-major ordering https://en.wikipedia.org/wiki/Row-_and_column-major_order#Address_calculation_in_general
  /// @param linearizedIndex the linear index (starting at zero) to be converted
  /// @param nElems 3D array containing the number of elements in each dimension
  /// @return A 3-component array consisting of the indices in each dimension (e0,e1,e2)
  static inline std::array<uint16_t, 3> LinearizedToVectorizedIndex(const uint64_t linearizedIndex, const std::array<uint16_t, 3>& nElems)
  {
    const uint16_t n0 = nElems[0];
    const uint16_t n1 = nElems[1];
    const uint16_t n2 = nElems[2];
    const uint64_t totalElems = static_cast<uint64_t>(n0)*n1*n2;
    if(linearizedIndex >= totalElems)
    {
      throw std::runtime_error("Index (" + std::to_string(linearizedIndex) + ") out of range (totalElems = " + std::to_string(totalElems) + ")" );
    }
    const uint16_t e0 = static_cast<uint16_t>((linearizedIndex/n2)/n1);
    const uint16_t e1 = static_cast<uint16_t>((linearizedIndex/n2)%n1);
    const uint16_t e2 = static_cast<uint16_t>( linearizedIndex%n2);
    return std::array<uint16_t,3>{e0, e1, e2};
  }

  //std::map<CoordinateSystemIdentifier, std::list<CoordinateSystemIdentifier>> GetCoordinateSystemsHierarchy()
  //{
  //  return CoordinateSystemsHierarchy;
  //}

  //TODO: Use the frame names instead (nobody knows these variables outside)
  //vtkGetObjectMacro(FixedReferenceToRasTransform, vtkTransform);
  //vtkGetObjectMacro(GantryToFixedReferenceTransform, vtkTransform);
  //vtkGetObjectMacro(CollimatorToGantryTransform, vtkTransform);
  //vtkGetObjectMacro(WedgeFilterToCollimatorTransform, vtkTransform);
  //vtkGetObjectMacro(AdditionalCollimatorDevicesToCollimatorTransform, vtkTransform);
  //vtkGetObjectMacro(LeftImagingPanelToGantryTransform, vtkTransform);
  //vtkGetObjectMacro(RightImagingPanelToGantryTransform, vtkTransform);
  //vtkGetObjectMacro(PatientSupportRotationToFixedReferenceTransform, vtkTransform);
  //vtkGetObjectMacro(PatientSupportToPatientSupportRotationTransform, vtkTransform);
  //vtkGetObjectMacro(TableTopEccentricRotationToPatientSupportRotationTransform, vtkTransform);
  //vtkGetObjectMacro(TableTopToTableTopEccentricRotationTransform, vtkTransform);
  //vtkGetObjectMacro(PatientToTableTopTransform, vtkTransform);
  //vtkGetObjectMacro(DICOMToPatientTransform, vtkTransform);
  //vtkGetObjectMacro(PatientImageRegularGridToDICOMTransform, vtkTransform);
  //vtkGetObjectMacro(RasToPatientTransform, vtkTransform);
  //vtkGetObjectMacro(FlatPanelToGantryTransform, vtkTransform);

protected:
  /// @brief Get coordinate system identifiers from frame system up to root system
  /// Root system = FixedReference system
  /// @see IEC 61217:2011 hierarchy
  bool GetPathToRoot(CoordinateSystemIdentifier frame, CoordinateSystemsList& path);

  /// @brief Get coordinate system identifiers from root system down to frame system
  /// Root system = FixedReference system
  /// @see IEC 61217:2011 hierarchy
  bool GetPathFromRoot(CoordinateSystemIdentifier frame, CoordinateSystemsList& path);

protected:
  /// @brief Map from \sa CoordinateSystemIdentifier to coordinate system name. Used for getting transforms
  std::map<CoordinateSystemIdentifier, std::string> CoordinateSystemsMap;

  /// @brief List of IEC transforms
  std::vector< std::pair<CoordinateSystemIdentifier, CoordinateSystemIdentifier> > IECTransforms;

  /// @todo for hierarchy use tree with nodes, something like graph
  /// @brief Map of IEC coordinate systems hierarchy
  std::map< CoordinateSystemIdentifier, std::list< CoordinateSystemIdentifier > > CoordinateSystemsHierarchy;

protected:
  vtkNew<vtkTransform> FixedReferenceToRasTransform;
  vtkNew<vtkTransform> GantryToFixedReferenceTransform;
  vtkNew<vtkTransform> CollimatorToGantryTransform;
  vtkNew<vtkTransform> WedgeFilterToCollimatorTransform;
  vtkNew<vtkTransform> LeftImagingPanelToGantryTransform;
  vtkNew<vtkTransform> RightImagingPanelToGantryTransform;
  vtkNew<vtkTransform> PatientSupportRotationToFixedReferenceTransform;
  vtkNew<vtkTransform> PatientSupportToPatientSupportRotationTransform;
  vtkNew<vtkTransform> TableTopEccentricRotationToPatientSupportRotationTransform;
  vtkNew<vtkTransform> TableTopToTableTopEccentricRotationTransform;
  vtkNew<vtkTransform> PatientToTableTopTransform;
  vtkNew<vtkTransform> DICOMToPatientTransform;
  vtkNew<vtkTransform> PatientImageRegularGridToDICOMTransform;
  vtkNew<vtkTransform> RasToPatientTransform;
  vtkNew<vtkTransform> FlatPanelToGantryTransform;

  vtkNew<vtkTransform> GantryToFixedReferenceConcatenatedTransform;
  vtkNew<vtkTransform> CollimatorToGantryConcatenatedTransform;
  vtkNew<vtkTransform> WedgeFilterToCollimatorConcatenatedTransform;
  vtkNew<vtkTransform> LeftImagingPanelToGantryConcatenatedTransform;
  vtkNew<vtkTransform> RightImagingPanelToGantryConcatenatedTransform;
  vtkNew<vtkTransform> FlatPanelToGantryConcatenatedTransform;
  vtkNew<vtkTransform> PatientSupportRotationToFixedReferenceConcatenatedTransform;
  vtkNew<vtkTransform> PatientSupportToPatientSupportRotationConcatenatedTransform;
  vtkNew<vtkTransform> TableTopEccentricRotationToPatientSupportRotationConcatenatedTransform;
  vtkNew<vtkTransform> TableTopToTableEccentricRotationConcatenatedTransform;
  vtkNew<vtkTransform> PatientToTableTopConcatenatedTransform;
  vtkNew<vtkTransform> DICOMToPatientConcatenatedTransform;
  vtkNew<vtkTransform> PatientImageRegularGridToDICOMConcatenatedTransform;
  vtkNew<vtkTransform> RasToPatientConcatenatedTransform;

protected:
  vtkIECTransformLogic();
  ~vtkIECTransformLogic() override;

private:
  vtkIECTransformLogic(const vtkIECTransformLogic&) = delete;
  void operator=(const vtkIECTransformLogic&) = delete;

private:
  std::vector<vtkTransform*> ElementaryTransforms;
};

#endif
