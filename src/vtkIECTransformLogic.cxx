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

// IEC Logic includes
#include "vtkIECTransformLogic.h"

// VTK includes
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkGeneralTransform.h>
#include <vtkTransform.h>

// STD includes
#include <algorithm>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkIECTransformLogic);

//-----------------------------------------------------------------------------
vtkIECTransformLogic::vtkIECTransformLogic()
{
  // Setup coordinate system ID to name map
  this->CoordinateSystemsMap.clear();
  this->CoordinateSystemsMap[RAS] = "Ras";
  this->CoordinateSystemsMap[FixedReference] = "FixedReference";
  this->CoordinateSystemsMap[Gantry] = "Gantry";
  this->CoordinateSystemsMap[Collimator] = "Collimator";
  this->CoordinateSystemsMap[LeftImagingPanel] = "LeftImagingPanel";
  this->CoordinateSystemsMap[RightImagingPanel] = "RightImagingPanel";
  this->CoordinateSystemsMap[PatientSupportRotation] = "PatientSupportRotation";
  this->CoordinateSystemsMap[PatientSupport] = "PatientSupport";
  this->CoordinateSystemsMap[TableTopEccentricRotation] = "TableTopEccentricRotation";
  this->CoordinateSystemsMap[TableTop] = "TableTop";
  this->CoordinateSystemsMap[FlatPanel] = "FlatPanel";
  this->CoordinateSystemsMap[WedgeFilter] = "WedgeFilter";
  this->CoordinateSystemsMap[Patient] = "Patient";
  this->CoordinateSystemsMap[DICOM] = "DICOM";
  this->CoordinateSystemsMap[PatientImageRegularGrid] = "PatientImageRegularGrid";

  this->IECTransforms.clear();
  this->IECTransforms.push_back(std::make_pair(FixedReference, RAS));
  this->IECTransforms.push_back(std::make_pair(Gantry, FixedReference));
  this->IECTransforms.push_back(std::make_pair(Collimator, Gantry));
  this->IECTransforms.push_back(std::make_pair(WedgeFilter, Collimator));
  this->IECTransforms.push_back(std::make_pair(LeftImagingPanel, Gantry));
  this->IECTransforms.push_back(std::make_pair(RightImagingPanel, Gantry));
  this->IECTransforms.push_back(std::make_pair(PatientSupportRotation, FixedReference)); // Rotation component of patient support transform
  this->IECTransforms.push_back(std::make_pair(PatientSupport, PatientSupportRotation)); // Scaling component of patient support transform
  this->IECTransforms.push_back(std::make_pair(TableTopEccentricRotation, PatientSupportRotation)); // NOTE: Currently not supported by REV
  this->IECTransforms.push_back(std::make_pair(TableTop, TableTopEccentricRotation));
  this->IECTransforms.push_back(std::make_pair(Patient, TableTop));
  this->IECTransforms.push_back(std::make_pair(DICOM, Patient));
  this->IECTransforms.push_back(std::make_pair(PatientImageRegularGrid, DICOM));
  this->IECTransforms.push_back(std::make_pair(RAS, Patient));
  this->IECTransforms.push_back(std::make_pair(FlatPanel, Gantry));

  // Set names for elementary transforms for discovery
  this->FixedReferenceToRasTransform->SetObjectName(this->GetTransformNameBetween(FixedReference, RAS).c_str());
  this->GantryToFixedReferenceTransform->SetObjectName(this->GetTransformNameBetween(Gantry, FixedReference).c_str());
  this->CollimatorToGantryTransform->SetObjectName(this->GetTransformNameBetween(Collimator, Gantry).c_str());
  this->WedgeFilterToCollimatorTransform->SetObjectName(this->GetTransformNameBetween(WedgeFilter, Collimator).c_str());
  this->LeftImagingPanelToGantryTransform->SetObjectName(this->GetTransformNameBetween(LeftImagingPanel, Gantry).c_str());
  this->RightImagingPanelToGantryTransform->SetObjectName(this->GetTransformNameBetween(RightImagingPanel, Gantry).c_str());
  this->PatientSupportRotationToFixedReferenceTransform->SetObjectName(this->GetTransformNameBetween(PatientSupportRotation, FixedReference).c_str());
  this->PatientSupportToPatientSupportRotationTransform->SetObjectName(this->GetTransformNameBetween(PatientSupport, PatientSupportRotation).c_str());
  this->TableTopEccentricRotationToPatientSupportRotationTransform->SetObjectName(this->GetTransformNameBetween(TableTopEccentricRotation, PatientSupportRotation).c_str());
  this->TableTopToTableTopEccentricRotationTransform->SetObjectName(this->GetTransformNameBetween(TableTop, TableTopEccentricRotation).c_str());
  this->PatientToTableTopTransform->SetObjectName(this->GetTransformNameBetween(Patient, TableTop).c_str());
  this->DICOMToPatientTransform->SetObjectName(this->GetTransformNameBetween(DICOM, Patient).c_str());
  this->PatientImageRegularGridToDICOMTransform->SetObjectName(this->GetTransformNameBetween(PatientImageRegularGrid, DICOM).c_str());
  this->RasToPatientTransform->SetObjectName(this->GetTransformNameBetween(RAS, Patient).c_str());
  this->FlatPanelToGantryTransform->SetObjectName(this->GetTransformNameBetween(FlatPanel, Gantry).c_str());

  // Build list of elementary transforms for discovery by name
  this->ElementaryTransforms.push_back(this->FixedReferenceToRasTransform);
  this->ElementaryTransforms.push_back(this->GantryToFixedReferenceTransform);
  this->ElementaryTransforms.push_back(this->CollimatorToGantryTransform);
  this->ElementaryTransforms.push_back(this->WedgeFilterToCollimatorTransform);
  this->ElementaryTransforms.push_back(this->LeftImagingPanelToGantryTransform);
  this->ElementaryTransforms.push_back(this->RightImagingPanelToGantryTransform);
  this->ElementaryTransforms.push_back(this->PatientSupportRotationToFixedReferenceTransform);
  this->ElementaryTransforms.push_back(this->PatientSupportToPatientSupportRotationTransform);
  this->ElementaryTransforms.push_back(this->TableTopEccentricRotationToPatientSupportRotationTransform);
  this->ElementaryTransforms.push_back(this->TableTopToTableTopEccentricRotationTransform);
  this->ElementaryTransforms.push_back(this->PatientToTableTopTransform);
  this->ElementaryTransforms.push_back(this->DICOMToPatientTransform);
  this->ElementaryTransforms.push_back(this->PatientImageRegularGridToDICOMTransform);
  this->ElementaryTransforms.push_back(this->RasToPatientTransform);
  this->ElementaryTransforms.push_back(this->FlatPanelToGantryTransform);

  // Define the transform hierarchy
  this->CoordinateSystemsHierarchy.clear();
  // key - parent, value - children
  this->CoordinateSystemsHierarchy[FixedReference] = { Gantry, PatientSupportRotation };
  this->CoordinateSystemsHierarchy[Gantry] = { Collimator, LeftImagingPanel, RightImagingPanel, FlatPanel };
  this->CoordinateSystemsHierarchy[Collimator] = { WedgeFilter };
  this->CoordinateSystemsHierarchy[PatientSupportRotation] = { PatientSupport, TableTopEccentricRotation };
  this->CoordinateSystemsHierarchy[TableTopEccentricRotation] = { TableTop };
  this->CoordinateSystemsHierarchy[TableTop] = { Patient };
  this->CoordinateSystemsHierarchy[Patient] = { DICOM, RAS };
  this->CoordinateSystemsHierarchy[DICOM] = { PatientImageRegularGrid };

  // Build transformations that are not identity by default
  // define transformation matrix from the DICOM patient frame(LPS) to IEC patient frame(LSA) which is equivalent to a rotation around the X-axis +90deg counter clockwise
  double dicomToPatientTransformationMatrix[16] = {1, 0,0,0,
                                                   0, 0,1,0,
                                                   0,-1,0,0,
                                                   0, 0,0,1};
  this->DICOMToPatientTransform->Concatenate(dicomToPatientTransformationMatrix);

  //
  // Build concatenated transform hierarchy
  //
  this->GantryToFixedReferenceConcatenatedTransform->Concatenate(this->FixedReferenceToRasTransform);
  this->GantryToFixedReferenceConcatenatedTransform->Concatenate(this->GantryToFixedReferenceTransform);

  this->CollimatorToGantryConcatenatedTransform->Concatenate(this->GantryToFixedReferenceConcatenatedTransform);
  this->CollimatorToGantryConcatenatedTransform->Concatenate(this->CollimatorToGantryTransform);

  this->WedgeFilterToCollimatorConcatenatedTransform->Concatenate(this->CollimatorToGantryConcatenatedTransform);
  this->WedgeFilterToCollimatorConcatenatedTransform->Concatenate(this->WedgeFilterToCollimatorTransform);

  this->LeftImagingPanelToGantryConcatenatedTransform->Concatenate(this->CollimatorToGantryConcatenatedTransform);
  this->LeftImagingPanelToGantryConcatenatedTransform->Concatenate(this->LeftImagingPanelToGantryTransform);

  this->RightImagingPanelToGantryConcatenatedTransform->Concatenate(this->CollimatorToGantryConcatenatedTransform);
  this->RightImagingPanelToGantryConcatenatedTransform->Concatenate(this->RightImagingPanelToGantryTransform);

  this->FlatPanelToGantryConcatenatedTransform->Concatenate(this->CollimatorToGantryConcatenatedTransform);
  this->FlatPanelToGantryConcatenatedTransform->Concatenate(this->FlatPanelToGantryTransform);

  this->PatientSupportRotationToFixedReferenceConcatenatedTransform->Concatenate(this->FixedReferenceToRasTransform);
  this->PatientSupportRotationToFixedReferenceConcatenatedTransform->Concatenate(this->PatientSupportRotationToFixedReferenceTransform);

  this->PatientSupportToPatientSupportRotationConcatenatedTransform->Concatenate(this->PatientSupportRotationToFixedReferenceConcatenatedTransform);
  this->PatientSupportToPatientSupportRotationConcatenatedTransform->Concatenate(this->PatientSupportToPatientSupportRotationTransform);

  this->TableTopEccentricRotationToPatientSupportRotationConcatenatedTransform->Concatenate(this->PatientSupportRotationToFixedReferenceConcatenatedTransform);
  this->TableTopEccentricRotationToPatientSupportRotationConcatenatedTransform->Concatenate(this->TableTopEccentricRotationToPatientSupportRotationTransform);

  this->TableTopToTableEccentricRotationConcatenatedTransform->Concatenate(this->TableTopEccentricRotationToPatientSupportRotationConcatenatedTransform);
  this->TableTopToTableEccentricRotationConcatenatedTransform->Concatenate(this->TableTopToTableTopEccentricRotationTransform);

  this->PatientToTableTopConcatenatedTransform->Concatenate(this->TableTopToTableEccentricRotationConcatenatedTransform);
  this->PatientToTableTopConcatenatedTransform->Concatenate(this->PatientToTableTopTransform);

  this->DICOMToPatientConcatenatedTransform->Concatenate(this->PatientToTableTopConcatenatedTransform);
  this->DICOMToPatientConcatenatedTransform->Concatenate(this->DICOMToPatientTransform);

  this->PatientImageRegularGridToDICOMConcatenatedTransform->Concatenate(this->DICOMToPatientConcatenatedTransform);
  this->PatientImageRegularGridToDICOMConcatenatedTransform->Concatenate(this->PatientImageRegularGridToDICOMTransform);

  this->RasToPatientConcatenatedTransform->Concatenate(this->PatientToTableTopConcatenatedTransform);
  this->RasToPatientConcatenatedTransform->Concatenate(this->RasToPatientTransform);
}

//-----------------------------------------------------------------------------
vtkIECTransformLogic::~vtkIECTransformLogic()
{
  this->CoordinateSystemsMap.clear();
  this->IECTransforms.clear();
  this->ElementaryTransforms.clear();
}

//----------------------------------------------------------------------------
void vtkIECTransformLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << std::endl << "Elementary tansforms:" << std::endl;
  os << indent << "FixedReferenceToRasTransform: " << this->FixedReferenceToRasTransform << std::endl;
  os << indent << "GantryToFixedReferenceTransform: " << this->GantryToFixedReferenceTransform << std::endl;
  os << indent << "CollimatorToGantryTransform: " << this->CollimatorToGantryTransform << std::endl;
  os << indent << "WedgeFilterToCollimatorTransform: " << this->WedgeFilterToCollimatorTransform << std::endl;
  os << indent << "LeftImagingPanelToGantryTransform: " << this->LeftImagingPanelToGantryTransform << std::endl;
  os << indent << "RightImagingPanelToGantryTransform: " << this->RightImagingPanelToGantryTransform << std::endl;
  os << indent << "FlatPanelToGantryTransform: " << this->FlatPanelToGantryTransform << std::endl;
  os << indent << "PatientSupportRotationToFixedReferenceTransform: " << this->PatientSupportRotationToFixedReferenceTransform << std::endl;
  os << indent << "PatientSupportToPatientSupportRotationTransform: " << this->PatientSupportToPatientSupportRotationTransform << std::endl;
  os << indent << "TableTopEccentricRotationToPatientSupportRotationTransform: " << this->TableTopEccentricRotationToPatientSupportRotationTransform << std::endl;
  os << indent << "TableTopToTableTopEccentricRotationTransform: " << this->TableTopToTableTopEccentricRotationTransform << std::endl;
  os << indent << "PatientToTableTopTransform: " << this->PatientToTableTopTransform << std::endl;
  os << indent << "DICOMToPatientTransform: " << this->DICOMToPatientTransform << std::endl;
  os << indent << "PatientImageRegularGridToDICOMTransform: " << this->PatientImageRegularGridToDICOMTransform << std::endl;
  os << indent << "RasToPatientTransform: " << this->RasToPatientTransform << std::endl;

  os << indent << std::endl << "Concatenated transforms:" << std::endl;
  os << indent << "GantryToFixedReferenceConcatenatedTransform: " << this->GantryToFixedReferenceConcatenatedTransform << std::endl;
  os << indent << "CollimatorToGantryConcatenatedTransform: " << this->CollimatorToGantryConcatenatedTransform << std::endl;
  os << indent << "WedgeFilterToCollimatorConcatenatedTransform: " << this->WedgeFilterToCollimatorConcatenatedTransform << std::endl;
  os << indent << "LeftImagingPanelToGantryConcatenatedTransform: " << this->LeftImagingPanelToGantryConcatenatedTransform << std::endl;
  os << indent << "RightImagingPanelToGantryConcatenatedTransform: " << this->RightImagingPanelToGantryConcatenatedTransform << std::endl;
  os << indent << "FlatPanelToGantryConcatenatedTransform: " << this->FlatPanelToGantryConcatenatedTransform << std::endl;
  os << indent << "PatientSupportRotationToFixedReferenceConcatenatedTransform: " << this->PatientSupportRotationToFixedReferenceConcatenatedTransform << std::endl;
  os << indent << "PatientSupportToPatientSupportRotationConcatenatedTransform: " << this->PatientSupportToPatientSupportRotationConcatenatedTransform << std::endl;
  os << indent << "TableTopEccentricRotationToPatientSupportRotationConcatenatedTransform: " << this->TableTopEccentricRotationToPatientSupportRotationConcatenatedTransform << std::endl;
  os << indent << "TableTopToTableTopEccentricRotationConcatenatedTransform: " << this->TableTopToTableEccentricRotationConcatenatedTransform << std::endl;
  os << indent << "PatientToTableTopConcatenatedTransform: " << this->PatientToTableTopConcatenatedTransform << std::endl;
  os << indent << "DICOMToPatientConcatenatedTransform: " << this->DICOMToPatientConcatenatedTransform << std::endl;
  os << indent << "PatientImageRegularGridToDICOMTransform: " << this->PatientImageRegularGridToDICOMTransform << std::endl;
  os << indent << "RasToPatientConcatenatedTransform: " << this->RasToPatientConcatenatedTransform << std::endl;
}

//----------------------------------------------------------------------------
void vtkIECTransformLogic::UpdateGantryToFixedReferenceTransform(double gantryRotationAngleDeg, double gantryPitchAngleDeg)
{
  this->GantryToFixedReferenceTransform->Identity();
  this->GantryToFixedReferenceTransform->RotateX(gantryPitchAngleDeg);
  this->GantryToFixedReferenceTransform->RotateY(gantryRotationAngleDeg);
}

//----------------------------------------------------------------------------
void vtkIECTransformLogic::UpdateCollimatorToGantryTransform(double collimatorRotationAngleDeg, double bz)
{
  this->CollimatorToGantryTransform->Identity();
  this->CollimatorToGantryTransform->Translate(0, 0, bz);
  this->CollimatorToGantryTransform->RotateZ(collimatorRotationAngleDeg);
}

//-----------------------------------------------------------------------------
void vtkIECTransformLogic::UpdateWedgeFilterToCollimatorTransform(double wedgefilterRotationAngleDeg, double wz)
{
  this->WedgeFilterToCollimatorTransform->Identity();
  this->WedgeFilterToCollimatorTransform->Translate(0, 0, wz);
  this->WedgeFilterToCollimatorTransform->RotateZ(wedgefilterRotationAngleDeg);
}

//-----------------------------------------------------------------------------
void vtkIECTransformLogic::UpdatePatientSupportRotationToFixedReferenceTransform(double patientSupportRotationAngleDeg)
{
  this->PatientSupportRotationToFixedReferenceTransform->Identity();
  this->PatientSupportRotationToFixedReferenceTransform->RotateZ(patientSupportRotationAngleDeg);
}

//-----------------------------------------------------------------------------
void vtkIECTransformLogic::UpdateTableTopEccentricRotationToPatientSupportRotationTransform(double tableTopEccentricRotationAngleDeg, double ey)
{
  this->TableTopEccentricRotationToPatientSupportRotationTransform->Identity();
  this->TableTopEccentricRotationToPatientSupportRotationTransform->Translate(0, ey, 0);
  this->TableTopEccentricRotationToPatientSupportRotationTransform->RotateZ(tableTopEccentricRotationAngleDeg);
}

//-----------------------------------------------------------------------------
void vtkIECTransformLogic::UpdateTableTopToTableTopEccentricRotationTransform(double tx, double ty, double tz, double tableTopPitchAngleDeg, double tableTopRollAngleDeg)
{
  this->TableTopToTableTopEccentricRotationTransform->Identity();
  this->TableTopToTableTopEccentricRotationTransform->Translate(tx, ty, tz);
  this->TableTopToTableTopEccentricRotationTransform->RotateX(tableTopPitchAngleDeg);
  this->TableTopToTableTopEccentricRotationTransform->RotateY(tableTopRollAngleDeg);
}

//-----------------------------------------------------------------------------
void vtkIECTransformLogic::UpdatePatientToTableTopTransform(double px, double py, double pz, double patientPsiAngleDeg, double patientPhiAngleDeg, double patientThetaAngleDeg)
{
  this->PatientToTableTopTransform->Identity();
  this->PatientToTableTopTransform->Translate(px, py, pz);
  this->PatientToTableTopTransform->RotateX(patientPsiAngleDeg);
  this->PatientToTableTopTransform->RotateY(patientPhiAngleDeg);
  this->PatientToTableTopTransform->RotateZ(patientThetaAngleDeg);
}

//-----------------------------------------------------------------------------
void vtkIECTransformLogic::UpdatePatientImageRegularGridToDICOMTransform(double columnPixelSpacing, double rowPixelSpacing, double sliceDistance, double sx, double sy, double sz,
                                                                         double directionCosineXx, double directionCosineXy, double directionCosineXz,
                                                                         double directionCosineYx, double directionCosineYy, double directionCosineYz)
{
  this->PatientImageRegularGridToDICOMTransform->Identity();

  double directionCosineZx = directionCosineXy*directionCosineYz - directionCosineXz*directionCosineYy;
  double directionCosineZy = directionCosineXz*directionCosineYx - directionCosineXx*directionCosineYz;
  double directionCosineZz = directionCosineXx*directionCosineYy - directionCosineXy*directionCosineYx;
  double m[16] = {directionCosineXx*columnPixelSpacing, directionCosineYx*rowPixelSpacing, directionCosineZx*sliceDistance, sx,
                  directionCosineXy*columnPixelSpacing, directionCosineYy*rowPixelSpacing, directionCosineZy*sliceDistance, sy,
                  directionCosineXz*columnPixelSpacing, directionCosineYz*rowPixelSpacing, directionCosineZz*sliceDistance, sz,
                  0, 0, 0, 1};
  this->PatientImageRegularGridToDICOMTransform->Concatenate(m);
}

//-----------------------------------------------------------------------------
vtkTransform* vtkIECTransformLogic::GetElementaryTransformBetween(
  CoordinateSystemIdentifier fromFrame, CoordinateSystemIdentifier toFrame)
{
  std::string requestedTransformName = this->GetTransformNameBetween(fromFrame, toFrame);
  for (auto& transform : this->ElementaryTransforms)
  {
    std::string currentTransformName(transform->GetObjectName());
    if (currentTransformName == requestedTransformName)
    {
      return transform;
    }
  }

  vtkErrorMacro("GetElementaryTransformBetween: Elementary transform not found: " << requestedTransformName);
  return nullptr;
}

//-----------------------------------------------------------------------------
std::string vtkIECTransformLogic::GetTransformNameBetween(
  CoordinateSystemIdentifier fromFrame, CoordinateSystemIdentifier toFrame)
{
  return this->CoordinateSystemsMap[fromFrame] + "To" + this->CoordinateSystemsMap[toFrame] + "Transform";
}

//-----------------------------------------------------------------------------
bool vtkIECTransformLogic::GetTransformBetween(vtkIECTransformLogic::CoordinateSystemIdentifier fromFrame, vtkIECTransformLogic::CoordinateSystemIdentifier toFrame,
  vtkGeneralTransform* outputTransform, bool transformForBeam/*=false*/)
{
  if (!outputTransform)
  {
    vtkErrorMacro("GetTransformBetween: Invalid output transform node");
    return false;
  }

  vtkIECTransformLogic::CoordinateSystemsList fromFramePath, toFramePath;
  if (this->GetPathToRoot(fromFrame, fromFramePath) && this->GetPathFromRoot(toFrame, toFramePath))
  {
    std::vector< vtkIECTransformLogic::CoordinateSystemIdentifier > toFrameVector(toFramePath.size());
    std::vector< vtkIECTransformLogic::CoordinateSystemIdentifier > fromFrameVector(fromFramePath.size());

    std::copy(toFramePath.begin(), toFramePath.end(), toFrameVector.begin());
    std::copy(fromFramePath.begin(), fromFramePath.end(), fromFrameVector.begin());

    outputTransform->Identity();
    outputTransform->PostMultiply();
    for (size_t i = 0; i < fromFrameVector.size() - 1; ++i)
    {
      vtkIECTransformLogic::CoordinateSystemIdentifier parent, child;
      child = fromFrameVector[i];
      parent = fromFrameVector[i + 1];

      if (child == parent)
      {
        continue;
      }

      vtkTransform* fromTransform = this->GetElementaryTransformBetween(child, parent);
      if (fromTransform)
      {
        outputTransform->Concatenate(fromTransform->GetMatrix());
      }
      else
      {
        vtkErrorMacro("GetTransformBetween: Transform node is invalid");
        return false;
      }
    }

    for (size_t i = 0; i < toFrameVector.size() - 1; ++i)
    {
      vtkIECTransformLogic::CoordinateSystemIdentifier parent, child;
      parent = toFrameVector[i];
      child = toFrameVector[i + 1];

      if (child == parent)
      {
        continue;
      }

      vtkTransform* toTransform = this->GetElementaryTransformBetween(child, parent);
      if (toTransform)
      {
        vtkNew<vtkMatrix4x4> mat;
        toTransform->GetMatrix(mat);
        if (!transformForBeam) // Do not invert for beam transformation
        {
          mat->Invert();
        }
        outputTransform->Concatenate(mat);
      }
      else
      {
        vtkErrorMacro("GetTransformBetween: Transform node is invalid");
        return false;
      }
    }

    outputTransform->Modified();
    return true;
  }

  vtkErrorMacro("GetTransformBetween: Failed to get transform " << this->GetTransformNameBetween(fromFrame, toFrame));
  return false;
}

//-----------------------------------------------------------------------------
bool vtkIECTransformLogic::GetPathToRoot(vtkIECTransformLogic::CoordinateSystemIdentifier frame, vtkIECTransformLogic::CoordinateSystemsList& path)
{
  if (frame == vtkIECTransformLogic::CoordinateSystemIdentifier::FixedReference)
  {
    path.push_back(vtkIECTransformLogic::FixedReference);
    return true;
  }

  bool found = false;
  do
  {
    for (auto& pair : this->CoordinateSystemsHierarchy)
    {
      vtkIECTransformLogic::CoordinateSystemIdentifier parent = pair.first;

      auto& children = pair.second;
      auto iter = std::find(children.begin(), children.end(), frame);
      if (iter != children.end())
      {
        vtkIECTransformLogic::CoordinateSystemIdentifier id = *iter;

        vtkDebugMacro("GetPathToRoot: Checking affine transformation "
          << "\"" << this->CoordinateSystemsMap[id] << "\" -> "
          << "\"" << this->CoordinateSystemsMap[parent] << "\"");

        frame = parent;
        path.push_back(id);
        if (frame != vtkIECTransformLogic::FixedReference)
        {
          found = true;
          break;
        }
        else
        {
          path.push_back(vtkIECTransformLogic::FixedReference);
        }
      }
      else
      {
        found = false;
      }
    }
  } while (found);

  return (path.size() > 0);
}

//-----------------------------------------------------------------------------
bool vtkIECTransformLogic::GetPathFromRoot(vtkIECTransformLogic::CoordinateSystemIdentifier frame, vtkIECTransformLogic::CoordinateSystemsList& path)
{
  if (this->GetPathToRoot(frame, path))
  {
    std::reverse(path.begin(), path.end());
    return true;
  }
  else
  {
    return false;
  }
}
