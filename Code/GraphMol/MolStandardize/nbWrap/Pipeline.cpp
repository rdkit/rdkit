//
//  Copyright (C) 2023 Novartis Biomedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Pipeline.h>

namespace RDKit {
namespace MolStandardize {

bool operator==(const PipelineLogEntry &lhs, const PipelineLogEntry &rhs) {
  return (lhs.status == rhs.status) && (lhs.detail == rhs.detail);
}

}  // namespace MolStandardize
}  // namespace RDKit

namespace python = boost::python;
using namespace RDKit;

void wrap_pipeline() {
  python::class_<MolStandardize::PipelineOptions>("PipelineOptions")
      .def_readwrite("strictParsing",
                     &MolStandardize::PipelineOptions::strictParsing)
      .def_readwrite("reportAllFailures",
                     &MolStandardize::PipelineOptions::reportAllFailures)
      .def_readwrite("allowEmptyMolecules",
                     &MolStandardize::PipelineOptions::allowEmptyMolecules)
      .def_readwrite("allowEnhancedStereo",
                     &MolStandardize::PipelineOptions::allowEnhancedStereo)
      .def_readwrite("allowAromaticBondType",
                     &MolStandardize::PipelineOptions::allowAromaticBondType)
      .def_readwrite("allowDativeBondType",
                     &MolStandardize::PipelineOptions::allowDativeBondType)
      .def_readwrite("is2DZeroThreshold",
                     &MolStandardize::PipelineOptions::is2DZeroThreshold)
      .def_readwrite("atomClashLimit",
                     &MolStandardize::PipelineOptions::atomClashLimit)
      .def_readwrite("minMedianBondLength",
                     &MolStandardize::PipelineOptions::minMedianBondLength)
      .def_readwrite("bondLengthLimit",
                     &MolStandardize::PipelineOptions::bondLengthLimit)
      .def_readwrite("allowLongBondsInRings",
                     &MolStandardize::PipelineOptions::allowLongBondsInRings)
      .def_readwrite(
          "allowAtomBondClashExemption",
          &MolStandardize::PipelineOptions::allowAtomBondClashExemption)
      .def_readwrite("metalNof", &MolStandardize::PipelineOptions::metalNof)
      .def_readwrite("metalNon", &MolStandardize::PipelineOptions::metalNon)
      .def_readwrite("normalizerData",
                     &MolStandardize::PipelineOptions::normalizerData)
      .def_readwrite("normalizerMaxRestarts",
                     &MolStandardize::PipelineOptions::normalizerMaxRestarts)
      .def_readwrite("scaledMedianBondLength",
                     &MolStandardize::PipelineOptions::scaledMedianBondLength)
      .def_readwrite("outputV2000",
                     &MolStandardize::PipelineOptions::outputV2000);

  python::enum_<MolStandardize::PipelineStatus>("PipelineStatus")
      .value("NO_EVENT", MolStandardize::PipelineStatus::NO_EVENT)
      .value("INPUT_ERROR", MolStandardize::PipelineStatus::INPUT_ERROR)
      .value("PREPARE_FOR_VALIDATION_ERROR",
             MolStandardize::PipelineStatus::PREPARE_FOR_VALIDATION_ERROR)
      .value("FEATURES_VALIDATION_ERROR",
             MolStandardize::PipelineStatus::FEATURES_VALIDATION_ERROR)
      .value("BASIC_VALIDATION_ERROR",
             MolStandardize::PipelineStatus::BASIC_VALIDATION_ERROR)
      .value("IS2D_VALIDATION_ERROR",
             MolStandardize::PipelineStatus::IS2D_VALIDATION_ERROR)
      .value("LAYOUT2D_VALIDATION_ERROR",
             MolStandardize::PipelineStatus::LAYOUT2D_VALIDATION_ERROR)
      .value("STEREO_VALIDATION_ERROR",
             MolStandardize::PipelineStatus::STEREO_VALIDATION_ERROR)
      .value("VALIDATION_ERROR",
             MolStandardize::PipelineStatus::VALIDATION_ERROR)
      .value("PREPARE_FOR_STANDARDIZATION_ERROR",
             MolStandardize::PipelineStatus::PREPARE_FOR_STANDARDIZATION_ERROR)
      .value("METAL_STANDARDIZATION_ERROR",
             MolStandardize::PipelineStatus::METAL_STANDARDIZATION_ERROR)
      .value("NORMALIZER_STANDARDIZATION_ERROR",
             MolStandardize::PipelineStatus::NORMALIZER_STANDARDIZATION_ERROR)
      .value("FRAGMENT_STANDARDIZATION_ERROR",
             MolStandardize::PipelineStatus::FRAGMENT_STANDARDIZATION_ERROR)
      .value("CHARGE_STANDARDIZATION_ERROR",
             MolStandardize::PipelineStatus::CHARGE_STANDARDIZATION_ERROR)
      .value("STANDARDIZATION_ERROR",
             MolStandardize::PipelineStatus::STANDARDIZATION_ERROR)
      .value("OUTPUT_ERROR", MolStandardize::PipelineStatus::OUTPUT_ERROR)
      .value("PIPELINE_ERROR", MolStandardize::PipelineStatus::PIPELINE_ERROR)
      .value("METALS_DISCONNECTED",
             MolStandardize::PipelineStatus::METALS_DISCONNECTED)
      .value("NORMALIZATION_APPLIED",
             MolStandardize::PipelineStatus::NORMALIZATION_APPLIED)
      .value("FRAGMENTS_REMOVED",
             MolStandardize::PipelineStatus::FRAGMENTS_REMOVED)
      .value("PROTONATION_CHANGED",
             MolStandardize::PipelineStatus::PROTONATION_CHANGED)
      .value("STRUCTURE_MODIFICATION",
             MolStandardize::PipelineStatus::STRUCTURE_MODIFICATION);

  python::enum_<MolStandardize::PipelineStage>("PipelineStage")
      .value("PARSING_INPUT", MolStandardize::PipelineStage::PARSING_INPUT)
      .value("PREPARE_FOR_VALIDATION",
             MolStandardize::PipelineStage::PREPARE_FOR_VALIDATION)
      .value("VALIDATION", MolStandardize::PipelineStage::VALIDATION)
      .value("PREPARE_FOR_STANDARDIZATION",
             MolStandardize::PipelineStage::PREPARE_FOR_STANDARDIZATION)
      .value("STANDARDIZATION", MolStandardize::PipelineStage::STANDARDIZATION)
      .value("SERIALIZING_OUTPUT",
             MolStandardize::PipelineStage::SERIALIZING_OUTPUT)
      .value("COMPLETED", MolStandardize::PipelineStage::COMPLETED);

  python::class_<MolStandardize::PipelineLogEntry>("PipelineLogEntry",
                                                   python::no_init)
      .def_readonly("status", &MolStandardize::PipelineLogEntry::status)
      .def_readonly("detail", &MolStandardize::PipelineLogEntry::detail);

  python::class_<MolStandardize::PipelineLog>("PipelineLog", python::no_init)
      .def(python::vector_indexing_suite<MolStandardize::PipelineLog>());

  python::class_<MolStandardize::PipelineResult>("PipelineResult",
                                                 python::no_init)
      .def_readonly("status", &MolStandardize::PipelineResult::status)
      .def_readonly("stage", &MolStandardize::PipelineResult::stage)
      .def_readonly("log", &MolStandardize::PipelineResult::log)
      .def_readonly("inputMolData",
                    &MolStandardize::PipelineResult::inputMolData)
      .def_readonly("outputMolData",
                    &MolStandardize::PipelineResult::outputMolData)
      .def_readonly("parentMolData",
                    &MolStandardize::PipelineResult::parentMolData);

  python::class_<MolStandardize::Pipeline>("Pipeline")
      .def(python::init<const MolStandardize::PipelineOptions &>())
      .def("run", &MolStandardize::Pipeline::run);
}
