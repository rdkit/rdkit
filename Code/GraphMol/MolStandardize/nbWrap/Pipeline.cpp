//
//  Copyright (C) 2023-2026 Novartis Biomedical Research and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/bind_vector.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Pipeline.h>

namespace RDKit {
namespace MolStandardize {

bool operator==(const PipelineLogEntry &lhs, const PipelineLogEntry &rhs) {
  return (lhs.status == rhs.status) && (lhs.detail == rhs.detail);
}

}  // namespace MolStandardize
}  // namespace RDKit

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

void wrap_pipeline(nb::module_ &m) {
  nb::class_<MolStandardize::PipelineOptions>(m, "PipelineOptions")
      .def(nb::init<>())
      .def_rw("strictParsing",
               &MolStandardize::PipelineOptions::strictParsing)
      .def_rw("reportAllFailures",
               &MolStandardize::PipelineOptions::reportAllFailures)
      .def_rw("allowEmptyMolecules",
               &MolStandardize::PipelineOptions::allowEmptyMolecules)
      .def_rw("allowEnhancedStereo",
               &MolStandardize::PipelineOptions::allowEnhancedStereo)
      .def_rw("allowAromaticBondType",
               &MolStandardize::PipelineOptions::allowAromaticBondType)
      .def_rw("allowDativeBondType",
               &MolStandardize::PipelineOptions::allowDativeBondType)
      .def_rw("is2DZeroThreshold",
               &MolStandardize::PipelineOptions::is2DZeroThreshold)
      .def_rw("atomClashLimit",
               &MolStandardize::PipelineOptions::atomClashLimit)
      .def_rw("minMedianBondLength",
               &MolStandardize::PipelineOptions::minMedianBondLength)
      .def_rw("bondLengthLimit",
               &MolStandardize::PipelineOptions::bondLengthLimit)
      .def_rw("allowLongBondsInRings",
               &MolStandardize::PipelineOptions::allowLongBondsInRings)
      .def_rw("allowAtomBondClashExemption",
               &MolStandardize::PipelineOptions::allowAtomBondClashExemption)
      .def_rw("metalNof", &MolStandardize::PipelineOptions::metalNof)
      .def_rw("metalNon", &MolStandardize::PipelineOptions::metalNon)
      .def_rw("normalizerData",
               &MolStandardize::PipelineOptions::normalizerData)
      .def_rw("normalizerMaxRestarts",
               &MolStandardize::PipelineOptions::normalizerMaxRestarts)
      .def_rw("scaledMedianBondLength",
               &MolStandardize::PipelineOptions::scaledMedianBondLength)
      .def_rw("outputV2000", &MolStandardize::PipelineOptions::outputV2000);

  nb::enum_<MolStandardize::PipelineStatus>(m, "PipelineStatus",
                                             nb::is_arithmetic(), nb::is_flag())
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
      .value(
          "NORMALIZER_STANDARDIZATION_ERROR",
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
             MolStandardize::PipelineStatus::STRUCTURE_MODIFICATION)
;

  nb::enum_<MolStandardize::PipelineStage>(m, "PipelineStage")
      .value("PARSING_INPUT", MolStandardize::PipelineStage::PARSING_INPUT)
      .value("PREPARE_FOR_VALIDATION",
             MolStandardize::PipelineStage::PREPARE_FOR_VALIDATION)
      .value("VALIDATION", MolStandardize::PipelineStage::VALIDATION)
      .value("PREPARE_FOR_STANDARDIZATION",
             MolStandardize::PipelineStage::PREPARE_FOR_STANDARDIZATION)
      .value("STANDARDIZATION",
             MolStandardize::PipelineStage::STANDARDIZATION)
      .value("SERIALIZING_OUTPUT",
             MolStandardize::PipelineStage::SERIALIZING_OUTPUT)
      .value("COMPLETED", MolStandardize::PipelineStage::COMPLETED);

  nb::class_<MolStandardize::PipelineLogEntry>(m, "PipelineLogEntry")
      .def_ro("status", &MolStandardize::PipelineLogEntry::status)
      .def_ro("detail", &MolStandardize::PipelineLogEntry::detail);

  nb::bind_vector<MolStandardize::PipelineLog>(m, "PipelineLog");

  nb::class_<MolStandardize::PipelineResult>(m, "PipelineResult")
      .def_ro("status", &MolStandardize::PipelineResult::status)
      .def_prop_ro("stage",
                   [](const MolStandardize::PipelineResult &self) {
                     return static_cast<MolStandardize::PipelineStage>(
                         self.stage);
                   })
      .def_ro("log", &MolStandardize::PipelineResult::log)
      .def_ro("inputMolData", &MolStandardize::PipelineResult::inputMolData)
      .def_ro("outputMolData",
               &MolStandardize::PipelineResult::outputMolData)
      .def_ro("parentMolData",
               &MolStandardize::PipelineResult::parentMolData);

  nb::class_<MolStandardize::Pipeline>(m, "Pipeline")
      .def(nb::init<>())
      .def(nb::init<const MolStandardize::PipelineOptions &>(), "options"_a)
      .def("run", &MolStandardize::Pipeline::run, "molData"_a);
}
