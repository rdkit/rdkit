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

bool operator==(const PipelineLogEntry & lhs, const PipelineLogEntry & rhs)
{
  return (lhs.status == rhs.status) && (lhs.detail == rhs.detail);
}

}
}


namespace python = boost::python;
using namespace RDKit;

void wrap_pipeline() {
    python::class_<MolStandardize::PipelineOptions>("PipelineOptions")
        .def_readwrite("strictParsing", &MolStandardize::PipelineOptions::strictParsing)
        .def_readwrite("reportAllFailures", &MolStandardize::PipelineOptions::reportAllFailures)
        .def_readwrite("allowEmptyMolecules", &MolStandardize::PipelineOptions::allowEmptyMolecules)
        .def_readwrite("allowEnhancedStereo", &MolStandardize::PipelineOptions::allowEnhancedStereo)
        .def_readwrite("allowAromaticBondType", &MolStandardize::PipelineOptions::allowAromaticBondType)
        .def_readwrite("is2DZeroThreshold", &MolStandardize::PipelineOptions::is2DZeroThreshold)
        .def_readwrite("atomClashLimit", &MolStandardize::PipelineOptions::atomClashLimit)
        .def_readwrite("minMedianBondLength", &MolStandardize::PipelineOptions::minMedianBondLength)
        .def_readwrite("bondLengthLimit", &MolStandardize::PipelineOptions::bondLengthLimit)
        .def_readwrite("allowLongBondsInRings", &MolStandardize::PipelineOptions::allowLongBondsInRings)
        .def_readwrite("allowAtomBondClashExemption", &MolStandardize::PipelineOptions::allowAtomBondClashExemption)
        .def_readwrite("metalNof", &MolStandardize::PipelineOptions::metalNof)
        .def_readwrite("metalNon", &MolStandardize::PipelineOptions::metalNon)
        .def_readwrite("normalizerData", &MolStandardize::PipelineOptions::normalizerData)
        .def_readwrite("normalizerMaxRestarts", &MolStandardize::PipelineOptions::normalizerMaxRestarts)
        .def_readwrite("scaledMedianBondLength", &MolStandardize::PipelineOptions::scaledMedianBondLength)
        .def_readwrite("outputV2000", &MolStandardize::PipelineOptions::outputV2000)
        ;

    python::enum_<MolStandardize::PipelineStatus>("PipelineStatus")
        .value("NO_EVENT", MolStandardize::NO_EVENT)
        .value("INPUT_ERROR", MolStandardize::INPUT_ERROR)
        .value("PREPARE_FOR_VALIDATION_ERROR", MolStandardize::PREPARE_FOR_VALIDATION_ERROR)
        .value("FEATURES_VALIDATION_ERROR", MolStandardize::FEATURES_VALIDATION_ERROR)
        .value("BASIC_VALIDATION_ERROR", MolStandardize::BASIC_VALIDATION_ERROR)
        .value("IS2D_VALIDATION_ERROR", MolStandardize::IS2D_VALIDATION_ERROR)
        .value("LAYOUT2D_VALIDATION_ERROR", MolStandardize::LAYOUT2D_VALIDATION_ERROR)
        .value("STEREO_VALIDATION_ERROR", MolStandardize::STEREO_VALIDATION_ERROR)
        .value("VALIDATION_ERROR", MolStandardize::VALIDATION_ERROR)
        .value("PREPARE_FOR_STANDARDIZATION_ERROR", MolStandardize::PREPARE_FOR_STANDARDIZATION_ERROR)
        .value("METAL_STANDARDIZATION_ERROR", MolStandardize::METAL_STANDARDIZATION_ERROR)
        .value("NORMALIZER_STANDARDIZATION_ERROR", MolStandardize::NORMALIZER_STANDARDIZATION_ERROR)
        .value("FRAGMENT_STANDARDIZATION_ERROR", MolStandardize::FRAGMENT_STANDARDIZATION_ERROR)
        .value("CHARGE_STANDARDIZATION_ERROR", MolStandardize::CHARGE_STANDARDIZATION_ERROR)
        .value("STANDARDIZATION_ERROR", MolStandardize::STANDARDIZATION_ERROR)
        .value("OUTPUT_ERROR", MolStandardize::OUTPUT_ERROR)
        .value("PIPELINE_ERROR", MolStandardize::PIPELINE_ERROR)
        .value("METALS_DISCONNECTED", MolStandardize::METALS_DISCONNECTED)
        .value("NORMALIZATION_APPLIED", MolStandardize::NORMALIZATION_APPLIED)
        .value("FRAGMENTS_REMOVED", MolStandardize::FRAGMENTS_REMOVED)
        .value("PROTONATION_CHANGED", MolStandardize::PROTONATION_CHANGED)
        .value("STRUCTURE_MODIFICATION", MolStandardize::STRUCTURE_MODIFICATION)
        ;

    python::enum_<MolStandardize::PipelineStage>("PipelineStage")
        .value("PARSING_INPUT", MolStandardize::PARSING_INPUT)
        .value("PREPARE_FOR_VALIDATION", MolStandardize::PREPARE_FOR_VALIDATION)
        .value("VALIDATION", MolStandardize::VALIDATION)
        .value("PREPARE_FOR_STANDARDIZATION", MolStandardize::PREPARE_FOR_STANDARDIZATION)
        .value("STANDARDIZATION", MolStandardize::STANDARDIZATION)
        .value("SERIALIZING_OUTPUT", MolStandardize::SERIALIZING_OUTPUT)
        .value("COMPLETED", MolStandardize::COMPLETED)
        ;

    python::class_<MolStandardize::PipelineLogEntry>("PipelineLogEntry", python::no_init)
        .def_readonly("status", &MolStandardize::PipelineLogEntry::status)
        .def_readonly("detail", &MolStandardize::PipelineLogEntry::detail)
        ;

    python::class_<MolStandardize::PipelineLog>("PipelineLog", python::no_init)
        .def(python::vector_indexing_suite<MolStandardize::PipelineLog>())
        ;

    python::class_<MolStandardize::PipelineResult>("PipelineResult", python::no_init)
        .def_readonly("status", &MolStandardize::PipelineResult::status)
        .def_readonly("stage", &MolStandardize::PipelineResult::stage)
        .def_readonly("log", &MolStandardize::PipelineResult::log)
        .def_readonly("inputMolBlock", &MolStandardize::PipelineResult::inputMolBlock)
        .def_readonly("outputMolBlock", &MolStandardize::PipelineResult::outputMolBlock)
        .def_readonly("parentMolBlock", &MolStandardize::PipelineResult::parentMolBlock)
        ;

    python::class_<MolStandardize::Pipeline>("Pipeline")
        .def(python::init<const MolStandardize::PipelineOptions &>())
        .def("run", &MolStandardize::Pipeline::run)
        ;
}
