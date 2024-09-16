//
//  Copyright (C) 2023 Novartis Biomedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_MOLSTANDARDIZE_PIPELINE_H
#define RD_MOLSTANDARDIZE_PIPELINE_H
#include <RDGeneral/export.h>
#include <GraphMol/RWMol.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace RDKit {

namespace MolStandardize {

struct RDKIT_MOLSTANDARDIZE_EXPORT PipelineOptions {
  // parsing
  bool strictParsing{false};

  // validation
  bool reportAllFailures{true};
  bool allowEmptyMolecules{false};
  bool allowEnhancedStereo{false};
  bool allowAromaticBondType{false};
  bool allowDativeBondType{false};
  double is2DZeroThreshold{1e-3};
  double atomClashLimit{0.03};
  double minMedianBondLength{1e-3};
  double bondLengthLimit{100.};
  bool allowLongBondsInRings{true};
  bool allowAtomBondClashExemption{true};

  // cleanup/standardization
  // metal disconnector options
  std::string metalNof{"[Li,Na,K,Rb,Cs,Fr]~[#7,#8,F]"};
  std::string metalNon{};
  // normalizer options
  std::string normalizerData{
      "// Name\tSMIRKS\n"
      "Nitro to N+(O-)=O\t[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]\n"
      "Sulfone to S(=O)(=O)\t[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])\n"
      "Pyridine oxide to n+O-\t[nH0+0:1]=[OH0+0:2]>>[n+:1][O-:2]\n"
      "Azide to N=N+=N-\t[*:1][N:2]=[N:3]#[N:4]>>[*:1][N:2]=[N+:3]=[N-:4]\n"
      "Diazo/azo to =N+=N-\t[*:1]=[N:2]#[N:3]>>[*:1]=[N+:2]=[N-:3]\n"
      // Note: the sulfoxide transformation by default included in the
      // Normalizer configuration was removed Note: the transformation below was
      // ported from STRUCHK and it's not part of the default Normalizer
      // configuration
      "[SH](=O)(=O) to S(=O)O\t[c,C,N,O,F,Cl,Br,I:1][SH+0:2](=[O:3])=[O:4]>>[*:1][*:2]([*:3])=[*:4]\n"
      // Note: the two transformations below replace the default Phosphate
      // normalization in order to ensure that, if an O is available, the double
      // bond is placed between P and O
      "Phosphate to P(O-)=O\t[O-:1][P+;D4:2][O,S,Se,Te;-1:3]>>[O+0:1]=[P+0;D5:2][*-1:3]\n"
      "Generalized phosphate to P(X-)=Y\t[S,Se,Te;-1:1][P+;D4:2][S,Se,Te;-1:3]>>[*+0:1]=[P+0;D5:2][*-1:3]\n"
      "C/S+N to C/S=N+\t[C,S&!$([S+]-[O-]);X3+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]\n"
      "P+N to P=N+\t[P;X4+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]\n"
      "Recombine 1,3-separated charges\t[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[N,P,As,Sb,O,S,Se,Te;+1:3]>>[*-0:1]=[*:2]-[*+0:3]\n"
      "Recombine 1,3-separated charges\t[n,o,p,s;-1:1]:[a:2]=[N,O,P,S;+1:3]>>[*-0:1]:[*:2]-[*+0:3]\n"
      "Recombine 1,3-separated charges\t[N,O,P,S;-1:1]-[a+0:2]:[n,o,p,s;+1:3]>>[*-0:1]=[*:2]:[*+0:3]\n"
      "Recombine 1,5-separated charges\t[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[A:3]-[A:4]=[N,P,As,Sb,O,S,Se,Te;+1:5]>>[*-0:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n"
      "Recombine 1,5-separated charges\t[n,o,p,s;-1:1]:[a:2]:[a:3]:[c:4]=[N,O,P,S;+1:5]>>[*-0:1]:[*:2]:[*:3]:[c:4]-[*+0:5]\n"
      "Recombine 1,5-separated charges\t[N,O,P,S;-1:1]-[c:2]:[a:3]:[a:4]:[n,o,p,s;+1:5]>>[*-0:1]=[c:2]:[*:3]:[*:4]:[*+0:5]\n"
      // Note: four transformations were added to the normalization of aliphatic
      // conjug cations in order to favor the positioning of new double bonds
      // within rings
      "Normalize 1,3 conjugated cation\t[N;+0!H0:1]@-[A:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]\n"
      "Normalize 1,5 conjugated cation\t[N;+0!H0:1]@-[A:2]=[A:3]@-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n"
      "Normalize 1,3 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]\n"
      "Normalize 1,3 conjugated cation\t[n;+0!H0:1]:[c:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]:[*:2]-[*+0:3]\n"
      "Normalize 1,5 conjugated cation\t[N;+0!H0:1]@-[A:2]=[A:3]-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n"
      "Normalize 1,5 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[A:3]@-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n"
      "Normalize 1,5 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[A:3]-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\n"
      "Normalize 1,5 conjugated cation\t[n;+0!H0:1]:[a:2]:[a:3]:[c:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]-[*+0:5]\n"
      "Charge normalization\t[F,Cl,Br,I,At;-1:1]=[O:2]>>[*-0:1][O-:2]\n"
      "Charge recombination\t[N,P,As,Sb;-1:1]=[C+;v3:2]>>[*+0:1]#[C+0:2]\n"};
  unsigned int normalizerMaxRestarts{200};
  double scaledMedianBondLength{1.};

  // serialization
  bool outputV2000{false};
};

enum RDKIT_MOLSTANDARDIZE_EXPORT PipelineStatus {
  NO_EVENT = 0,
  INPUT_ERROR = (1 << 0),
  PREPARE_FOR_VALIDATION_ERROR = (1 << 1),
  FEATURES_VALIDATION_ERROR = (1 << 2),
  BASIC_VALIDATION_ERROR = (1 << 3),
  IS2D_VALIDATION_ERROR = (1 << 4),
  LAYOUT2D_VALIDATION_ERROR = (1 << 5),
  STEREO_VALIDATION_ERROR = (1 << 6),
  VALIDATION_ERROR = (FEATURES_VALIDATION_ERROR | BASIC_VALIDATION_ERROR |
                      IS2D_VALIDATION_ERROR | LAYOUT2D_VALIDATION_ERROR |
                      STEREO_VALIDATION_ERROR),
  PREPARE_FOR_STANDARDIZATION_ERROR = (1 << 7),
  METAL_STANDARDIZATION_ERROR = (1 << 8),
  NORMALIZER_STANDARDIZATION_ERROR = (1 << 9),
  FRAGMENT_STANDARDIZATION_ERROR = (1 << 10),
  CHARGE_STANDARDIZATION_ERROR = (1 << 11),
  STANDARDIZATION_ERROR =
      (METAL_STANDARDIZATION_ERROR | NORMALIZER_STANDARDIZATION_ERROR |
       FRAGMENT_STANDARDIZATION_ERROR | CHARGE_STANDARDIZATION_ERROR),
  OUTPUT_ERROR = (1 << 12),
  PIPELINE_ERROR = (INPUT_ERROR | PREPARE_FOR_VALIDATION_ERROR |
                    VALIDATION_ERROR | PREPARE_FOR_STANDARDIZATION_ERROR |
                    STANDARDIZATION_ERROR | OUTPUT_ERROR),
  METALS_DISCONNECTED = (1 << 23),
  NORMALIZATION_APPLIED = (1 << 24),
  FRAGMENTS_REMOVED = (1 << 25),
  PROTONATION_CHANGED = (1 << 26),
  STRUCTURE_MODIFICATION = (METALS_DISCONNECTED | NORMALIZATION_APPLIED |
                            FRAGMENTS_REMOVED | PROTONATION_CHANGED)
};

enum class RDKIT_MOLSTANDARDIZE_EXPORT PipelineStage : std::uint32_t {
  NOT_STARTED = 0,
  PARSING_INPUT,
  PREPARE_FOR_VALIDATION,
  VALIDATION,
  PREPARE_FOR_STANDARDIZATION,
  STANDARDIZATION,
  REAPPLY_WEDGING,
  CLEANUP_2D,
  MAKE_PARENT,
  SERIALIZING_OUTPUT,
  COMPLETED
};

struct RDKIT_MOLSTANDARDIZE_EXPORT PipelineLogEntry {
  PipelineStatus status;
  std::string detail;
};

using PipelineLog = std::vector<PipelineLogEntry>;

struct RDKIT_MOLSTANDARDIZE_EXPORT PipelineResult {
  PipelineStatus status;
  std::uint32_t stage;
  PipelineLog log;
  std::string inputMolData;
  std::string outputMolData;
  std::string parentMolData;

  void append(PipelineStatus newStatus, const std::string &info);
};

using RWMOL_SPTR_PAIR = std::pair<RWMOL_SPTR, RWMOL_SPTR>;

namespace Operations {
RDKIT_MOLSTANDARDIZE_EXPORT RWMOL_SPTR prepareForValidation(
    RWMOL_SPTR mol, PipelineResult &result, const PipelineOptions &options);
RDKIT_MOLSTANDARDIZE_EXPORT RWMOL_SPTR validate(RWMOL_SPTR mol,
                                                PipelineResult &result,
                                                const PipelineOptions &options);
RDKIT_MOLSTANDARDIZE_EXPORT RWMOL_SPTR prepareForStandardization(
    RWMOL_SPTR mol, PipelineResult &result, const PipelineOptions &options);
RDKIT_MOLSTANDARDIZE_EXPORT RWMOL_SPTR standardize(
    RWMOL_SPTR mol, PipelineResult &result, const PipelineOptions &options);
RDKIT_MOLSTANDARDIZE_EXPORT RWMOL_SPTR reapplyWedging(
    RWMOL_SPTR mol, PipelineResult &result, const PipelineOptions &options);
RDKIT_MOLSTANDARDIZE_EXPORT RWMOL_SPTR cleanup2D(
    RWMOL_SPTR mol, PipelineResult &result, const PipelineOptions &options);
RDKIT_MOLSTANDARDIZE_EXPORT RWMOL_SPTR_PAIR makeParent(
    RWMOL_SPTR mol, PipelineResult &result, const PipelineOptions &options);

RDKIT_MOLSTANDARDIZE_EXPORT RWMOL_SPTR parse(const std::string &molblock,
                                             PipelineResult &result,
                                             const PipelineOptions &options);
RDKIT_MOLSTANDARDIZE_EXPORT void serialize(RWMOL_SPTR_PAIR output,
                                           PipelineResult &result,
                                           const PipelineOptions &options);

using ParseOperation = decltype(&parse);
using SerializeOperation = decltype(&serialize);
using Operation = decltype(&prepareForValidation);
using ParentOperation = decltype(&makeParent);
using PipelineVector = std::vector<std::pair<std::uint32_t, Operation>>;

const PipelineVector validationSteps{
    // input sanitization and cleanup
    {static_cast<uint32_t>(PipelineStage::PREPARE_FOR_VALIDATION),
     &prepareForValidation},
    //  validate the structure
    {static_cast<uint32_t>(PipelineStage::VALIDATION), &validate}};

const PipelineVector standardizationSteps{
    {static_cast<uint32_t>(PipelineStage::PREPARE_FOR_STANDARDIZATION),
     &prepareForStandardization},
    {static_cast<uint32_t>(PipelineStage::STANDARDIZATION), &standardize},
    {static_cast<uint32_t>(PipelineStage::REAPPLY_WEDGING), &reapplyWedging},
    {static_cast<uint32_t>(PipelineStage::CLEANUP_2D), &cleanup2D}};
}  // namespace Operations

class RDKIT_MOLSTANDARDIZE_EXPORT Pipeline {
 private:
  PipelineOptions options;
  Operations::ParseOperation parse = Operations::parse;
  Operations::SerializeOperation serialize = Operations::serialize;
  Operations::PipelineVector validationSteps = Operations::validationSteps;
  Operations::PipelineVector standardizationSteps =
      Operations::standardizationSteps;
  Operations::ParentOperation makeParent = Operations::makeParent;

 public:
  Pipeline() = default;
  explicit Pipeline(const PipelineOptions &o) : options(o){};
  ~Pipeline() = default;

  PipelineResult run(const std::string &molblock) const;

  void setValidationSteps(const Operations::PipelineVector &steps) {
    validationSteps = steps;
  }
  void setStandardizationSteps(const Operations::PipelineVector &steps) {
    standardizationSteps = steps;
  }
  void setMakeParent(Operations::ParentOperation op) { makeParent = op; }
  void setParse(Operations::ParseOperation op) { parse = op; }
  void setSerialize(Operations::SerializeOperation op) { serialize = op; }

 private:
};

}  // namespace MolStandardize
}  // namespace RDKit

#endif
