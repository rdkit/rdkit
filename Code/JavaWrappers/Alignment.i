/*
 *
 *  Copyright (c) 2025, Greg Landrum and T5 Informatics GmbH
 *  All rights reserved.
 *
 *  This file is part of the RDKit.
 *  The contents are covered by the terms of the BSD license
 *  which is included in the file license.txt, found at the root
 *  of the RDKit source tree.
 *
 */

%{
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolAlign/O3AAlignMolecules.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <RDGeneral/Exceptions.h>
%}

// internals that can't be used from the wrappers (function pointers,
// void * property blocks, boost containers)
%ignore RDKit::MolAlign::MolAlignException;
%ignore RDKit::MolAlign::details::symmetrizeTerminalAtoms;
%ignore RDKit::MolAlign::O3AFuncData;
%ignore RDKit::MolAlign::isDoubleZero;
%ignore RDKit::MolAlign::O3AConstraint;
%ignore RDKit::MolAlign::O3AConstraintVect;
%ignore RDKit::MolAlign::MolHistogram;
%ignore RDKit::MolAlign::LAP;
%ignore RDKit::MolAlign::SDM;
// the default arguments are needed so that all generated overloads are ignored
%ignore RDKit::MolAlign::O3A::O3A(
    ROMol &prbMol, const ROMol &refMol, void *prbProp, void *refProp,
    AtomTypeScheme atomTypes = MMFF94, const int prbCid = -1,
    const int refCid = -1, const bool reflect = false,
    const unsigned int maxIters = 50, unsigned int options = 0,
    const MatchVectType *constraintMap = nullptr,
    const RDNumeric::DoubleVector *constraintWeights = nullptr,
    LAP *extLAP = nullptr, MolHistogram *extPrbHist = nullptr,
    MolHistogram *extRefHist = nullptr);
%ignore RDKit::MolAlign::O3A::O3A(
    int (*costFunc)(const unsigned int, const unsigned int, double, void *),
    double (*weightFunc)(const unsigned int, const unsigned int, void *),
    double (*scoringFunc)(const unsigned int, const unsigned int, void *),
    void *data, ROMol &prbMol, const ROMol &refMol, const int prbCid,
    const int refCid, const boost::dynamic_bitset<> &prbHvyAtoms,
    const boost::dynamic_bitset<> &refHvyAtoms, const bool reflect = false,
    const unsigned int maxIters = 50, unsigned int options = 0,
    O3AConstraintVect *o3aConstraintVect = nullptr,
    ROMol *extWorkPrbMol = nullptr, LAP *extLAP = nullptr,
    MolHistogram *extPrbHist = nullptr, MolHistogram *extRefHist = nullptr);
%ignore RDKit::MolAlign::reflect;
%ignore RDKit::MolAlign::o3aMMFFCostFunc;
%ignore RDKit::MolAlign::o3aMMFFWeightFunc;
%ignore RDKit::MolAlign::o3aMMFFScoringFunc;
%ignore RDKit::MolAlign::o3aCrippenCostFunc;
%ignore RDKit::MolAlign::o3aCrippenWeightFunc;
%ignore RDKit::MolAlign::o3aCrippenScoringFunc;
%ignore RDKit::MolAlign::getO3AForProbeConfs;
// these return pointers into the O3A; copies are returned below instead
%ignore RDKit::MolAlign::O3A::matches;
%ignore RDKit::MolAlign::O3A::weights;
%rename(matches) RDKit::MolAlign::O3A::matchesCopy;
%rename(weights) RDKit::MolAlign::O3A::weightsCopy;

// O3A keeps raw pointers to the probe and reference molecules, which are
// used again by align() and trans(). Hold references to their proxies so
// they can't be garbage collected while the O3A is still alive. prbMol and
// refMol are the argument names of the constructor defined below.
#ifdef SWIGJAVA
%typemap(javacode) RDKit::MolAlign::O3A %{
  private ROMol prbMolRef;
  private ROMol refMolRef;
%}
%typemap(javaconstruct) RDKit::MolAlign::O3A {
    this($imcall, true);
    prbMolRef = prbMol;
    refMolRef = refMol;
  }
#endif
#ifdef SWIGCSHARP
%typemap(cscode) RDKit::MolAlign::O3A %{
  private ROMol prbMolRef;
  private ROMol refMolRef;
%}
%typemap(csconstruct, excode=SWIGEXCODE) RDKit::MolAlign::O3A %{: this($imcall, true) {$excode
    prbMolRef = prbMol;
    refMolRef = refMol;
  }
%}
#endif

%include <GraphMol/MolAlign/AlignMolecules.h>
%include <GraphMol/MolAlign/O3AAlignMolecules.h>

// The C++ constructor takes the MMFF properties (or Crippen logP
// contributions) as void pointers; compute them here instead.
%extend RDKit::MolAlign::O3A {
  O3A(RDKit::ROMol &prbMol, RDKit::ROMol &refMol,
      RDKit::MolAlign::O3A::AtomTypeScheme atomTypes =
          RDKit::MolAlign::O3A::MMFF94,
      int prbCid = -1, int refCid = -1, bool reflect = false,
      unsigned int maxIters = 50, unsigned int options = 0,
      const RDKit::MatchVectType *constraintMap = nullptr,
      const RDNumeric::DoubleVector *constraintWeights = nullptr) {
    if (atomTypes == RDKit::MolAlign::O3A::CRIPPEN) {
      std::vector<double> prbLogp(prbMol.getNumAtoms());
      std::vector<double> prbMR(prbMol.getNumAtoms());
      RDKit::Descriptors::getCrippenAtomContribs(prbMol, prbLogp, prbMR, true);
      std::vector<double> refLogp(refMol.getNumAtoms());
      std::vector<double> refMR(refMol.getNumAtoms());
      RDKit::Descriptors::getCrippenAtomContribs(refMol, refLogp, refMR, true);
      return new RDKit::MolAlign::O3A(
          prbMol, refMol, &prbLogp, &refLogp, atomTypes, prbCid, refCid,
          reflect, maxIters, options, constraintMap, constraintWeights);
    }
    RDKit::MMFF::MMFFMolProperties prbMP(prbMol);
    RDKit::MMFF::MMFFMolProperties refMP(refMol);
    if (!prbMP.isValid() || !refMP.isValid()) {
      throw ValueErrorException(
          "missing MMFF94 parameters for probe or reference molecule");
    }
    return new RDKit::MolAlign::O3A(prbMol, refMol, &prbMP, &refMP, atomTypes,
                                    prbCid, refCid, reflect, maxIters, options,
                                    constraintMap, constraintWeights);
  }

  RDKit::MatchVectType matchesCopy() {
    const auto matchVect = $self->matches();
    return matchVect ? *matchVect : RDKit::MatchVectType();
  }
  RDNumeric::DoubleVector weightsCopy() {
    const auto weights = $self->weights();
    return weights ? *weights : RDNumeric::DoubleVector(0);
  }
}
