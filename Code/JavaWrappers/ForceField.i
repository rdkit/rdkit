/* 
* $Id$
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
%{
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/Contrib.h>
#include <ForceField/UFF/Params.h>
#include <ForceField/UFF/AngleBend.h>
#include <ForceField/UFF/BondStretch.h>
#include <ForceField/UFF/Contribs.h>
#include <ForceField/UFF/DistanceConstraint.h>
#include <ForceField/UFF/Nonbonded.h>
#include <ForceField/UFF/TorsionAngle.h>
#include <ForceField/UFF/Inversion.h>
#include <boost/tuple/tuple.hpp>
#include <GraphMol/ROMol.h>
%}

// We have trouble with the definition as it is -- make it a Point3D_Vect instead
%ignore ForceFields::ForceField::positions(); 
%ignore ForceFields::ForceField::positions() const; 
%rename(get) ForceFields::UFF::ParamCollection::operator();

%include <ForceField/ForceField.h>
%include <ForceField/Contrib.h>
%include <ForceField/UFF/Params.h>
%include <ForceField/UFF/AngleBend.h>
%include <ForceField/UFF/BondStretch.h>
%include <ForceField/UFF/DistanceConstraint.h>
%include <ForceField/UFF/Nonbonded.h>
%include <ForceField/UFF/TorsionAngle.h>
%include <ForceField/UFF/Inversion.h>

%template(FF_Contrib_Vect) std::vector<ForceFields::ContribPtr>;
%extend ForceFields::ForceField {
  std::vector<RDGeom::Point3D *> &positions3D () {
    return ((std::vector<RDGeom::Point3D *> &) ($self)->positions());
  }

  static int UFFOptimizeMolecule(RDKit::ROMol &mol, int maxIters=200,
    double vdwThresh=10.0, int confId=-1,
    bool ignoreInterfragInteractions=true ) {

    ForceFields::ForceField *ff=RDKit::UFF::constructForceField(mol,vdwThresh, confId, ignoreInterfragInteractions);
    ff->initialize();
    int res=ff->minimize(maxIters);
    delete ff;
    return res;
  }

  %newobject UFFGetMoleculeForceField;
  static ForceFields::ForceField *UFFGetMoleculeForceField(RDKit::ROMol &mol,
    double vdwThresh=10.0,
    int confId=-1,
    bool ignoreInterfragInteractions=true ) {

    ForceFields::ForceField *ff=RDKit::UFF::constructForceField(mol,vdwThresh, confId, ignoreInterfragInteractions);
    ff->initialize();
    return ff;
  }

  static bool UFFHasAllMoleculeParams(const RDKit::ROMol &mol){
    RDKit::UFF::AtomicParamVect types;
    bool foundAll;
    boost::tie(types,foundAll)=RDKit::UFF::getAtomTypes(mol);
    return foundAll;
  }

  /* From GraphMol/ForceFieldHelpers/UFF/AtomTyper.h */
  static void UFFAddAtomChargeFlags(const RDKit::Atom *atom, std::string &atomKey, bool tolerateChargeMismatch=true) {
    RDKit::UFF::Tools::addAtomChargeFlags(atom, atomKey, tolerateChargeMismatch);
  }

  static std::string UFFGetAtomLabel(const RDKit::Atom *atom) {
    return RDKit::UFF::Tools::getAtomLabel(atom);
  }

  /* From GraphMol/ForceFieldHelpers/UFF/AtomTyper.h */
  static std::pair<std::vector<const ForceFields::UFF::AtomicParams *>,bool> UFFGetAtomTypes(const RDKit::ROMol &mol, const std::string &paramData="") {
    return RDKit::UFF::getAtomTypes(mol, paramData);
  }

  static int MMFFOptimizeMolecule(RDKit::ROMol &mol,
                                  std::string mmffVariant="MMFF94",
                                  int maxIters=200,
                                  double nonBondedThresh=100.0, int confId=-1,
                                  bool ignoreInterfragInteractions=true ) {
    int res=1;
    RDKit::MMFF::MMFFMolProperties mmffMolProperties(mol, mmffVariant);

    if (mmffMolProperties.isValid()) {
      ForceFields::ForceField *ff = RDKit::MMFF::constructForceField(mol,
        &mmffMolProperties, nonBondedThresh, confId, ignoreInterfragInteractions);
      ff->initialize();
      res = ff->minimize(maxIters);
      delete ff;
    } else {
      BOOST_LOG(rdErrorLog) << "Could not construct MMFF force field"<<std::endl;
    }

    return res;
  }

  %newobject MMFFGetMoleculeForceField;
  static ForceFields::ForceField *MMFFGetMoleculeForceField(RDKit::ROMol &mol,
                                                            std::string mmffVariant="MMFF94",
                                                            double nonBondedThresh=100.0,
                                                            int confId=-1,
                                                            bool ignoreInterfragInteractions=true ) {

    ForceFields::ForceField *ff = 0;
    RDKit::MMFF::MMFFMolProperties mmffMolProperties(mol, mmffVariant);
    if (mmffMolProperties.isValid()) {
      ff = RDKit::MMFF::constructForceField(mol,
                                            &mmffMolProperties, nonBondedThresh,
                                            confId, ignoreInterfragInteractions);
      ff->initialize();
    } else {
      BOOST_LOG(rdErrorLog) << "Could not construct MMFF force field"<<std::endl;
    }
    return ff;
  }
  static bool MMFFHasAllMoleculeParams(RDKit::ROMol &mol){
    RDKit::MMFF::MMFFMolProperties mmffMolProperties(mol);
    return mmffMolProperties.isValid();
  }

 }


