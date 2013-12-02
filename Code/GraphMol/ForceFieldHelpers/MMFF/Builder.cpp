// $Id$
//
//  Copyright (C) 2013 Paolo Tosco
//
//  Copyright (C) 2004-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <cmath>

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/Contribs.h>
#include <ForceField/MMFF/BondStretch.h>
#include "AtomTyper.h"
#include "Builder.h"
#include <stdio.h>


namespace RDKit {
  namespace MMFF {
    using namespace ForceFields::MMFF;

    namespace Tools {
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addBonds(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field)
      {
        PRECONDITION(field, "bad ForceField");
        PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
        PRECONDITION(mmffMolProperties->isValid(), "missing atom types - invalid force-field");

        std::ostream &oStream = mmffMolProperties->getMMFFOStream();
        MMFFBondCollection *mmffBond = MMFFBondCollection::getMMFFBond();
        double totalBondStretchEnergy = 0.0;
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream <<
              "\n"
              "B O N D   S T R E T C H I N G\n\n"
              "------ATOMS------   ATOM TYPES   FF     BOND     IDEAL                        FORCE\n"
              "  I        J          I    J   CLASS   LENGTH   LENGTH    DIFF.    ENERGY   CONSTANT\n"
              "-------------------------------------------------------------------------------------"
              << std::endl;
          }
          field->initialize();
        }
        for (ROMol::ConstBondIterator bi = mol.beginBonds();
          bi != mol.endBonds(); ++bi) {
          unsigned int idx1 = (*bi)->getBeginAtomIdx();
          unsigned int idx2 = (*bi)->getEndAtomIdx();
          unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx1);
          unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx2);
          unsigned int bondType = mmffMolProperties->getMMFFBondType(*bi);
          BondStretchContrib *contrib;
          bool areMMFFBondParamsEmpirical = false;
          const MMFFBond *mmffBondParams = (*mmffBond)(bondType, iAtomType, jAtomType);
          if (!mmffBondParams) {
            mmffBondParams = mmffMolProperties->getMMFFBondStretchEmpiricalRuleParams(mol, *bi);
            areMMFFBondParamsEmpirical = true;
          }
          contrib = new BondStretchContrib(field, idx1, idx2, mmffBondParams);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
          if (mmffMolProperties->getMMFFVerbosity()) {
            const Atom *iAtom = (*bi)->getBeginAtom();
            const Atom *jAtom = (*bi)->getEndAtom();
            const double r0 = Utils::calcBondRestLength(mmffBondParams);
            const double kb = Utils::calcBondForceConstant(mmffBondParams);
            const double dist = field->distance(idx1, idx2);
            const double bondStretchEnergy = Utils::calcBondStretchEnergy(r0, kb, dist);
            if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
              oStream << std::left
                << std::setw(2) << iAtom->getSymbol()
                << " #" << std::setw(5) << idx1 + 1
                << std::setw(2) << jAtom->getSymbol()
                << " #" << std::setw(5) << idx2 + 1
                << std::right << std::setw(5) << iAtomType
                << std::setw(5) << jAtomType
                << std::setw(6) << bondType << "  "
                << std::fixed << std::setprecision(3)
                << std::setw(9) << dist
                << std::setw(9) << r0
                << std::setw(9) << dist - r0
                << std::setw(10) << bondStretchEnergy
                << std::setw(10) << kb
                << std::endl;
            }
            totalBondStretchEnergy += bondStretchEnergy;
          }
          if (areMMFFBondParamsEmpirical) {
            delete mmffBondParams;
          }
        }
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::endl;
          }
          oStream << "TOTAL BOND STRETCH ENERGY      ="
            << std::right << std::setw(16) << std::fixed << std::setprecision(4)
            << totalBondStretchEnergy << std::endl;
        }
      }
      
      
      void setTwoBitCell(boost::shared_array<boost::uint8_t> &res,
        unsigned int pos, boost::uint8_t value)
      {
        unsigned int twoBitPos = pos / 4;
        unsigned int shift = 2 * (pos % 4);
        boost::uint8_t twoBitMask = 3 << shift;
        res[twoBitPos] = ((res[twoBitPos] & (~twoBitMask)) | (value << shift));
      }
      
      
      boost::uint8_t getTwoBitCell
        (boost::shared_array<boost::uint8_t> &res, unsigned int pos)
      {
        unsigned int twoBitPos = pos / 4;
        unsigned int shift = 2 * (pos % 4);
        boost::uint8_t twoBitMask = 3 << shift;
        
        return ((res[twoBitPos] & twoBitMask) >> shift);
      }
      

      // ------------------------------------------------------------------------
      //
      // the two-bit matrix returned by this contains:
      //   0: if atoms i and j are directly connected
      //   1: if atoms i and j are connected via an atom
      //   2: if atoms i and j are in a 1,4 relationship
      //   3: otherwise
      //
      //  NOTE: the caller is responsible for calling delete []
      //  on the result
      //
      // ------------------------------------------------------------------------
      boost::shared_array<boost::uint8_t> buildNeighborMatrix(const ROMol &mol)
      {
        unsigned int nAtoms = mol.getNumAtoms();
        unsigned nTwoBitCells = (nAtoms * nAtoms - 1) / 4 + 1;
        boost::shared_array<boost::uint8_t> res(new boost::uint8_t[nTwoBitCells]);
        for (unsigned int i = 0; i < nTwoBitCells; ++i) {
          res[i] = 0;
        }
        for (unsigned int i = 0; i < nAtoms; ++i) {
          unsigned int iTab = i * nAtoms;
          for (unsigned int j = i; j < nAtoms; ++j) {
            setTwoBitCell(res, iTab + j, RELATION_1_X);
            setTwoBitCell(res, i + j * nAtoms, RELATION_1_X);
          }
        }
        for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
          const Bond *bondi = mol.getBondWithIdx(i);

          setTwoBitCell(res, bondi->getBeginAtomIdx() * nAtoms + bondi->getEndAtomIdx(), RELATION_1_2);
          setTwoBitCell(res, bondi->getEndAtomIdx() * nAtoms + bondi->getBeginAtomIdx(), RELATION_1_2);

          for (unsigned int j = i + 1; j < mol.getNumBonds(); ++j) {
            const Bond *bondj = mol.getBondWithIdx(j);
            int idx1 = -1;
            int idx3 = -1;
            if (bondi->getBeginAtomIdx() == bondj->getBeginAtomIdx()) {
              idx1 = bondi->getEndAtomIdx();
              idx3 = bondj->getEndAtomIdx();
            }
            else if (bondi->getBeginAtomIdx() == bondj->getEndAtomIdx()) {
              idx1 = bondi->getEndAtomIdx();
              idx3 = bondj->getBeginAtomIdx();
            }
            else if (bondi->getEndAtomIdx() == bondj->getBeginAtomIdx()) {
              idx1 = bondi->getBeginAtomIdx();
              idx3 = bondj->getEndAtomIdx();
            }
            else if (bondi->getEndAtomIdx() == bondj->getEndAtomIdx()) {
              idx1 = bondi->getBeginAtomIdx();
              idx3 = bondj->getBeginAtomIdx();
            }
            else {
              // check if atoms i and j are in a 1,4-relationship
              if ((mol.getBondBetweenAtoms
                (bondi->getBeginAtomIdx(), bondj->getBeginAtomIdx()))
                && (getTwoBitCell(res, bondi->getEndAtomIdx() * nAtoms + bondj->getEndAtomIdx()) == RELATION_1_X)) {
                setTwoBitCell(res, bondi->getEndAtomIdx() * nAtoms + bondj->getEndAtomIdx(), RELATION_1_4);
                setTwoBitCell(res, bondj->getEndAtomIdx() * nAtoms + bondi->getEndAtomIdx(), RELATION_1_4);
              }
              else if ((mol.getBondBetweenAtoms
                (bondi->getBeginAtomIdx(), bondj->getEndAtomIdx()))
                && (getTwoBitCell(res, bondi->getEndAtomIdx() * nAtoms + bondj->getBeginAtomIdx()) == RELATION_1_X)) {
                setTwoBitCell(res, bondi->getEndAtomIdx() * nAtoms + bondj->getBeginAtomIdx(), RELATION_1_4);
                setTwoBitCell(res, bondj->getBeginAtomIdx() * nAtoms + bondi->getEndAtomIdx(), RELATION_1_4);
              }
              else if ((mol.getBondBetweenAtoms
                (bondi->getEndAtomIdx(), bondj->getBeginAtomIdx()))
                && (getTwoBitCell(res, bondi->getBeginAtomIdx() * nAtoms + bondj->getEndAtomIdx()) == RELATION_1_X)) {
                setTwoBitCell(res, bondi->getBeginAtomIdx() * nAtoms + bondj->getEndAtomIdx(), RELATION_1_4);
                setTwoBitCell(res, bondj->getEndAtomIdx() * nAtoms + bondi->getBeginAtomIdx(), RELATION_1_4);
              }
              else if ((mol.getBondBetweenAtoms
                (bondi->getEndAtomIdx(), bondj->getEndAtomIdx()))
                && (getTwoBitCell(res, bondi->getBeginAtomIdx() * nAtoms + bondj->getBeginAtomIdx()) == RELATION_1_X)) {
                setTwoBitCell(res, bondi->getBeginAtomIdx() * nAtoms + bondj->getBeginAtomIdx(), RELATION_1_4);
                setTwoBitCell(res, bondj->getBeginAtomIdx() * nAtoms + bondi->getBeginAtomIdx(), RELATION_1_4);
              }
            }
            if (idx1 > -1) {
              setTwoBitCell(res, idx1 * nAtoms + idx3, RELATION_1_3);
              setTwoBitCell(res, idx3 * nAtoms + idx1, RELATION_1_3);
            }
          }
        }
        return res;
      }
      
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addAngles(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field)
      {
        PRECONDITION(field, "bad ForceField");
        PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
        PRECONDITION(mmffMolProperties->isValid(), "missing atom types - invalid force-field");

        std::ostream &oStream = mmffMolProperties->getMMFFOStream();
        unsigned int idx[3];
        MMFFAngleCollection *mmffAngle = MMFFAngleCollection::getMMFFAngle();
        MMFFBondCollection *mmffBond = MMFFBondCollection::getMMFFBond();
        MMFFPropCollection *mmffProp = MMFFPropCollection::getMMFFProp();
        AngleBendContrib *contrib;
        ROMol::ADJ_ITER nbr1Idx;
        ROMol::ADJ_ITER end1Nbrs;
        ROMol::ADJ_ITER nbr2Idx;
        ROMol::ADJ_ITER end2Nbrs;

        unsigned int nAtoms = mol.getNumAtoms();
        double totalAngleBendEnergy = 0.0;
        RDGeom::PointPtrVect points;
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream <<
              "\n"
              "A N G L E   B E N D I N G\n\n"
              "-----------ATOMS-----------    ATOM TYPES      FF    VALENCE    IDEAL                          FORCE\n"
              "  I        J        K          I    J    K   CLASS    ANGLE     ANGLE      DIFF.    ENERGY   CONSTANT\n"
              "------------------------------------------------------------------------------------------------------"
              << std::endl;
          }
          field->initialize();
          points = field->positions();
        }
        for (idx[1] = 0; idx[1] < nAtoms; ++idx[1]) {
          const Atom *jAtom = mol.getAtomWithIdx(idx[1]);
          if (jAtom->getDegree() == 1) {
            continue;
          }
          unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
          const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(jAtomType);
          boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(jAtom);
          for (; nbr1Idx != end1Nbrs; ++nbr1Idx) {
            const Atom *iAtom = mol[*nbr1Idx].get();
            idx[0] = iAtom->getIdx();
            boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(jAtom);
            for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
              if (nbr2Idx < (nbr1Idx + 1)) {
                continue;
              }
              const Atom *kAtom = mol[*nbr2Idx].get();
              idx[2] = kAtom->getIdx();
              unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[0]);
              unsigned int kAtomType = mmffMolProperties->getMMFFAtomType(idx[2]);
              unsigned int angleType = mmffMolProperties->getMMFFAngleType
                (mol, idx[0], idx[1], idx[2]);
              bool areMMFFBondParamsEmpirical[2] = { false, false };
              bool areMMFFAngleParamsEmpirical = false;
              const MMFFAngle *mmffAngleParams = (*mmffAngle)(angleType, iAtomType, jAtomType, kAtomType);
              const MMFFBond *mmffBondParams[2] = { NULL, NULL };
              if ((!mmffAngleParams) || (isDoubleZero(mmffAngleParams->ka))) {
                areMMFFAngleParamsEmpirical = true;
                for (unsigned int i = 0; i < 2; ++i) {
                  const Bond *bond = mol.getBondBetweenAtoms(idx[i], idx[i + 1]);
                  unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[i]);
                  unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[i + 1]);
                  unsigned int bondType = mmffMolProperties->getMMFFBondType(bond);
                  mmffBondParams[i] = (*mmffBond)(bondType, iAtomType, jAtomType);
                  if (!(mmffBondParams[i])) {
                    mmffBondParams[i] = mmffMolProperties->getMMFFBondStretchEmpiricalRuleParams(mol, bond);
                    areMMFFBondParamsEmpirical[i] = true;
                  }
                }
                mmffAngleParams = getMMFFAngleBendEmpiricalRuleParams
                  (mol, mmffAngleParams, mmffPropParamsCentralAtom,
                  mmffBondParams[0], mmffBondParams[1], idx[0], idx[1], idx[2]);
              }
              contrib = new AngleBendContrib(field, idx[0], idx[1], idx[2],
                mmffAngleParams, mmffPropParamsCentralAtom);
              field->contribs().push_back(ForceFields::ContribPtr(contrib));
              if (mmffMolProperties->getMMFFVerbosity()) {
                const double theta0 = Utils::calcAngleRestValue(mmffAngleParams);
                const double ka = Utils::calcAngleForceConstant(mmffAngleParams);
                const RDGeom::Point3D p1((*(points[idx[0]]))[0],
                  (*(points[idx[0]]))[1], (*(points[idx[0]]))[2]);
                const RDGeom::Point3D p2((*(points[idx[1]]))[0],
                  (*(points[idx[1]]))[1], (*(points[idx[1]]))[2]);
                const RDGeom::Point3D p3((*(points[idx[2]]))[0],
                  (*(points[idx[2]]))[1], (*(points[idx[2]]))[2]);
                const double cosTheta = Utils::calcCosTheta(p1, p2, p3,
                  field->distance(idx[0], idx[1]), field->distance(idx[1], idx[2]));
                const double theta = RAD2DEG * acos(cosTheta);
                const double angleBendEnergy = Utils::calcAngleBendEnergy
                  (theta0, ka, mmffPropParamsCentralAtom->linh, cosTheta);
                if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
                  oStream << std::left
                    << std::setw(2) << iAtom->getSymbol()
                    << " #" << std::setw(5) << idx[0] + 1
                    << std::setw(2) << jAtom->getSymbol()
                    << " #" << std::setw(5) << idx[1] + 1
                    << std::setw(2) << kAtom->getSymbol()
                    << " #" << std::setw(5) << idx[2] + 1
                    << std::right << std::setw(5) << iAtomType
                    << std::setw(5) << jAtomType
                    << std::setw(5) << kAtomType
                    << std::setw(6) << angleType << "  "
                    << std::fixed << std::setprecision(3)
                    << std::setw(10) << theta
                    << std::setw(10) << theta0
                    << std::setw(10) << theta - theta0
                    << std::setw(10) << angleBendEnergy
                    << std::setw(10) << ka
                    << std::endl;
                }
                totalAngleBendEnergy += angleBendEnergy;
              }
              if (areMMFFAngleParamsEmpirical) {
                delete mmffAngleParams;
              }
              if (areMMFFBondParamsEmpirical[0]) {
                delete mmffBondParams[0];
              }
              if (areMMFFBondParamsEmpirical[1]) {
                delete mmffBondParams[1];
              }
            }
          }
        }
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::endl;
          }
          oStream << "TOTAL ANGLE BEND ENERGY        ="
            << std::right << std::setw(16) << std::fixed << std::setprecision(4)
            << totalAngleBendEnergy << std::endl;
        }
      }

      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addStretchBend(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field)
      {
        PRECONDITION(field, "bad ForceField");
        PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
        PRECONDITION(mmffMolProperties->isValid(), "missing atom types - invalid force-field");

        std::ostream &oStream = mmffMolProperties->getMMFFOStream();
        unsigned int idx[3];
        unsigned int bondType[2];
        MMFFStbnCollection *mmffStbn = MMFFStbnCollection::getMMFFStbn();
        MMFFDfsbCollection *mmffDfsb = MMFFDfsbCollection::getMMFFDfsb();
        MMFFAngleCollection *mmffAngle = MMFFAngleCollection::getMMFFAngle();
        MMFFBondCollection *mmffBond = MMFFBondCollection::getMMFFBond();
        MMFFPropCollection *mmffProp = MMFFPropCollection::getMMFFProp();
        std::pair<bool, const MMFFStbn *> mmffStbnParams;
        const MMFFAngle *mmffAngleParams = NULL;
        StretchBendContrib *contrib;
        ROMol::ADJ_ITER nbr1Idx;
        ROMol::ADJ_ITER end1Nbrs;
        ROMol::ADJ_ITER nbr2Idx;
        ROMol::ADJ_ITER end2Nbrs;

        unsigned int nAtoms = mol.getNumAtoms();
        double totalStretchBendEnergy = 0.0;
        RDGeom::PointPtrVect points;
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream <<
              "\n"
              "S T R E T C H   B E N D I N G\n\n"
              "-----------ATOMS-----------    ATOM TYPES      FF    VALENCE     DELTA     DELTA     DELTA               F CON\n"
              "  I        J        K          I    J    K   CLASS    ANGLE      ANGLE     R(I,J)    R(J,K)   ENERGY    I-J (J-K)\n"
              "------------------------------------------------------------------------------------------------------------------"
              << std::endl;
          }
          field->initialize();
          points = field->positions();
        }
        for (idx[1] = 0; idx[1] < nAtoms; ++idx[1]) {
          const Atom *jAtom = mol.getAtomWithIdx(idx[1]);
          if (jAtom->getDegree() == 1) {
            continue;
          }
          unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
          const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(jAtomType);
          if (mmffPropParamsCentralAtom->linh) {
            continue;
          }
          boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(jAtom);
          unsigned int i = 0;
          for (; nbr1Idx != end1Nbrs; ++nbr1Idx) {
            const Atom *iAtom = mol[*nbr1Idx].get();
            boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(jAtom);
            unsigned int j = 0;
            for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
              const Atom *kAtom = mol[*nbr2Idx].get();
              if (j < (i + 1)) {
                ++j;
                continue;
              }
              idx[0] = iAtom->getIdx();
              idx[2] = kAtom->getIdx();
              unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[0]);
              unsigned int kAtomType = mmffMolProperties->getMMFFAtomType(idx[2]);
              unsigned int angleType = mmffMolProperties->getMMFFAngleType
                (mol, idx[0], idx[1], idx[2]);
              bool areMMFFBondParamsEmpirical[2] = { false, false };
              bool areMMFFAngleParamsEmpirical = false;
              mmffAngleParams = (*mmffAngle)(angleType, iAtomType, jAtomType, kAtomType);
              const MMFFBond *mmffBondParams[2] = { NULL, NULL };
              for (unsigned int i = 0; i < 2; ++i) {
                const Bond *bond = mol.getBondBetweenAtoms(idx[i], idx[i + 1]);
                unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[i]);
                unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[i + 1]);
                bondType[i] = mmffMolProperties->getMMFFBondType(bond);
                mmffBondParams[i] = (*mmffBond)(bondType[i], iAtomType, jAtomType);
                if (!(mmffBondParams[i])) {
                  mmffBondParams[i] = mmffMolProperties->getMMFFBondStretchEmpiricalRuleParams(mol, bond);
                  areMMFFBondParamsEmpirical[i] = true;
                }
              }
              if ((!mmffAngleParams) || isDoubleZero(mmffAngleParams->ka)) {
                mmffAngleParams = getMMFFAngleBendEmpiricalRuleParams
                  (mol, mmffAngleParams, mmffPropParamsCentralAtom,
                  mmffBondParams[0], mmffBondParams[1], idx[0], idx[1], idx[2]);
                areMMFFAngleParamsEmpirical = true;
              }
              unsigned int stretchBendType = MMFF::getMMFFStretchBendType
                (angleType, (iAtomType <= kAtomType) ? bondType[0] : bondType[1],
                (iAtomType < kAtomType) ? bondType[1] : bondType[0]);
              mmffStbnParams = mmffStbn->getMMFFStbnParams(stretchBendType,
                bondType[0], bondType[1], iAtomType, jAtomType, kAtomType);
              if (!(mmffStbnParams.second)) {
                mmffStbnParams = mmffDfsb->getMMFFDfsbParams(getPeriodicTableRow(iAtom->getAtomicNum()),
                  getPeriodicTableRow(jAtom->getAtomicNum()), getPeriodicTableRow(kAtom->getAtomicNum()));
              }
              if (!(isDoubleZero((mmffStbnParams.second)->kbaIJK)
                && isDoubleZero((mmffStbnParams.second)->kbaKJI))) {
                contrib = new StretchBendContrib(field, idx[0], idx[1], idx[2],
                  mmffStbnParams, mmffAngleParams, mmffBondParams[0], mmffBondParams[1]);
                field->contribs().push_back(ForceFields::ContribPtr(contrib));
                if (mmffMolProperties->getMMFFVerbosity()) {
                  const double theta0 = Utils::calcAngleRestValue(mmffAngleParams);
                  const double restLen1 = Utils::calcBondRestLength(mmffBondParams[0]);
                  const double restLen2 = Utils::calcBondRestLength(mmffBondParams[1]);
                  const double dist1 = field->distance(idx[0], idx[1]);
                  const double dist2 = field->distance(idx[1], idx[2]);
                  const RDGeom::Point3D p1((*(points[idx[0]]))[0],
                    (*(points[idx[0]]))[1], (*(points[idx[0]]))[2]);
                  const RDGeom::Point3D p2((*(points[idx[1]]))[0],
                    (*(points[idx[1]]))[1], (*(points[idx[1]]))[2]);
                  const RDGeom::Point3D p3((*(points[idx[2]]))[0],
                    (*(points[idx[2]]))[1], (*(points[idx[2]]))[2]);
                  const double cosTheta = Utils::calcCosTheta(p1, p2, p3, dist1, dist2);
                  const double theta = RAD2DEG * acos(cosTheta);
                  const std::pair<double, double> forceConstants =
                    Utils::calcStbnForceConstants(mmffStbnParams);
                  const std::pair<double, double> stretchBendEnergies
                    = Utils::calcStretchBendEnergy(dist1 - restLen1,
                    dist2 - restLen2, theta - theta0, forceConstants);
                  if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
                    if (!isDoubleZero(forceConstants.first)) {
                      oStream << std::left
                        << std::setw(2) << iAtom->getSymbol()
                        << " #" << std::setw(5) << idx[0] + 1
                        << std::setw(2) << jAtom->getSymbol()
                        << " #" << std::setw(5) << idx[1] + 1
                        << std::setw(2) << kAtom->getSymbol()
                        << " #" << std::setw(5) << idx[2] + 1
                        << std::right << std::setw(5) << iAtomType
                        << std::setw(5) << jAtomType
                        << std::setw(5) << kAtomType
                        << std::setw(6) << stretchBendType << "  "
                        << std::fixed << std::setprecision(3)
                        << std::setw(10) << theta
                        << std::setw(10) << theta - theta0
                        << std::setw(10) << dist1 - restLen1
                        << std::setw(10) << dist2 - restLen2
                        << std::setw(10) << stretchBendEnergies.first
                        << std::setw(10) << forceConstants.first
                        << std::endl;
                    }
                    if (!isDoubleZero(forceConstants.second)) {
                      oStream << std::left
                        << std::setw(2) << kAtom->getSymbol()
                        << " #" << std::setw(5) << idx[2] + 1
                        << std::setw(2) << jAtom->getSymbol()
                        << " #" << std::setw(5) << idx[1] + 1
                        << std::setw(2) << iAtom->getSymbol()
                        << " #" << std::setw(5) << idx[0] + 1
                        << std::right << std::setw(5) << kAtomType
                        << std::setw(5) << jAtomType
                        << std::setw(5) << iAtomType
                        << std::setw(6) << stretchBendType << "  "
                        << std::fixed << std::setprecision(3)
                        << std::setw(10) << theta
                        << std::setw(10) << theta - theta0
                        << std::setw(10) << dist1 - restLen1
                        << std::setw(10) << dist2 - restLen2
                        << std::setw(10) << stretchBendEnergies.second
                        << std::setw(10) << forceConstants.second
                        << std::endl;
                    }
                  }
                  totalStretchBendEnergy += (stretchBendEnergies.first + stretchBendEnergies.second);
                }
              }
              if (areMMFFAngleParamsEmpirical) {
                delete mmffAngleParams;
              }
              if (areMMFFBondParamsEmpirical[0]) {
                delete mmffBondParams[0];
              }
              if (areMMFFBondParamsEmpirical[1]) {
                delete mmffBondParams[1];
              }
              ++j;
            }
            ++i;
          }
        }
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::endl;
          }
          oStream << "TOTAL STRETCH-BEND ENERGY      ="
            << std::right << std::setw(16) << std::fixed << std::setprecision(4)
            << totalStretchBendEnergy << std::endl;
        }
      }

      void addOop(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field)
      {
        PRECONDITION(field, "bad ForceField");
        PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
        PRECONDITION(mmffMolProperties->isValid(), "missing atom types - invalid force-field");
        
        std::ostream &oStream = mmffMolProperties->getMMFFOStream();
        unsigned int idx[4];
        unsigned int atomType[4];
        unsigned int n[4];
        const Atom *atom[4];
        ROMol::ADJ_ITER nbrIdx;
        ROMol::ADJ_ITER endNbrs;

        MMFFOopCollection *mmffOop = MMFFOopCollection::getMMFFOop
          (mmffMolProperties->getMMFFVariant() == "MMFF94s");
        OopBendContrib *contrib;
        const MMFFOop *mmffOopParams;
        double totalOopBendEnergy = 0.0;
        RDGeom::PointPtrVect points;
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream <<
              "\n"
              "O U T - O F - P L A N E   B E N D I N G\n\n"
              "--------------ATOMS---------------         ATOM TYPES         OOP                FORCE\n"
              "  I        J        K        L          I    J    K    L     ANGLE    ENERGY   CONSTANT\n"
              "------------------------------------------------------------------------------------------"
              << std::endl;
          }
          field->initialize();
          points = field->positions();
        }
        for (idx[1] = 0; idx[1] < mol.getNumAtoms(); ++idx[1]) {
          atom[1] = mol.getAtomWithIdx(idx[1]);
          if (atom[1]->getDegree() != 3) {
            continue;
          }
          atomType[1] = mmffMolProperties->getMMFFAtomType(idx[1]);
          boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom[1]);
          unsigned int i = 0;
          for (; nbrIdx != endNbrs; ++nbrIdx) {
            atom[i] = mol[*nbrIdx].get();
            idx[i] = atom[i]->getIdx();
            atomType[i] = mmffMolProperties->getMMFFAtomType(idx[i]);
            if (!i) {
              ++i;
            }
            ++i;
          }
          mmffOopParams = (*mmffOop)(atomType[0], atomType[1], atomType[2], atomType[3]);
          // if no parameters could be found, we exclude this term (SURDOX02)
          if (!mmffOopParams) {
            continue;
          }
          for (unsigned int i = 0; i < 3; ++i) {
            n[1] = 1;
            switch (i) {
              case 0:
                n[0] = 0;
                n[2] = 2;
                n[3] = 3;
              break;
              
              case 1:
                n[0] = 0;
                n[2] = 3;
                n[3] = 2;
              break;
              
              case 2:
                n[0] = 2;
                n[2] = 3;
                n[3] = 0;
              break;
            }
            contrib = new OopBendContrib
              (field, idx[n[0]], idx[n[1]], idx[n[2]], idx[n[3]], mmffOopParams);
            field->contribs().push_back(ForceFields::ContribPtr(contrib));
            if (mmffMolProperties->getMMFFVerbosity()) {
              const RDGeom::Point3D p1((*(points[idx[n[0]]]))[0],
                (*(points[idx[n[0]]]))[1], (*(points[idx[n[0]]]))[2]);
              const RDGeom::Point3D p2((*(points[idx[n[1]]]))[0],
                (*(points[idx[n[1]]]))[1], (*(points[idx[n[1]]]))[2]);
              const RDGeom::Point3D p3((*(points[idx[n[2]]]))[0],
                (*(points[idx[n[2]]]))[1], (*(points[idx[n[2]]]))[2]);
              const RDGeom::Point3D p4((*(points[idx[n[3]]]))[0],
                (*(points[idx[n[3]]]))[1], (*(points[idx[n[3]]]))[2]);
              const double chi = Utils::calcOopChi(p1, p2, p3, p4);
              const double koop = Utils::calcOopBendForceConstant(mmffOopParams);
              const double oopBendEnergy = Utils::calcOopBendEnergy(chi, koop);
              if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
                oStream << std::left
                  << std::setw(2) << atom[0]->getSymbol()
                  << " #" << std::setw(5) << idx[n[0]] + 1
                  << std::setw(2) << atom[1]->getSymbol()
                  << " #" << std::setw(5) << idx[n[1]] + 1
                  << std::setw(2) << atom[2]->getSymbol()
                  << " #" << std::setw(5) << idx[n[2]] + 1
                  << std::setw(2) << atom[3]->getSymbol()
                  << " #" << std::setw(5) << idx[n[3]] + 1
                  << std::right << std::setw(5) << atomType[n[0]]
                  << std::setw(5) << atomType[n[1]]
                  << std::setw(5) << atomType[n[2]]
                  << std::setw(5) << atomType[n[3]]
                  << std::fixed << std::setprecision(3)
                  << std::setw(10) << chi
                  << std::setw(10) << oopBendEnergy
                  << std::setw(10) << koop
                  << std::endl;
              }
              totalOopBendEnergy += oopBendEnergy;
            }
          }
        }
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::endl;
          }
          oStream << "TOTAL OUT-OF-PLANE BEND ENERGY ="
            << std::right << std::setw(16) << std::fixed << std::setprecision(4)
            << totalOopBendEnergy << std::endl;
        }
      }
                       
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addTorsions(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field,
        std::string torsionBondSmarts)
      {
        PRECONDITION(field, "bad ForceField");
        PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
        PRECONDITION(mmffMolProperties->isValid(), "missing atom types - invalid force-field");

        std::ostream &oStream = mmffMolProperties->getMMFFOStream();
        ROMol::ADJ_ITER nbr1Idx;
        ROMol::ADJ_ITER end1Nbrs;
        ROMol::ADJ_ITER nbr2Idx;
        ROMol::ADJ_ITER end2Nbrs;
        MMFFTorCollection *mmffTor = MMFFTorCollection::getMMFFTor
          (mmffMolProperties->getMMFFVariant() == "MMFF94s");
        double totalTorsionEnergy = 0.0;
        RDGeom::PointPtrVect points;
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream <<
              "\n"
              "T O R S I O N A L\n\n"
              "--------------ATOMS---------------      ---ATOM TYPES---     FF     TORSION              -----FORCE CONSTANTS-----\n"
              "  I        J        K        L          I    J    K    L   CLASS     ANGLE    ENERGY       V1        V2        V3\n"
              "------------------------------------------------------------------------------------------------------------------"
              << std::endl;
          }
          field->initialize();
          points = field->positions();
        }
        std::vector<MatchVectType> matchVect;
        ROMol *query = SmartsToMol(torsionBondSmarts);
        TEST_ASSERT(query);
        unsigned int nHits = SubstructMatch(mol, *query, matchVect);
        delete query;

        for (unsigned int i = 0; i < nHits; ++i) {
          MatchVectType match = matchVect[i];
          TEST_ASSERT(match.size() == 2);
          int idx2 = match[0].second;
          int idx3 = match[1].second;
          const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
          TEST_ASSERT(bond);
          const Atom *jAtom = mol.getAtomWithIdx(idx2);
          const Atom *kAtom = mol.getAtomWithIdx(idx3);
          if (((jAtom->getHybridization() == Atom::SP2) || (jAtom->getHybridization() == Atom::SP3))
            && ((kAtom->getHybridization() == Atom::SP2) || (kAtom->getHybridization() == Atom::SP3))) {
            ROMol::OEDGE_ITER beg1, end1;
            boost::tie(beg1, end1) = mol.getAtomBonds(jAtom);
            while (beg1 != end1) {
              const Bond *tBond1 = mol[*beg1].get();
              if (tBond1 != bond) {
                int idx1 = tBond1->getOtherAtomIdx(idx2);
                ROMol::OEDGE_ITER beg2, end2;
                boost::tie(beg2, end2) = mol.getAtomBonds(kAtom);
                while (beg2 != end2) {
                  const Bond *tBond2 = mol[*beg2].get();
                  if ((tBond2 != bond) && (tBond2 != tBond1)) {
                    int idx4 = tBond2->getOtherAtomIdx(idx3);
                    // make sure this isn't a three-membered ring:
                    if (idx4 != idx1) {
                      // we now have a torsion involving atoms (bonds):
                      //  bIdx - (tBond1) - idx1 - (bond) - idx2 - (tBond2) - eIdx
                      const Atom *iAtom = mol.getAtomWithIdx(idx1);
                      const Atom *lAtom = mol.getAtomWithIdx(idx4);
                      TorsionAngleContrib *contrib;
                      unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx1);
                      unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx2);
                      unsigned int kAtomType = mmffMolProperties->getMMFFAtomType(idx3);
                      unsigned int lAtomType = mmffMolProperties->getMMFFAtomType(idx4);
                      const std::pair<unsigned int, unsigned int> torTypePair =
                        mmffMolProperties->getMMFFTorsionType(mol, idx1, idx2, idx3, idx4);
                      bool areMMFFTorParamsEmpirical = false;
                      const std::pair<const unsigned int, const MMFFTor *> mmffTorPair =
                        mmffTor->getMMFFTorParams(torTypePair, iAtomType, jAtomType, kAtomType, lAtomType);
                      unsigned int torType = (mmffTorPair.first ? mmffTorPair.first : torTypePair.first);
                      const MMFFTor *mmffTorParams = mmffTorPair.second;
                      if (!mmffTorParams) {
                        torType = torTypePair.first;
                        mmffTorParams = mmffMolProperties->getMMFFTorsionEmpiricalRuleParams(mol, idx2, idx3);
                        areMMFFTorParamsEmpirical = true;
                      }
                      if (!(isDoubleZero(mmffTorParams->V1) && isDoubleZero(mmffTorParams->V2)
                        && isDoubleZero(mmffTorParams->V3))) {
                        contrib = new TorsionAngleContrib(field, idx1, idx2, idx3, idx4, mmffTorParams);
                        field->contribs().push_back(ForceFields::ContribPtr(contrib));
                        if (mmffMolProperties->getMMFFVerbosity()) {
                          const RDGeom::Point3D p1((*(points[idx1]))[0],
                            (*(points[idx1]))[1], (*(points[idx1]))[2]);
                          const RDGeom::Point3D p2((*(points[idx2]))[0],
                            (*(points[idx2]))[1], (*(points[idx2]))[2]);
                          const RDGeom::Point3D p3((*(points[idx3]))[0],
                            (*(points[idx3]))[1], (*(points[idx3]))[2]);
                          const RDGeom::Point3D p4((*(points[idx4]))[0],
                            (*(points[idx4]))[1], (*(points[idx4]))[2]);
                          const double cosPhi = Utils::calcTorsionCosPhi(p1, p2, p3, p4);
                          const boost::tuple<double, double, double> forceConstants =
                            Utils::calcTorsionForceConstant(mmffTorParams);
                          const double torsionEnergy = Utils::calcTorsionEnergy
                            (boost::tuples::get<0>(forceConstants), boost::tuples::get<1>(forceConstants),
                            boost::tuples::get<2>(forceConstants), cosPhi);
                          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
                            oStream << std::left
                              << std::setw(2) << iAtom->getSymbol()
                              << " #" << std::setw(5) << idx1 + 1
                              << std::setw(2) << jAtom->getSymbol()
                              << " #" << std::setw(5) << idx2 + 1
                              << std::setw(2) << kAtom->getSymbol()
                              << " #" << std::setw(5) << idx3 + 1
                              << std::setw(2) << lAtom->getSymbol()
                              << " #" << std::setw(5) << idx4 + 1
                              << std::right << std::setw(5) << iAtomType
                              << std::setw(5) << jAtomType
                              << std::setw(5) << kAtomType
                              << std::setw(5) << lAtomType
                              << std::setw(6) << torType << "  "
                              << std::fixed << std::setprecision(3)
                              << std::setw(10) << RAD2DEG * acos(cosPhi)
                              << std::setw(10) << torsionEnergy
                              << std::setw(10) << boost::tuples::get<0>(forceConstants)
                              << std::setw(10) << boost::tuples::get<1>(forceConstants)
                              << std::setw(10) << boost::tuples::get<2>(forceConstants)
                              << std::endl;
                          }
                          totalTorsionEnergy += torsionEnergy;
                        }
                      }
                      if(areMMFFTorParamsEmpirical){
                        delete mmffTorParams;
                      }
                    }
                  }
                  beg2++;
                }
              }
              beg1++;
            }
          }
        }
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::endl;
          }
          oStream << "TOTAL TORSIONAL ENERGY         ="
            << std::right << std::setw(16) << std::fixed << std::setprecision(4)
            << totalTorsionEnergy << std::endl;
        }
      }
                       
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addVdW(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
        ForceFields::ForceField *field, boost::shared_array<boost::uint8_t> neighborMatrix,
        double nonBondedThresh, bool ignoreInterfragInteractions)
      {
        PRECONDITION(field, "bad ForceField");
        PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
        PRECONDITION(mmffMolProperties->isValid(), "missing atom types - invalid force-field");

        std::ostream &oStream = mmffMolProperties->getMMFFOStream();
        MMFFVdWCollection *mmffVdW = MMFFVdWCollection::getMMFFVdW();
        INT_VECT fragMapping;
        if (ignoreInterfragInteractions) {
          std::vector<ROMOL_SPTR> molFrags = MolOps::getMolFrags(mol, true, &fragMapping);
        }

        unsigned int nAtoms = mol.getNumAtoms();
        double totalVdWEnergy = 0.0;
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream <<
              "\n"
              "V A N   D E R   W A A L S\n\n"
              "------ATOMS------   ATOM TYPES                                 WELL\n"
              "  I        J          I    J    DISTANCE   ENERGY     R*      DEPTH\n"
              "--------------------------------------------------------------------"
              << std::endl;
          }
        }
        const Conformer &conf = mol.getConformer(confId);
        for (unsigned int i = 0; i < nAtoms; ++i) {
          for (unsigned int j = i + 1; j < nAtoms; ++j) {
            if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
              continue;
            }
            if (getTwoBitCell(neighborMatrix, i * nAtoms + j) >= RELATION_1_4) {
              double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
              if (dist > nonBondedThresh) {
                continue;
              }
              const unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(i);
              const unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(j);
              const MMFFVdW *mmffVdWParamsIAtom = (*mmffVdW)(iAtomType);
              const MMFFVdW *mmffVdWParamsJAtom = (*mmffVdW)(jAtomType);
              VdWContrib *contrib;
              contrib = new VdWContrib(field, i, j,
                mmffVdW, mmffVdWParamsIAtom, mmffVdWParamsJAtom);
              field->contribs().push_back(ForceFields::ContribPtr(contrib));
              if (mmffMolProperties->getMMFFVerbosity()) {
                const Atom *iAtom = mol.getAtomWithIdx(i);
                const Atom *jAtom = mol.getAtomWithIdx(j);
                double R_star_ij = Utils::calcUnscaledVdWMinimum
                  (mmffVdW, mmffVdWParamsIAtom, mmffVdWParamsJAtom);
                double wellDepth = Utils::calcUnscaledVdWWellDepth
                  (R_star_ij, mmffVdWParamsIAtom, mmffVdWParamsJAtom);
                Utils::scaleVdWParams(R_star_ij, wellDepth,
                  mmffVdW, mmffVdWParamsIAtom, mmffVdWParamsJAtom);
                const double vdWEnergy = Utils::calcVdWEnergy(dist, R_star_ij, wellDepth);
                if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
                  oStream << std::left
                    << std::setw(2) << iAtom->getSymbol()
                    << " #" << std::setw(5) << i + 1
                    << std::setw(2) << jAtom->getSymbol()
                    << " #" << std::setw(5) << j + 1
                    << std::right << std::setw(5) << iAtomType
                    << std::setw(5) << jAtomType << "  "
                    << std::fixed << std::setprecision(3)
                    << std::setw(9) << dist
                    << std::setw(10) << vdWEnergy
                    << std::setw(9) << R_star_ij
                    << std::setw(9) << wellDepth
                    << std::endl;
                }
                totalVdWEnergy += vdWEnergy;
              }
            }
          }
        }
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::endl;
          }
          oStream << "TOTAL VAN DER WAALS ENERGY     ="
            << std::right << std::setw(16) << std::fixed << std::setprecision(4)
            << totalVdWEnergy << std::endl;
        }
      }

      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addEle(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
        ForceFields::ForceField *field, boost::shared_array<boost::uint8_t> neighborMatrix,
        double nonBondedThresh, bool ignoreInterfragInteractions)
      {
        PRECONDITION(field, "bad ForceField");
        PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
        PRECONDITION(mmffMolProperties->isValid(), "missing atom types - invalid force-field");

        std::ostream &oStream = mmffMolProperties->getMMFFOStream();
        INT_VECT fragMapping;
        if (ignoreInterfragInteractions) {
          std::vector<ROMOL_SPTR> molFrags = MolOps::getMolFrags(mol, true, &fragMapping);
        }
        unsigned int nAtoms = mol.getNumAtoms();
        double totalEleEnergy = 0.0;
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream <<
              "\n"
              "E L E C T R O S T A T I C\n\n"
              "------ATOMS------   ATOM TYPES\n"
              "  I        J          I    J    DISTANCE   ENERGY\n"
              "--------------------------------------------------"
              << std::endl;
          }
        }
        const Conformer &conf = mol.getConformer(confId);
        double dielConst = mmffMolProperties->getMMFFDielectricConstant();
        boost::uint8_t dielModel = mmffMolProperties->getMMFFDielectricModel();
        for (unsigned int i = 0; i < nAtoms; ++i) {
          for (unsigned int j = i + 1; j < nAtoms; ++j) {
            if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
              continue;
            }
            if ((getTwoBitCell(neighborMatrix, i * nAtoms + j) >= RELATION_1_4)
              && (!isDoubleZero(mmffMolProperties->getMMFFPartialCharge(i)))
              && (!isDoubleZero(mmffMolProperties->getMMFFPartialCharge(j)))) {
              double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
              if (dist > nonBondedThresh) {
                continue;
              }
              EleContrib *contrib;
              double chargeTerm = mmffMolProperties->getMMFFPartialCharge(i)
                * mmffMolProperties->getMMFFPartialCharge(j) / dielConst;
              contrib = new EleContrib(field, i, j, chargeTerm, dielModel,
                getTwoBitCell(neighborMatrix, i * nAtoms + j) == RELATION_1_4);
              field->contribs().push_back(ForceFields::ContribPtr(contrib));
              if (mmffMolProperties->getMMFFVerbosity()) {
                const unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(i);
                const unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(j);
                const Atom *iAtom = mol.getAtomWithIdx(i);
                const Atom *jAtom = mol.getAtomWithIdx(j);
                const double eleEnergy = Utils::calcEleEnergy
                  (i, j, dist, chargeTerm, dielModel,
                  getTwoBitCell(neighborMatrix, i * nAtoms + j) == RELATION_1_4);
                if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
                  oStream << std::left
                    << std::setw(2) << iAtom->getSymbol()
                    << " #" << std::setw(5) << i + 1
                    << std::setw(2) << jAtom->getSymbol()
                    << " #" << std::setw(5) << j + 1
                    << std::right << std::setw(5) << iAtomType
                    << std::setw(5) << jAtomType << "  "
                    << std::fixed << std::setprecision(3)
                    << std::setw(9) << dist
                    << std::setw(10) << eleEnergy
                    << std::endl;
                }
                totalEleEnergy += eleEnergy;
              }
            }
          }
        }
        if (mmffMolProperties->getMMFFVerbosity()) {
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::endl;
          }
          oStream << "TOTAL ELECTROSTATIC ENERGY     ="
            << std::right << std::setw(16) << std::fixed << std::setprecision(4)
            << totalEleEnergy << std::endl;
        }
      }

    } // end of namespace Tools
    
    // ------------------------------------------------------------------------
    //
    //
    //
    // ------------------------------------------------------------------------
    ForceFields::ForceField *constructForceField(ROMol &mol,
      double nonBondedThresh, int confId, bool ignoreInterfragInteractions)
    {
      MMFFMolProperties mmffMolProperties(mol);
      PRECONDITION(mmffMolProperties.isValid(), "missing atom types - invalid force-field");
      ForceFields::ForceField *res = constructForceField(mol,
        &mmffMolProperties, nonBondedThresh, confId, ignoreInterfragInteractions);
        
      return res;
    }
    
    // ------------------------------------------------------------------------
    //
    //
    //
    // ------------------------------------------------------------------------
    ForceFields::ForceField *constructForceField(ROMol &mol,
      MMFFMolProperties *mmffMolProperties, double nonBondedThresh,
      int confId, bool ignoreInterfragInteractions)
    {
      PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
      PRECONDITION(mmffMolProperties->isValid(), "missing atom types - invalid force-field");

      ForceFields::ForceField *res = new ForceFields::ForceField();
      // add the atomic positions:
      Conformer &conf = mol.getConformer(confId);
      for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        res->positions().push_back(&(conf.getAtomPos(i)));
      }
      
      if (mmffMolProperties->getMMFFBondTerm()) {
        Tools::addBonds(mol, mmffMolProperties, res);
      }
      if (mmffMolProperties->getMMFFAngleTerm()) {
        Tools::addAngles(mol, mmffMolProperties, res);
      }
      if (mmffMolProperties->getMMFFStretchBendTerm()) {
        Tools::addStretchBend(mol, mmffMolProperties, res);
      }
      if (mmffMolProperties->getMMFFOopTerm()) {
        Tools::addOop(mol, mmffMolProperties, res);
      }
      if (mmffMolProperties->getMMFFTorsionTerm()) {
        Tools::addTorsions(mol, mmffMolProperties, res);
      }
      if (mmffMolProperties->getMMFFVdWTerm()
        || mmffMolProperties->getMMFFEleTerm()) {
        boost::shared_array<boost::uint8_t> neighborMat = Tools::buildNeighborMatrix(mol);
        if (mmffMolProperties->getMMFFVdWTerm()) {
          Tools::addVdW(mol, confId, mmffMolProperties, res,
            neighborMat, nonBondedThresh, ignoreInterfragInteractions);
        }
        if (mmffMolProperties->getMMFFEleTerm()) {
          Tools::addEle(mol, confId, mmffMolProperties, res,
            neighborMat, nonBondedThresh, ignoreInterfragInteractions);
        }
      }

      return res;
    }
  }
}
