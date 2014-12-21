// $Id$
//
// Copyright (C)  2013 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//


#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/DistanceConstraint.h>
#include <ForceField/MMFF/AngleConstraint.h>
#include <ForceField/MMFF/TorsionConstraint.h>
#include <ForceField/MMFF/PositionConstraint.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <ForceField/ForceField.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "testMMFFForceField.h"


#define FCON_TOLERANCE   0.05
#define ENERGY_TOLERANCE   0.025


using namespace RDKit;

bool fexist(std::string filename)
{
  std::ifstream ifStream(filename.c_str());
  
  return (ifStream ? true : false);
}


bool fgrep(std::fstream &fStream, std::string key, std::string &line)
{
  bool found = false;
  std::streampos current = fStream.tellg();
  while ((!found) && (!(std::getline
    (fStream, line).rdstate() & std::ifstream::failbit))) {
    found = (line.find(key) != std::string::npos);
  }
  if (!found) {
    fStream.seekg(current);
  }
  
  return found;
}


bool fgrep(std::fstream &rdkFStream, std::string key)
{
  std::string line;
  
  return fgrep(rdkFStream, key, line);
}


bool getLineByNum(std::fstream &stream,
  std::streampos startPos, unsigned int n, std::string& line)
{
  std::streampos current = stream.tellg();
  stream.seekg(startPos);
  bool fail;
  unsigned int i = 0;
  while ((i <= n) && (!(fail = std::getline
    (stream, line).rdstate() & std::ifstream::failbit))) {
    ++i;
  }
  stream.seekg(current);
  
  return (!fail);
}
  

bool sortBondStretchInstanceVec(BondStretchInstance *a, BondStretchInstance *b)
{
  bool crit = false;
  
  if (a->iAtomType == b->iAtomType) {
    if (a->jAtomType == b->jAtomType) {
      crit = (a->ffType < b->ffType);
    }
    else {
      crit = (a->jAtomType < b->jAtomType);
    }
  }
  else {
    crit = (a->iAtomType < b->iAtomType);
  }
  
  return crit;
}


bool sortAngleBendInstanceVec(AngleBendInstance *a, AngleBendInstance *b)
{
  bool crit = false;
  
  if (a->jAtomType == b->jAtomType) {
    if (a->iAtomType == b->iAtomType) {
      if (a->kAtomType == b->kAtomType) {
        crit = (a->ffType < b->ffType);
      }
      else {
        crit = (a->kAtomType < b->kAtomType);
      }
    }
    else {
      crit = (a->iAtomType < b->iAtomType);
    }
  }
  else {
    crit = (a->jAtomType < b->jAtomType);
  }
  
  return crit;
}


bool sortStretchBendInstanceVec(StretchBendInstance *a, StretchBendInstance *b)
{
  bool crit = false;
  
  if (a->jAtomType == b->jAtomType) {
    if (a->iAtomType == b->iAtomType) {
      if (a->kAtomType == b->kAtomType) {
        crit = (a->ffType < b->ffType);
      }
      else {
        crit = (a->kAtomType < b->kAtomType);
      }
    }
    else {
      crit = (a->iAtomType < b->iAtomType);
    }
  }
  else {
    crit = (a->jAtomType < b->jAtomType);
  }
  
  return crit;
}


bool sortOopBendInstanceVec(OopBendInstance *a, OopBendInstance *b)
{
  bool crit = false;
  
  if (a->jAtomType == b->jAtomType) {
    if (a->iAtomType == b->iAtomType) {
      if (a->kAtomType == b->kAtomType) {
        crit = (a->lAtomType < b->lAtomType);
      }
      else {
        crit = (a->kAtomType < b->kAtomType);
      }
    }
    else {
      crit = (a->iAtomType < b->iAtomType);
    }
  }
  else {
    crit = (a->jAtomType < b->jAtomType);
  }
  
  return crit;
}


bool sortTorsionInstanceVec(TorsionInstance *a, TorsionInstance *b)
{
  bool crit = false;
  
  if (a->jAtomType == b->jAtomType) {
    if (a->kAtomType == b->kAtomType) {
      if (a->iAtomType == b->iAtomType) {
        if (a->lAtomType == b->lAtomType) {
          crit = (a->ffType < b->ffType);
        }
        else {
          crit = (a->lAtomType < b->lAtomType);
        }
      }
      else {
        crit = (a->iAtomType < b->iAtomType);
      }
    }
    else {
      crit = (a->kAtomType < b->kAtomType);
    }
  }
  else {
    crit = (a->jAtomType < b->jAtomType);
  }
  
  return crit;
}


void fixBondStretchInstance(BondStretchInstance *bondStretchInstance)
{
  if (bondStretchInstance->iAtomType > bondStretchInstance->jAtomType) {
    unsigned int temp;
    temp = bondStretchInstance->iAtomType;
    bondStretchInstance->iAtomType = bondStretchInstance->jAtomType;
    bondStretchInstance->jAtomType = temp;
  }
}


void fixAngleBendInstance(AngleBendInstance *angleBendInstance)
{
  if (angleBendInstance->iAtomType > angleBendInstance->kAtomType) {
    unsigned int temp;
    temp = angleBendInstance->iAtomType;
    angleBendInstance->iAtomType = angleBendInstance->kAtomType;
    angleBendInstance->kAtomType = temp;
  }
}


void fixOopBendInstance(OopBendInstance *oopBendInstance)
{
  if (oopBendInstance->iAtomType > oopBendInstance->kAtomType) {
    unsigned int temp;
    temp = oopBendInstance->iAtomType;
    oopBendInstance->iAtomType = oopBendInstance->kAtomType;
    oopBendInstance->kAtomType = temp;
  }
}


void fixTorsionInstance(TorsionInstance *torsionInstance)
{
  if (torsionInstance->jAtomType > torsionInstance->kAtomType) {
    unsigned int temp;
    temp = torsionInstance->jAtomType;
    torsionInstance->jAtomType = torsionInstance->kAtomType;
    torsionInstance->kAtomType = temp;
    temp = torsionInstance->iAtomType;
    torsionInstance->iAtomType = torsionInstance->lAtomType;
    torsionInstance->lAtomType = temp;
  }
  else if (torsionInstance->jAtomType == torsionInstance->kAtomType) {
    unsigned int temp;
    if (torsionInstance->iAtomType > torsionInstance->lAtomType) {
      temp = torsionInstance->iAtomType;
      torsionInstance->iAtomType = torsionInstance->lAtomType;
      torsionInstance->lAtomType = temp;
    }
  }
}


int mmffValidationSuite(int argc, char *argv[])
{
  std::string arg;
  std::string ffVariant = "";
  std::vector<std::string> ffVec;
  std::string verbosity;
  std::string molFile = "";
  std::string molType = "";
  std::vector<std::string> molFileVec;
  std::vector<std::string> molTypeVec;
  std::string rdk = "";
  std::string ref = "";
  std::streampos rdkPos;
  std::streampos rdkCurrent;
  std::streampos refPos;
  std::streampos refCurrent;
  unsigned int i = 1;
  unsigned int n = 0;
  bool error = false;
  bool fullTest = false;
  int inc = 4;
  bool testFailure = false;
  while (i < argc) {
    arg = argv[i];
    error = (arg.at(0) != '-');
    if (error) {
      break;
    }
    error = ((arg == "-h") || (arg == "--help"));
    if (error) {
      break;
    }
    else if (arg == "-f") {
      error = ((i + 1) >= argc);
      if (error) {
        break;
      }
      ffVariant = argv[i + 1];
      error = ((ffVariant != "MMFF94") && (ffVariant != "MMFF94s"));
      if (error) {
        break;
      }
      ++i;
    }
    else if ((arg == "-sdf") || (arg == "-smi")) {
      molType = arg.substr(1);
      if ((i + 1) < argc) {
        molFile = argv[i + 1];
        ++i;
      }
    }
    else if (arg == "-L") {
      fullTest = true;
      inc = 1;
    }
    else if (arg == "-l") {
      error = ((i + 1) >= argc);
      if (error) {
        break;
      }
      rdk = argv[i + 1];
      ++i;
    }
    ++i;
  }
  if (error) {
    std::cerr <<
      "Usage: testMMFFForceField [-L] [-f {MMFF94|MMFF94s}] "
      "[{-sdf [<sdf_file>] | -smi [<smiles_file>]}] [-l <log_file>]"
      << std::endl;
    
    return -1;
  }
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Official MMFF validation suite." << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/ForceField/MMFF/test_data/";
  if (molFile != "") {
    ffVariant = ((molFile.substr(0, 7) == "MMFF94s") ? "MMFF94s" : "MMFF94");
  }
  if (ffVariant == "") {
    ffVec.push_back("MMFF94");
    if (fullTest) {
      ffVec.push_back("MMFF94s");
    }
  }
  else {
    ffVec.push_back(ffVariant);
  }
  std::fstream rdkFStream;
  std::fstream refFStream;
  if (rdk == "") {
    rdk = pathName + "testMMFFForceField.log";
  }
  bool firstTest = true;
  if (molType == "") {
    molTypeVec.push_back("sdf");
    molTypeVec.push_back("smi");
  }
  else {
    molTypeVec.push_back(molType);
  }
  boost::logging::disable_logs("rdApp.error");
  for (std::vector<std::string>::iterator molTypeIt = molTypeVec.begin();
    molTypeIt != molTypeVec.end(); ++molTypeIt) {
    bool checkEnergy = (*molTypeIt == "sdf");
    for (std::vector<std::string>::iterator ffIt = ffVec.begin();
      ffIt != ffVec.end(); ++ffIt) {
      ref = pathName + (*ffIt) + "_reference.log";
      molFileVec.clear();
      if (molFile == "") {
        molFileVec.push_back(pathName + (*ffIt) + "_dative." + (*molTypeIt));
        molFileVec.push_back(pathName + (*ffIt) + "_hypervalent." + (*molTypeIt));
      }
      else {
        molFileVec.push_back(molFile);
      }
      for (std::vector<std::string>::iterator molFileIt = molFileVec.begin();
        molFileIt != molFileVec.end(); ++molFileIt) {
        std::vector<ROMol *> molVec;
        molVec.clear();
        if (*molTypeIt == "sdf") {
          SDMolSupplier sdfSupplier(*molFileIt, false, false);
          for (i = 0; i < sdfSupplier.length(); ++i) {
            molVec.push_back(sdfSupplier[i]);
          }
          sdfSupplier.reset();
        }
        else {
          SmilesMolSupplier smiSupplier(*molFileIt);
          for (i = 0; i < smiSupplier.length(); ++i) {
            molVec.push_back(smiSupplier[i]);
          }
          smiSupplier.reset();
        }
        SDWriter *sdfWriter = new SDWriter
          ((*molFileIt).substr(0, (*molFileIt).length() - 4) + "_min"
          + ((*molTypeIt == "smi") ? "_from_SMILES" : "") + ".sdf");
        ROMol *mol;
        std::string molName;
        std::vector<std::string> nameArray;

        rdkFStream.open(rdk.c_str(), (firstTest
          ? std::fstream::out | std::fstream::binary
          : (std::fstream::out | std::fstream::binary | std::fstream::app)));
        std::string computingKey = (*ffIt) + " energies for " + (*molFileIt);
        std::cerr << std::endl << "Computing " << computingKey << "..." << std::endl;
        for (i = 0; i < computingKey.length(); ++i) {
          rdkFStream << "*";
        }
        rdkFStream << std::endl << computingKey << std::endl;
        for (i = 0; i < computingKey.length(); ++i) {
          rdkFStream << "*";
        }
        rdkFStream << std::endl << std::endl;
        for (i = 0; i < molVec.size(); i += inc) {
          if (*molTypeIt == "sdf") {
            mol = molVec[i];
          }
          else {
            mol = MolOps::addHs(*(molVec[i]));
            DGeomHelpers::EmbedMolecule(*mol);
          }
          MMFF::sanitizeMMFFMol((RWMol &)(*mol));
          if (mol->hasProp("_Name")) {
            mol->getProp("_Name", molName);
            rdkFStream << molName << std::endl;
            nameArray.push_back(molName);
          }
          MMFF::MMFFMolProperties *mmffMolProperties = new
            MMFF::MMFFMolProperties(*mol, *ffIt, MMFF::MMFF_VERBOSITY_HIGH, rdkFStream);
          if ((!mmffMolProperties) || (!(mmffMolProperties->isValid()))) {
            std::cerr << molName + ": error setting up force-field" << std::endl;
            continue;
          }
          else {
            ForceFields::ForceField *field = MMFF::constructForceField
              (*mol, mmffMolProperties, 100.0, -1, false);
            field->initialize();
            field->minimize();
            sdfWriter->write(*mol);
            rdkFStream << "\nTOTAL MMFF ENERGY              ="
              << std::right << std::setw(16) << std::fixed
              << std::setprecision(4) << field->calcEnergy() << "\n\n" << std::endl;
          }
        }
        sdfWriter->close();
        rdkFStream.close();
        std::cerr << "Checking against " << ref << "..." << std::endl;
        rdkFStream.open(rdk.c_str(), std::fstream::in | std::fstream::binary);
        refFStream.open(ref.c_str(), std::fstream::in | std::fstream::binary);
        
        bool found = false;
        std::string skip;
        std::string errorMsg;
        std::string corruptedMsg;
        std::string energyString;
        std::vector<BondStretchInstance *> rdkBondStretchInstanceVec;
        std::vector<BondStretchInstance *> refBondStretchInstanceVec;
        std::vector<AngleBendInstance *> rdkAngleBendInstanceVec;
        std::vector<AngleBendInstance *> refAngleBendInstanceVec;
        std::vector<StretchBendInstance *> rdkStretchBendInstanceVec;
        std::vector<StretchBendInstance *> refStretchBendInstanceVec;
        std::vector<OopBendInstance *> rdkOopBendInstanceVec;
        std::vector<OopBendInstance *> refOopBendInstanceVec;
        std::vector<TorsionInstance *> rdkTorsionInstanceVec;
        std::vector<TorsionInstance *> refTorsionInstanceVec;
        double rdkEnergy;
        double rdkEleEnergy;
        double refEnergy;
        double refVdWEnergy;
        double refEleEnergy;
        unsigned int n_failed = 0;
        bool failed;
        rdkFStream.seekg(0);
        fgrep(rdkFStream, computingKey);
        for (i = 0; i < nameArray.size(); i += inc) {
          error = false;
          failed = false;
          errorMsg = rdk + ": Molecule not found";
          found = fgrep(rdkFStream, nameArray[i]);
          if (!found) {
            break;
          }
          errorMsg = ref + ": Molecule not found";
          found = fgrep(refFStream, nameArray[i]);
          if (!found) {
            break;
          }
          refPos = refFStream.tellg();
          found = false;
          bool refHaveBondStretch = false;
          bool refHaveAngleBend = false;
          bool refHaveStretchBend = false;
          bool refHaveOopBend = false;
          bool refHaveTorsion = false;
          bool refHaveVdW = false;
          bool refHaveEle = false;
          std::string rdkLine;
          std::string refLine;
          while ((!found) && (!(std::getline
            (refFStream, refLine).rdstate() & std::ifstream::failbit))) {
            std::size_t occ = refLine.find("****");
            found = (occ != std::string::npos);
            if (!found) {
              if (!refHaveVdW) {
                occ = refLine.find("Net vdW");
                refHaveVdW = (occ != std::string::npos);
                if (refHaveVdW) {
                  std::stringstream refVdWStringStream(refLine);
                  refVdWStringStream >> skip >> skip >> refVdWEnergy;
                }
              }
              if (!refHaveEle) {
                occ = refLine.find("Electrostatic");
                refHaveEle = (occ != std::string::npos);
                if (refHaveEle) {
                  std::stringstream refEleStringStream(refLine);
                  refEleStringStream >> skip >> refEleEnergy;
                }
              }
              if (!refHaveBondStretch) {
                occ = refLine.find("B O N D   S T R E T C H I N G");
                refHaveBondStretch = (occ != std::string::npos);
              }
              if (!refHaveAngleBend) {
                occ = refLine.find("A N G L E   B E N D I N G");
                refHaveAngleBend = (occ != std::string::npos);
              }
              if (!refHaveStretchBend) {
                occ = refLine.find("S T R E T C H   B E N D I N G");
                refHaveStretchBend = (occ != std::string::npos);
              }
              if (!refHaveOopBend) {
                occ = refLine.find("O U T - O F - P L A N E    B E N D I N G");
                refHaveOopBend = (occ != std::string::npos);
              }
              if (!refHaveTorsion) {
                occ = refLine.find("T O R S I O N A L");
                refHaveTorsion = (occ != std::string::npos);
              }
            }
          }
          refFStream.seekg(refPos);
          errorMsg = "";


          // --------------------------------------------------
          // B O N D   S T R E T C H I N G
          // --------------------------------------------------
          if (refHaveBondStretch) {
            corruptedMsg = ": Bond stretching: corrupted input data\n";
            found = fgrep(rdkFStream, "B O N D   S T R E T C H I N G");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = fgrep(rdkFStream, "----------------");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = false;
            rdkBondStretchInstanceVec.clear();
            rdkPos = rdkFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (rdkFStream, rdkLine).rdstate() & std::ifstream::failbit))) {
              found = (!(rdkLine.length()));
              if (found) {
                break;
              }
              std::stringstream rdkLineStream(rdkLine);
              BondStretchInstance *bondStretchInstance = new BondStretchInstance();
              bondStretchInstance->idx = n;
              rdkLineStream >> skip >> skip >> skip >> skip
                >> bondStretchInstance->iAtomType
                >> bondStretchInstance->jAtomType
                >> bondStretchInstance->ffType
                >> skip >> skip >> skip >> skip
                >> bondStretchInstance->kb;
              fixBondStretchInstance(bondStretchInstance);
              rdkBondStretchInstanceVec.push_back(bondStretchInstance);
              ++n;
            }
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::sort(rdkBondStretchInstanceVec.begin(),
              rdkBondStretchInstanceVec.end(), sortBondStretchInstanceVec);
            found = fgrep(rdkFStream, "TOTAL BOND STRETCH ENERGY", energyString);
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::stringstream rdkBondStretchStringStream(energyString);
            rdkBondStretchStringStream >> skip >> skip >> skip >> skip >> skip >> rdkEnergy;
            found = fgrep(refFStream, "B O N D   S T R E T C H I N G");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = fgrep(refFStream, "----------------");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = false;
            refBondStretchInstanceVec.clear();
            refPos = refFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (refFStream, refLine).rdstate() & std::ifstream::failbit))) {
              found = (!(refLine.length()));
              if (found) {
                break;
              }
              std::stringstream refLineStream(refLine);
              BondStretchInstance *bondStretchInstance = new BondStretchInstance();
              bondStretchInstance->idx = n;
              refLineStream >> skip >> skip >> skip >> skip
                >> bondStretchInstance->iAtomType
                >> bondStretchInstance->jAtomType
                >> bondStretchInstance->ffType
                >> skip >> skip >> skip >> skip
                >> bondStretchInstance->kb;
              fixBondStretchInstance(bondStretchInstance);
              refBondStretchInstanceVec.push_back(bondStretchInstance);
              ++n;
            }
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::sort(refBondStretchInstanceVec.begin(),
              refBondStretchInstanceVec.end(), sortBondStretchInstanceVec);
            found = fgrep(refFStream, "TOTAL BOND STRAIN ENERGY", energyString);
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::stringstream refBondStretchStringStream(energyString);
            refBondStretchStringStream >> skip >> skip >> skip >> skip >> skip >> refEnergy;
            error = (rdkBondStretchInstanceVec.size() != refBondStretchInstanceVec.size());
            if (error) {
              failed = true;
              std::stringstream diff;
              diff << "Bond stretching: expected " << refBondStretchInstanceVec.size()
                << " interactions, found " << rdkBondStretchInstanceVec.size() << "\n";
              errorMsg += diff.str();
            }
            for (unsigned j = 0; (!error) & (j < rdkBondStretchInstanceVec.size()); ++j) {
              error = ((rdkBondStretchInstanceVec[j]->iAtomType != refBondStretchInstanceVec[j]->iAtomType)
                || (rdkBondStretchInstanceVec[j]->jAtomType != refBondStretchInstanceVec[j]->jAtomType)
                || (rdkBondStretchInstanceVec[j]->ffType != refBondStretchInstanceVec[j]->ffType)
                || (fabs(rdkBondStretchInstanceVec[j]->kb
                - refBondStretchInstanceVec[j]->kb) > FCON_TOLERANCE));
              if (error) {
                failed = true;
                bool haveLine;
                haveLine = getLineByNum(rdkFStream, rdkPos, rdkBondStretchInstanceVec[j]->idx, rdkLine);
                if (haveLine) {
                  haveLine = getLineByNum(refFStream, refPos, refBondStretchInstanceVec[j]->idx, refLine);
                }
                std::stringstream diff;
                diff << "Bond stretching: found a difference\n";
                if (haveLine) {
                  diff << "Expected:\n" << rdkLine << "\nFound:\n" << refLine << "\n";
                }
                errorMsg += diff.str();
              }
              delete rdkBondStretchInstanceVec[j];
              delete refBondStretchInstanceVec[j];
            }
            error = (fabs(rdkEnergy - refEnergy) > ENERGY_TOLERANCE);
            if (error && checkEnergy) {
              failed = true;
              std::stringstream diff;
              diff << "Bond stretching: energies differ\n"
                "Expected " << std::fixed << std::setprecision(4) << refEnergy
                << ", found " << std::fixed << std::setprecision(4) << rdkEnergy << "\n";
              errorMsg += diff.str();
            }
          }

          // --------------------------------------------------
          // A N G L E   B E N D I N G
          // --------------------------------------------------
          if (refHaveAngleBend) {
            corruptedMsg = ": Angle bending: corrupted input data\n";
            found = fgrep(rdkFStream, "A N G L E   B E N D I N G");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = fgrep(rdkFStream, "----------------");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = false;
            rdkAngleBendInstanceVec.clear();
            rdkPos = rdkFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (rdkFStream, rdkLine).rdstate() & std::ifstream::failbit))) {
              found = (!(rdkLine.length()));
              if (found) {
                break;
              }
              std::stringstream rdkLineStream(rdkLine);
              AngleBendInstance *angleBendInstance = new AngleBendInstance();
              angleBendInstance->idx = n;
              rdkLineStream >> skip >> skip >> skip >> skip >> skip >> skip
                >> angleBendInstance->iAtomType
                >> angleBendInstance->jAtomType
                >> angleBendInstance->kAtomType
                >> angleBendInstance->ffType
                >> skip >> skip >> skip >> skip
                >> angleBendInstance->ka;
              fixAngleBendInstance(angleBendInstance);
              rdkAngleBendInstanceVec.push_back(angleBendInstance);
              ++n;
            }
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::sort(rdkAngleBendInstanceVec.begin(),
              rdkAngleBendInstanceVec.end(), sortAngleBendInstanceVec);
            found = fgrep(rdkFStream, "TOTAL ANGLE BEND ENERGY", energyString);
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::stringstream rdkAngleBendStringStream(energyString);
            rdkAngleBendStringStream >> skip >> skip >> skip >> skip >> skip >> rdkEnergy;
            found = fgrep(refFStream, "A N G L E   B E N D I N G");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = fgrep(refFStream, "----------------");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = false;
            refAngleBendInstanceVec.clear();
            refPos = refFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (refFStream, refLine).rdstate() & std::ifstream::failbit))) {
              found = (!(refLine.length()));
              if (found) {
                break;
              }
              std::stringstream refLineStream(refLine);
              AngleBendInstance *angleBendInstance = new AngleBendInstance();
              angleBendInstance->idx = n;
              refLineStream >> skip >> skip >> skip >> skip
                >> angleBendInstance->iAtomType
                >> angleBendInstance->jAtomType
                >> angleBendInstance->kAtomType
                >> angleBendInstance->ffType
                >> skip >> skip >> skip >> skip
                >> angleBendInstance->ka;
              fixAngleBendInstance(angleBendInstance);
              refAngleBendInstanceVec.push_back(angleBendInstance);
              ++n;
            }
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::sort(refAngleBendInstanceVec.begin(),
              refAngleBendInstanceVec.end(), sortAngleBendInstanceVec);
            found = fgrep(refFStream, "TOTAL ANGLE STRAIN ENERGY", energyString);
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::stringstream refAngleBendStringStream(energyString);
            refAngleBendStringStream >> skip >> skip >> skip >> skip >> skip >> refEnergy;
            error = (rdkAngleBendInstanceVec.size() != refAngleBendInstanceVec.size());
            if (error) {
              failed = true;
              std::stringstream diff;
              diff << "Angle bending: expected " << refAngleBendInstanceVec.size()
                << " interactions, found " << rdkAngleBendInstanceVec.size() << "\n";
              errorMsg += diff.str();
            }
            for (unsigned j = 0; (!error) & (j < rdkAngleBendInstanceVec.size()); ++j) {
              error = ((rdkAngleBendInstanceVec[j]->iAtomType != refAngleBendInstanceVec[j]->iAtomType)
                || (rdkAngleBendInstanceVec[j]->jAtomType != refAngleBendInstanceVec[j]->jAtomType)
                || (rdkAngleBendInstanceVec[j]->kAtomType != refAngleBendInstanceVec[j]->kAtomType)
                || (rdkAngleBendInstanceVec[j]->ffType != refAngleBendInstanceVec[j]->ffType)
                || (fabs(rdkAngleBendInstanceVec[j]->ka
                - refAngleBendInstanceVec[j]->ka) > FCON_TOLERANCE));
              if (error) {
                failed = true;
                bool haveLine;
                haveLine = getLineByNum(rdkFStream, rdkPos, rdkAngleBendInstanceVec[j]->idx, rdkLine);
                if (haveLine) {
                  haveLine = getLineByNum(refFStream, refPos, refAngleBendInstanceVec[j]->idx, refLine);
                }
                std::stringstream diff;
                diff << "Angle bending: found a difference\n";
                if (haveLine) {
                  diff << "Expected:\n" << rdkLine << "\nFound:\n" << refLine << "\n";
                }
                errorMsg += diff.str();
              }
              delete rdkAngleBendInstanceVec[j];
              delete refAngleBendInstanceVec[j];
            }
            error = (fabs(rdkEnergy - refEnergy) > ENERGY_TOLERANCE);
            if (error && checkEnergy) {
              failed = true;
              std::stringstream diff;
              diff << "Angle bending: energies differ\n"
                "Expected " << std::fixed << std::setprecision(4) << refEnergy
                << ", found " << std::fixed << std::setprecision(4) << rdkEnergy << "\n";
              errorMsg += diff.str();
            }
          }

          // --------------------------------------------------
          // S T R E T C H   B E N D I N G
          // --------------------------------------------------
          if (refHaveStretchBend) {
            corruptedMsg = ": Stretch-bending: corrupted input data\n";
            found = fgrep(rdkFStream, "S T R E T C H   B E N D I N G");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = fgrep(rdkFStream, "----------------");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = false;
            rdkStretchBendInstanceVec.clear();
            rdkPos = rdkFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (rdkFStream, rdkLine).rdstate() & std::ifstream::failbit))) {
              found = (!(rdkLine.length()));
              if (found) {
                break;
              }
              std::stringstream rdkLineStream(rdkLine);
              StretchBendInstance *stretchBendInstance = new StretchBendInstance();
              stretchBendInstance->idx = n;
              rdkLineStream >> skip >> skip >> skip >> skip >> skip >> skip
                >> stretchBendInstance->iAtomType
                >> stretchBendInstance->jAtomType
                >> stretchBendInstance->kAtomType
                >> stretchBendInstance->ffType
                >> skip >> skip >> skip >> skip >> skip
                >> stretchBendInstance->kba;
              rdkStretchBendInstanceVec.push_back(stretchBendInstance);
              ++n;
            }
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::sort(rdkStretchBendInstanceVec.begin(),
              rdkStretchBendInstanceVec.end(), sortStretchBendInstanceVec);
            found = fgrep(rdkFStream, "TOTAL STRETCH-BEND ENERGY", energyString);
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::stringstream rdkStretchBendStringStream(energyString);
            rdkStretchBendStringStream >> skip >> skip >> skip >> skip >> rdkEnergy;
            found = fgrep(refFStream, "S T R E T C H   B E N D I N G");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = fgrep(refFStream, "----------------");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = false;
            refStretchBendInstanceVec.clear();
            refPos = refFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (refFStream, refLine).rdstate() & std::ifstream::failbit))) {
              found = (!(refLine.length()));
              if (found) {
                break;
              }
              std::stringstream refLineStream(refLine);
              StretchBendInstance *stretchBendInstance = new StretchBendInstance();
              stretchBendInstance->idx = n;
              refLineStream >> skip >> skip >> skip >> skip
                >> stretchBendInstance->iAtomType
                >> stretchBendInstance->jAtomType
                >> stretchBendInstance->kAtomType
                >> stretchBendInstance->ffType
                >> skip >> skip >> skip >> skip
                >> stretchBendInstance->kba;
              refStretchBendInstanceVec.push_back(stretchBendInstance);
              ++n;
            }
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::sort(refStretchBendInstanceVec.begin(),
              refStretchBendInstanceVec.end(), sortStretchBendInstanceVec);
            found = fgrep(refFStream, "TOTAL STRETCH-BEND STRAIN ENERGY", energyString);
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::stringstream refStretchBendStringStream(energyString);
            refStretchBendStringStream >> skip >> skip >> skip >> skip >> skip >> refEnergy;
            error = (rdkStretchBendInstanceVec.size() != refStretchBendInstanceVec.size());
            if (error) {
              failed = true;
              std::stringstream diff;
              diff << "Stretch-bending: expected " << refStretchBendInstanceVec.size()
                << " interactions, found " << rdkStretchBendInstanceVec.size() << "\n";
              errorMsg += diff.str();
            }
            for (unsigned j = 0; (!error) & (j < rdkStretchBendInstanceVec.size()); ++j) {
              error = ((rdkStretchBendInstanceVec[j]->iAtomType != refStretchBendInstanceVec[j]->iAtomType)
                || (rdkStretchBendInstanceVec[j]->jAtomType != refStretchBendInstanceVec[j]->jAtomType)
                || (rdkStretchBendInstanceVec[j]->kAtomType != refStretchBendInstanceVec[j]->kAtomType)
                || (rdkStretchBendInstanceVec[j]->ffType != refStretchBendInstanceVec[j]->ffType)
                || (fabs(rdkStretchBendInstanceVec[j]->kba
                - refStretchBendInstanceVec[j]->kba) > FCON_TOLERANCE));
              if (error) {
                failed = true;
                bool haveLine;
                haveLine = getLineByNum(rdkFStream, rdkPos, rdkStretchBendInstanceVec[j]->idx, rdkLine);
                if (haveLine) {
                  haveLine = getLineByNum(refFStream, refPos, refStretchBendInstanceVec[j]->idx, refLine);
                }
                std::stringstream diff;
                diff << "Stretch-bending: found a difference\n";
                if (haveLine) {
                  diff << "Expected:\n" << rdkLine << "\nFound:\n" << refLine << "\n";
                }
                errorMsg += diff.str();
              }
              delete rdkStretchBendInstanceVec[j];
              delete refStretchBendInstanceVec[j];
            }
            error = (fabs(rdkEnergy - refEnergy) > ENERGY_TOLERANCE);
            if (error && checkEnergy) {
              failed = true;
              std::stringstream diff;
              diff << "Stretch-bending: energies differ\n"
                "Expected " << std::fixed << std::setprecision(4) << refEnergy
                << ", found " << std::fixed << std::setprecision(4) << rdkEnergy << "\n";
              errorMsg += diff.str();
            }
          }

          // --------------------------------------------------
          // O U T - O F - P L A N E   B E N D I N G
          // --------------------------------------------------
          if (refHaveOopBend) {
            corruptedMsg = ": Out-of-plane bending: corrupted input data\n";
            found = fgrep(rdkFStream, "O U T - O F - P L A N E   B E N D I N G");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = fgrep(rdkFStream, "----------------");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = false;
            rdkOopBendInstanceVec.clear();
            rdkPos = rdkFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (rdkFStream, rdkLine).rdstate() & std::ifstream::failbit))) {
              found = (!(rdkLine.length()));
              if (found) {
                break;
              }
              std::stringstream rdkLineStream(rdkLine);
              OopBendInstance *oopBendInstance = new OopBendInstance();
              oopBendInstance->idx = n;
              rdkLineStream >> skip >> skip >> skip
                >> skip >> skip >> skip >> skip >> skip
                >> oopBendInstance->iAtomType
                >> oopBendInstance->jAtomType
                >> oopBendInstance->kAtomType
                >> oopBendInstance->lAtomType
                >> skip >> skip
                >> oopBendInstance->koop;
              fixOopBendInstance(oopBendInstance);
              rdkOopBendInstanceVec.push_back(oopBendInstance);
              ++n;
            }
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::sort(rdkOopBendInstanceVec.begin(),
              rdkOopBendInstanceVec.end(), sortOopBendInstanceVec);
            found = fgrep(rdkFStream, "TOTAL OUT-OF-PLANE BEND ENERGY", energyString);
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::stringstream rdkOopBendStringStream(energyString);
            rdkOopBendStringStream >> skip >> skip >> skip >> skip >> skip >> rdkEnergy;
            found = fgrep(refFStream, "O U T - O F - P L A N E    B E N D I N G");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = fgrep(refFStream, "----------------");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = false;
            refOopBendInstanceVec.clear();
            refPos = refFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (refFStream, refLine).rdstate() & std::ifstream::failbit))) {
              found = (!(refLine.length()));
              if (found) {
                break;
              }
              std::stringstream refLineStream(refLine);
              OopBendInstance *oopBendInstance = new OopBendInstance();
              oopBendInstance->idx = n;
              refLineStream >> skip >> skip >> skip >> skip >> skip
                >> oopBendInstance->iAtomType
                >> oopBendInstance->jAtomType
                >> oopBendInstance->kAtomType
                >> oopBendInstance->lAtomType
                >> skip >> skip
                >> oopBendInstance->koop;
              fixOopBendInstance(oopBendInstance);
              refOopBendInstanceVec.push_back(oopBendInstance);
              ++n;
            }
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::sort(refOopBendInstanceVec.begin(),
              refOopBendInstanceVec.end(), sortOopBendInstanceVec);
            found = fgrep(refFStream, "TOTAL OUT-OF-PLANE STRAIN ENERGY", energyString);
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::stringstream refOopBendStringStream(energyString);
            refOopBendStringStream >> skip >> skip >> skip >> skip >> skip >> refEnergy;
            error = (rdkOopBendInstanceVec.size() != refOopBendInstanceVec.size());
            if (error) {
              failed = true;
              std::stringstream diff;
              diff << "Out-of-plane bending: expected " << refOopBendInstanceVec.size()
                << " interactions, found " << rdkOopBendInstanceVec.size() << "\n";
              errorMsg += diff.str();
            }
            for (unsigned j = 0; (!error) & (j < rdkOopBendInstanceVec.size()); ++j) {
              error = ((rdkOopBendInstanceVec[j]->iAtomType != refOopBendInstanceVec[j]->iAtomType)
                || (rdkOopBendInstanceVec[j]->jAtomType != refOopBendInstanceVec[j]->jAtomType)
                || (rdkOopBendInstanceVec[j]->kAtomType != refOopBendInstanceVec[j]->kAtomType)
                || (rdkOopBendInstanceVec[j]->lAtomType != refOopBendInstanceVec[j]->lAtomType)
                || (fabs(rdkOopBendInstanceVec[j]->koop
                - refOopBendInstanceVec[j]->koop) > FCON_TOLERANCE));
              if (error) {
                failed = true;
                bool haveLine;
                haveLine = getLineByNum(rdkFStream, rdkPos, rdkOopBendInstanceVec[j]->idx, rdkLine);
                if (haveLine) {
                  haveLine = getLineByNum(refFStream, refPos, refOopBendInstanceVec[j]->idx, refLine);
                }
                std::stringstream diff;
                diff << "Out-of-plane bending: found a difference\n";
                if (haveLine) {
                  diff << "Expected:\n" << rdkLine << "\nFound:\n" << refLine << "\n";
                }
                errorMsg += diff.str();
              }
              delete rdkOopBendInstanceVec[j];
              delete refOopBendInstanceVec[j];
            }
            error = (fabs(rdkEnergy - refEnergy) > ENERGY_TOLERANCE);
            if (error && checkEnergy) {
              failed = true;
              std::stringstream diff;
              diff << "Out-of-plane bending: energies differ\n"
                "Expected " << std::fixed << std::setprecision(4) << refEnergy
                << ", found " << std::fixed << std::setprecision(4) << rdkEnergy << "\n";
              errorMsg += diff.str();
            }
          }

          // --------------------------------------------------
          // T O R S I O N A L
          // --------------------------------------------------
          if (refHaveTorsion) {
            corruptedMsg = ": Torsional: corrupted input data\n";
            found = fgrep(rdkFStream, "T O R S I O N A L");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = fgrep(rdkFStream, "----------------");
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            found = false;
            rdkTorsionInstanceVec.clear();
            rdkPos = rdkFStream.tellg();
            n = 0;
            while ((!found) && (!(std::getline
              (rdkFStream, rdkLine).rdstate() & std::ifstream::failbit))) {
              found = (!(rdkLine.length()));
              if (found) {
                break;
              }
              std::stringstream rdkLineStream(rdkLine);
              TorsionInstance *torsionInstance = new TorsionInstance();
              torsionInstance->idx = n;
              rdkLineStream >> skip >> skip >> skip
                >> skip >> skip >> skip >> skip >> skip
                >> torsionInstance->iAtomType
                >> torsionInstance->jAtomType
                >> torsionInstance->kAtomType
                >> torsionInstance->lAtomType
                >> torsionInstance->ffType
                >> skip >> skip
                >> torsionInstance->V1
                >> torsionInstance->V2
                >> torsionInstance->V3;
              fixTorsionInstance(torsionInstance);
              rdkTorsionInstanceVec.push_back(torsionInstance);
              ++n;
            }
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::sort(rdkTorsionInstanceVec.begin(),
              rdkTorsionInstanceVec.end(), sortTorsionInstanceVec);
            found = fgrep(rdkFStream, "TOTAL TORSIONAL ENERGY", energyString);
            if (!found) {
              errorMsg += rdk + corruptedMsg;
              break;
            }
            std::stringstream rdkTorsionStringStream(energyString);
            rdkTorsionStringStream >> skip >> skip >> skip >> skip >> rdkEnergy;
            found = fgrep(refFStream, "T O R S I O N A L");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = fgrep(refFStream, "----------------");
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            found = false;
            refTorsionInstanceVec.clear();
            refPos = refFStream.tellg();
            n = 0;
            double refTorsionEnergy;
            while ((!found) && (!(std::getline
              (refFStream, refLine).rdstate() & std::ifstream::failbit))) {
              found = (!(refLine.length()));
              if (found) {
                break;
              }
              std::stringstream refLineStream(refLine);
              TorsionInstance *torsionInstance = new TorsionInstance();
              torsionInstance->idx = n;
              refLineStream >> skip >> skip
                >> skip >> skip >> skip >> skip
                >> torsionInstance->iAtomType
                >> torsionInstance->jAtomType
                >> torsionInstance->kAtomType
                >> torsionInstance->lAtomType
                >> torsionInstance->ffType
                >> skip >> refTorsionEnergy
                >> torsionInstance->V1
                >> torsionInstance->V2
                >> torsionInstance->V3;
              if (!(ForceFields::MMFF::isDoubleZero(torsionInstance->V1)
                && ForceFields::MMFF::isDoubleZero(torsionInstance->V2)
                && ForceFields::MMFF::isDoubleZero(torsionInstance->V3)
                && ForceFields::MMFF::isDoubleZero(refTorsionEnergy))) {
                fixTorsionInstance(torsionInstance);
                refTorsionInstanceVec.push_back(torsionInstance);
              }
              ++n;
            }
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::sort(refTorsionInstanceVec.begin(),
              refTorsionInstanceVec.end(), sortTorsionInstanceVec);
            found = fgrep(refFStream, "TOTAL TORSION STRAIN ENERGY", energyString);
            if (!found) {
              errorMsg += ref + corruptedMsg;
              break;
            }
            std::stringstream refTorsionStringStream(energyString);
            refTorsionStringStream >> skip >> skip >> skip >> skip >> skip >> refEnergy;
            error = (rdkTorsionInstanceVec.size() != refTorsionInstanceVec.size());
            if (error) {
              failed = true;
              std::stringstream diff;
              diff << "Torsional: expected " << refTorsionInstanceVec.size()
                << " interactions, found " << rdkTorsionInstanceVec.size() << "\n";
              errorMsg += diff.str();
            }
            for (unsigned j = 0; (!error) & (j < rdkTorsionInstanceVec.size()); ++j) {
              error = ((rdkTorsionInstanceVec[j]->iAtomType != refTorsionInstanceVec[j]->iAtomType)
                || (rdkTorsionInstanceVec[j]->jAtomType != refTorsionInstanceVec[j]->jAtomType)
                || (rdkTorsionInstanceVec[j]->kAtomType != refTorsionInstanceVec[j]->kAtomType)
                || (rdkTorsionInstanceVec[j]->lAtomType != refTorsionInstanceVec[j]->lAtomType)
                || (rdkTorsionInstanceVec[j]->ffType != refTorsionInstanceVec[j]->ffType)
                || (fabs(rdkTorsionInstanceVec[j]->V1 - refTorsionInstanceVec[j]->V1) > FCON_TOLERANCE)
                || (fabs(rdkTorsionInstanceVec[j]->V2 - refTorsionInstanceVec[j]->V2) > FCON_TOLERANCE)
                || (fabs(rdkTorsionInstanceVec[j]->V3 - refTorsionInstanceVec[j]->V3) > FCON_TOLERANCE));
              if (error) {
                failed = true;
                bool haveLine;
                haveLine = getLineByNum(rdkFStream, rdkPos, rdkTorsionInstanceVec[j]->idx, rdkLine);
                if (haveLine) {
                  haveLine = getLineByNum(refFStream, refPos, refTorsionInstanceVec[j]->idx, refLine);
                }
                std::stringstream diff;
                diff << "Torsional: found a difference\n";
                if (haveLine) {
                  diff << "Expected:\n" << rdkLine << "\nFound:\n" << refLine << "\n";
                }
                errorMsg += diff.str();
              }
              delete rdkTorsionInstanceVec[j];
              delete refTorsionInstanceVec[j];
            }
            error = (fabs(rdkEnergy - refEnergy) > ENERGY_TOLERANCE);
            if (error && checkEnergy) {
              failed = true;
              std::stringstream diff;
              diff << "Torsional: energies differ\n"
                "Expected " << std::fixed << std::setprecision(4) << refEnergy
                << ", found " << std::fixed << std::setprecision(4) << rdkEnergy << "\n";
              errorMsg += diff.str();
            }
          }

          
          // --------------------------------------------------
          // V A N   D E R   W A A L S
          // E L E C T R O S T A T I C
          // --------------------------------------------------
          corruptedMsg = ": Van der Waals: corrupted input data\n";
          found = fgrep(rdkFStream, "V A N   D E R   W A A L S");
          if (!found) {
            errorMsg += rdk + corruptedMsg;
            break;
          }
          found = fgrep(rdkFStream, "----------------");
          if (!found) {
            errorMsg += rdk + corruptedMsg;
            break;
          }
          found = fgrep(rdkFStream, "TOTAL VAN DER WAALS ENERGY", energyString);
          if (!found) {
            errorMsg += rdk + corruptedMsg;
            break;
          }
          if (!refHaveVdW) {
            errorMsg += ref + corruptedMsg;
            break;
          }
          std::stringstream rdkVdWStringStream(energyString);
          rdkVdWStringStream >> skip >> skip >> skip >> skip >> skip >> skip >> rdkEnergy;
          corruptedMsg = ": Electrostatic: corrupted input data";
          found = fgrep(rdkFStream, "E L E C T R O S T A T I C");
          if (!found) {
            errorMsg += rdk + corruptedMsg;
            break;
          }
          found = fgrep(rdkFStream, "----------------");
          if (!found) {
            errorMsg += rdk + corruptedMsg;
            break;
          }
          found = fgrep(rdkFStream, "TOTAL ELECTROSTATIC ENERGY", energyString);
          if (!found) {
            errorMsg += rdk + corruptedMsg;
            break;
          }
          if (!refHaveEle) {
            errorMsg += ref + corruptedMsg;
            break;
          }
          std::stringstream rdkEleStringStream(energyString);
          rdkEleStringStream >> skip >> skip >> skip >> skip >> rdkEleEnergy;
          error = (fabs(rdkEnergy - refVdWEnergy) > ENERGY_TOLERANCE);
          if (error && checkEnergy) {
            failed = true;
            std::stringstream diff;
            diff << "Van der Waals: energies differ\n"
              "Expected " << std::fixed << std::setprecision(4) << refVdWEnergy
              << ", found " << std::fixed << std::setprecision(4) << rdkEnergy << "\n";
            errorMsg += diff.str();
          }
          error = (fabs(rdkEleEnergy - refEleEnergy) > ENERGY_TOLERANCE);
          if (error && checkEnergy) {
            failed = true;
            std::stringstream diff;
            diff << "Electrostatic: energies differ\n"
              "Expected " << std::fixed << std::setprecision(4) << refEleEnergy
              << ", found " << std::fixed << std::setprecision(4) << rdkEnergy << "\n";
            errorMsg += diff.str();
          }
          if (failed) {
            std::cerr << nameArray[i] << "\n\n" << errorMsg << std::endl;
            ++n_failed;
          }
        }
        if (!found) {
          if (!failed) {
            std::cerr << nameArray[i] << "\n\n";
          }
          std::cerr << errorMsg << std::endl;
        }
        else {
          unsigned int n_passed = nameArray.size() - n_failed;
          std::cerr << (n_failed ? "" : "All ") << n_passed << " test"
            << ((n_passed == 1) ? "" : "s") << " passed" << std::endl;
          if (n_failed) {
            std::cerr << n_failed << " test"
              << ((n_failed == 1) ? "" : "s") << " failed" << std::endl;
          }
        }
        if (n_failed) {
          testFailure = true;
        }
        rdkFStream.close();
        refFStream.close();
        firstTest = false;
      }
    }
  }
   
  TEST_ASSERT(!testFailure);
  std::cerr << "  done" << std::endl;
  return 0;
}

void testMMFFAllConstraints(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for all MMFF constraint terms." << std::endl;

  std::string molBlock =
    "butane\n"
    "     RDKit          3D\n"
    "butane\n"
    " 17 16  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.4280    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.7913   -0.2660    0.9927 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.9040    1.3004   -0.3485 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.5407    2.0271    0.3782 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.5407    1.5664   -1.3411 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.3320    1.3004   -0.3485 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.6953    1.5162   -1.3532 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.8080    0.0192    0.0649 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.4447   -0.7431   -0.6243 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.4447   -0.1966    1.0697 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    4.8980    0.0192    0.0649 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.6954    2.0627    0.3408 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.7913   -0.7267   -0.7267 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3633    0.7267    0.7267 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3633   -0.9926    0.2660 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3633    0.2660   -0.9926 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "  1 15  1  0  0  0  0\n"
    "  1 16  1  0  0  0  0\n"
    "  1 17  1  0  0  0  0\n"
    "  2  3  1  0  0  0  0\n"
    "  2  4  1  0  0  0  0\n"
    "  2 14  1  0  0  0  0\n"
    "  4  5  1  0  0  0  0\n"
    "  4  6  1  0  0  0  0\n"
    "  4  7  1  0  0  0  0\n"
    "  7  8  1  0  0  0  0\n"
    "  7  9  1  0  0  0  0\n"
    "  7 13  1  0  0  0  0\n"
    "  9 10  1  0  0  0  0\n"
    "  9 11  1  0  0  0  0\n"
    "  9 12  1  0  0  0  0\n"
    "M  END\n";
  RDKit::RWMol *mol;
  ForceFields::ForceField *field;
  MMFF::MMFFMolProperties *mmffMolProperties;

  // distance constraints
  ForceFields::MMFF::DistanceConstraintContrib *dc;
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  TEST_ASSERT(field);
  field->initialize();
  dc = new ForceFields::MMFF::DistanceConstraintContrib(field, 1, 3, 2.0, 2.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(dc));
  field->minimize();
  TEST_ASSERT(MolTransforms::getBondLength(mol->getConformer(), 1, 3) > 1.99);
  delete field;
  delete mmffMolProperties;
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  field->initialize();
  dc = new ForceFields::MMFF::DistanceConstraintContrib(field, 1, 3, true, -0.2, 0.2, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(dc));
  field->minimize();
  TEST_ASSERT(MolTransforms::getBondLength(mol->getConformer(), 1, 3) > 1.79);
  delete field;
  delete mmffMolProperties;
  delete mol;
  
  // angle constraints
  ForceFields::MMFF::AngleConstraintContrib *ac;
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  TEST_ASSERT(field);
  field->initialize();
  ac = new ForceFields::MMFF::AngleConstraintContrib(field, 1, 3, 6, 90.0, 90.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(ac));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getAngleDeg(mol->getConformer(), 1, 3, 6) == 90);
  delete field;
  delete mmffMolProperties;
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  field->initialize();
  ac = new ForceFields::MMFF::AngleConstraintContrib(field, 1, 3, 6, true, -10.0, 10.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(ac));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getAngleDeg(mol->getConformer(), 1, 3, 6) == 100);
  delete field;
  delete mmffMolProperties;
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  field->initialize();
  ac = new ForceFields::MMFF::AngleConstraintContrib(field, 1, 3, 6, false, -10.0, 10.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(ac));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getAngleDeg(mol->getConformer(), 1, 3, 6) == 10);
  delete field;
  delete mmffMolProperties;
  delete mol;
  
  // torsion constraints
  ForceFields::MMFF::TorsionConstraintContrib *tc;
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  MolTransforms::setDihedralDeg(mol->getConformer(), 1, 3, 6, 8, 50.0);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  TEST_ASSERT(field);
  field->initialize();
  tc = new ForceFields::MMFF::TorsionConstraintContrib(field, 1, 3, 6, 8, 10.0, 20.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(tc));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getDihedralDeg(mol->getConformer(), 1, 3, 6, 8) == 20);
  delete field;
  delete mmffMolProperties;
  MolTransforms::setDihedralDeg(mol->getConformer(), 1, 3, 6, 8, -30.0);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  field->initialize();
  tc = new ForceFields::MMFF::TorsionConstraintContrib(field, 1, 3, 6, 8, true, -10.0, -8.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(tc));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getDihedralDeg(mol->getConformer(), 1, 3, 6, 8) == -40);
  delete field;
  delete mmffMolProperties;
  MolTransforms::setDihedralDeg(mol->getConformer(), 1, 3, 6, 8, 30.0);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  field->initialize();
  tc = new ForceFields::MMFF::TorsionConstraintContrib(field, 1, 3, 6, 8, false, -10.0, -8.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(tc));
  field->minimize(500);
  TEST_ASSERT((int)MolTransforms::getDihedralDeg(mol->getConformer(), 1, 3, 6, 8) == -10);
  delete field;
  delete mmffMolProperties;
  delete mol;
  
  // position constraints
  ForceFields::MMFF::PositionConstraintContrib *pc;
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
  TEST_ASSERT(field);
  field->initialize();
  RDGeom::Point3D p = mol->getConformer().getAtomPos(1);
  pc = new ForceFields::MMFF::PositionConstraintContrib(field, 1, 0.3, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(pc));
  field->minimize();
  RDGeom::Point3D q = mol->getConformer().getAtomPos(1);
  TEST_ASSERT((p - q).length() < 0.3);
  delete field;
  delete mmffMolProperties;
  delete mol;
  
  // fixed atoms
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  field = RDKit::MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  RDGeom::Point3D fp = mol->getConformer().getAtomPos(1);
  field->fixedPoints().push_back(1);
  field->minimize();
  RDGeom::Point3D fq = mol->getConformer().getAtomPos(1);
  TEST_ASSERT((fp - fq).length() < 0.01);
  delete field;
  delete mol;

  std::cerr << "  done" << std::endl;
}

int main(int argc, char *argv[])
{
  if (!mmffValidationSuite(argc, argv)) {
    testMMFFAllConstraints();
  }
  
  return 0;
}
