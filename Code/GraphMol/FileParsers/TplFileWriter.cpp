// $Id$
//
//  Copyright (C) 2007-2008 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include <RDGeneral/RDLog.h>
#include <sstream>
#include <fstream>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

namespace RDKit{
  namespace TPLWriter {
    void writeAtom(const ROMol &mol,unsigned int atomId,
                   ROMol::ConstConformerIterator confIt,
                   std::ostringstream &dest,
                   std::string partialChargeProp){
      const Atom *atom=mol.getAtomWithIdx(atomId);
      dest << atomId+1;
      dest << " " << atom->getSymbol();
      dest << " " << atom->getFormalCharge();
      std::string propVal;
      if(atom->hasProp(partialChargeProp)){
        atom->getProp(partialChargeProp,propVal);
      } else {
        propVal = "0.0";
      }
      dest << " " << propVal;
      
      const RDGeom::Point3D &pos=(*confIt)->getAtomPos(atomId);
      dest  << " " << 100.*pos.x << " " << 100.*pos.y << " " << 100.*pos.z;

      ROMol::ADJ_ITER nbrIdx,endNbrs;
      boost::tie(nbrIdx,endNbrs)=mol.getAtomNeighbors(atom);
      dest << " " << (endNbrs-nbrIdx);
      while(nbrIdx!=endNbrs){
        dest << " " << (*nbrIdx+1);
        ++nbrIdx;
      }

      // FIX: get this right:
      dest << " " << "U";

      dest << std::endl;
    }

    void writeBond(const ROMol &mol,unsigned int bondId,
                   ROMol::ConstConformerIterator confIt,
                   std::ostringstream &dest){
      const Bond *bond=mol.getBondWithIdx(bondId);
      dest << bondId+1;
      std::string bondLabel;

      switch(bond->getBondType()){
      case Bond::SINGLE:
        if(bond->getIsAromatic()){
          bondLabel="1.5";
        } else {
          bondLabel="1.0";
        }
        break;
      case Bond::DOUBLE:
        if(bond->getIsAromatic()){
          bondLabel="1.5";
        } else {
          bondLabel="2.0";
        }
        break;
      case Bond::AROMATIC:
        bondLabel="1.5";break;
      case Bond::TRIPLE:
        bondLabel="3.0";break;
      default:
        BOOST_LOG(rdWarningLog)<<"TPL files only support single, double, aromatic, and triple bonds." << std::endl;
        BOOST_LOG(rdWarningLog)<<"Bond of with type " << bond->getBondType() << " written as single in output." << std::endl;
        bondLabel="1.0";
      }
      dest << " " << bondLabel;

      dest << " " << bond->getBeginAtomIdx()+1 << " " << bond->getEndAtomIdx()+1;

      // FIX: add these
      dest << " " << "0"<< " " << "0";

      
      dest << std::endl;
    }
  }

  
  std::string MolToTPLText(const ROMol &mol,std::string partialChargeProp,bool writeFirstConfTwice){
    if(!mol.getNumConformers()){
      BOOST_LOG(rdErrorLog)<<"Cannot write molecules with no conformers to TPL files\n";
      return "";
    }
    std::ostringstream res;
    std::string tempStr;
    res << "BioCAD format, all rights reserved"<<std::endl;
    res << "Output from RDKit"<<std::endl;
    if(!mol.hasProp("_Name")){
      BOOST_LOG(rdWarningLog)<<"Molecule has no name; arbitrary name assigned.\n";
      tempStr = "Unnamed molecule";
    } else {
      mol.getProp("_Name",tempStr);
    }
    res << "NAME "<<tempStr<<std::endl;
    res << "PROP 7 1"<<std::endl;
    res << mol.getNumAtoms() << " " << mol.getNumBonds() << std::endl;
    

    ROMol::ConstConformerIterator confIt=mol.beginConformers();
    // write the atoms:
    for(unsigned int i=0;i<mol.getNumAtoms();++i){
      TPLWriter::writeAtom(mol,i,confIt,res,partialChargeProp);
    }

    // write the bonds:
    for(unsigned int i=0;i<mol.getNumBonds();++i){
      TPLWriter::writeBond(mol,i,confIt,res);
    }

    // write the additional conformations:
    res << "CONFS "<< mol.getNumConformers()-1<<std::endl;
    if(!writeFirstConfTwice) ++confIt;
    while(confIt!=mol.endConformers()){
      std::stringstream tmpStrm;
      std::string confName;

      tmpStrm<<"conformer_"<<(*confIt)->getId();
      confName = tmpStrm.str();
      res << "NAME " << confName << std::endl;
      
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        const RDGeom::Point3D &pos=(*confIt)->getAtomPos(i);
        res << " " << 100.*pos.x << " " << 100.*pos.y << " " << 100.*pos.z << std::endl;
      }
      ++confIt;
      if(confIt!=mol.endConformers()){
        res << std::endl;
      }
    }

    return res.str();
  }


  void MolToTPLFile(const ROMol &mol,std::string fName,
                    std::string partialChargeProp,
                    bool writeFirstConfTwice){
    std::ofstream *outStream = new std::ofstream(fName.c_str());
    if(!outStream||!(*outStream)||outStream->bad()){
      std::ostringstream errout;
      errout << "Bad output file " << fName;
      throw BadFileException(errout.str());
    }

    std::string outString = MolToTPLText(mol,partialChargeProp,
                                         writeFirstConfTwice);
    *outStream  << outString;
    delete outStream;
  }    

}

