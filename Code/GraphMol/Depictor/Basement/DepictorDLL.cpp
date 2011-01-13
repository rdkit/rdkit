//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>

#include "Depictor.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/RDLog.h>
#include <stdlib.h>
namespace RDKit {
#ifdef WIN32_DLLBUILD
  //*************************************************************************************
  //
  //  Uses the Depict DLL to generate a Mol file (on disk) for a
  //  SMILES.
  //
  //  ARGUMENTS:
  //    smi:     the SMILES string to be converted
  //    fName:   the name of the output file
  //    dllName: (OPTIONAL) the name of the actual DLL to be opened.  The path
  //             to the DLL is taken from the environment (see below).
  //
  //  RETURNS:
  //   1 on success, 0 otherwise
  //
  //  NOTES:
  //   - To use the DLL, it's essential that the COMBICHEM_ROOT and COMBICHEM_RELEASE
  //     environment variables be set.  If this isn't done, this whole process
  //     will fail.
  //   - The Depict DLL is not re-entrant, which means it can't be used from multiple
  //     threads at the same time.  We're going to try and minimize this problem by
  //     leaving it open for as short a period of time as possible.
  //     Still, opening the DLL may fail, in which case we'll return 0 for failure.
  //     You can always try again.
  //
  //*************************************************************************************
  int SmilesToMolFileDLL(std::string smi,std::string fName,
			 std::string dllName){
    //std::string dllName="c:/glandrum/dgen/combi/bin/depict32-0.dll";
    const char *root=getenv("COMBICHEM_ROOT");
    const char *release=getenv("COMBICHEM_RELEASE");
    if(!root || !release){
      BOOST_LOG(rdErrorLog) << "ERROR: COMBICHEM_ROOT and COMBICHEM_RELEASE must be set to use depictor\n";
      return 0;
    }
    std::string fullName=root;
    fullName += std::string("/")+release;
    fullName += std::string("/bin/")+dllName;

    Depictor_TwoArgFunc func1;
    HINSTANCE hDLL;
    hDLL = LoadLibrary(fullName.c_str());
    if(!hDLL){
      BOOST_LOG(rdErrorLog) << "ERROR: could not load depict32.dll\n" << std::endl;
      return 0;
    }
    func1 = (Depictor_TwoArgFunc)GetProcAddress(hDLL,
						   "SMILESSTRINGTOMOLFILE");
    if(!func1){
      FreeLibrary(hDLL);
      BOOST_LOG(rdErrorLog) << "ERROR: could not find SmilesToMolFile function in DLL\n" << std::endl;
      return 0;
    }
    func1(smi.c_str(),fName.c_str());

    FreeLibrary(hDLL);
    return 1;
  }

  //*************************************************************************************
  //
  //  Adds 2D coordinates to a molecule using the Depict.dll
  //
  //  ARGUMENTS:
  //    mol:          the molecule to be altered
  //    tempFilename: (OPTIONAL) the name of the temporary file 
  //
  //  RETURNS:
  //   1 on success, 0 otherwise
  //
  //  Here's the process by which this works (it's kind of contorted):
  //   1) convert the mol to SMILES
  //   2) use the DLL to convert the SMILES to a mol file (on disk)
  //   3) parse the mol file into a temporary molecule
  //   4) do a substructure match from the old molecule to the
  //      temp one (which may have a different atom numbering or additional
  //      atoms added).
  //   5) Update the positions of the atoms on the old molecule.
  //   6) Free the temp molecule.
  //
  //  NOTES:
  //   - *FIX:* at the moment we're not doing anything to clear up the
  //     temp file created in this process.  It'll always have the same
  //     name unless the user explicitly asks that we do something different.
  //   - To use the DLL, it's essential that the COMBICHEM_ROOT and COMBICHEM_RELEASE
  //     environment variables be set.  If this isn't done, this whole process
  //     will fail.
  //   - See the notes above about failures when opening the DLL.
  //
  //*************************************************************************************
  int Add2DCoordsToMolDLL(ROMol &mol,std::string tempFilename){
    std::string smi=MolToSmiles(mol,true);
    int tmp = SmilesToMolFileDLL(smi,tempFilename);
    int res = -1;
    if(tmp){
      // build another mol from that mol file:
      RWMol *tmpMol = MolFileToMol(tempFilename,false);
      // match it up with the starting mol:
      //  (We need to do this because the depict.dll conversion
      //   to a mol file may have added Hs)
      MatchVectType matchVect;
      bool hasMatch=SubstructMatch(tmpMol,&mol,matchVect);
      if(hasMatch){
        const Conformer &conf = tmpMol->getCoformer(0);
        Coformer *nconf = new Coformer(mol.getNumAtoms());
	for(MatchVectType::const_iterator mvi=matchVect.begin();
	    mvi!=matchVect.end();mvi++){
	  
          nconf->setAtomPos(conf.getAtomPos(mvi->first));
	}
        confId = (int)mol.addConformer(nconf, true);
      }
      delete tmpMol;
    }
    return res;
  }

#endif  
}
