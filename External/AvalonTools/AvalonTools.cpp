// $Id$
//
// Created by Greg Landrum, July 2008
//

#include <DataStructs/ExplicitBitVect.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <Geometry/point.h>
#include "AvalonTools.h"

extern "C" {
#include "local.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "utilities.h"
#include "ssmatch.h"
#include "smi2mol.h"
#include "canonizer.h"
#include "layout.h"
#include "struchk.h"

extern int RunStruchk(struct reaccs_molecule_t **mpp,
                      struct data_line_t *data_list);
}

// already defined in struchk.c
// FILE *log_file=NULL;

namespace AvalonTools {
  using namespace RDKit;
  namespace {
    char *getFp(struct reaccs_molecule_t *molPtr,unsigned int bitFlags,
                bool isQuery,unsigned int nBytes){
      PRECONDITION(molPtr,"bad molecule");
      while(nBytes%4) ++nBytes;
      char *fingerprint = TypeAlloc(nBytes, char);
      SetFingerprintBits(molPtr,fingerprint,static_cast<int>(nBytes),
                         static_cast<int>(bitFlags),
                         static_cast<int>(isQuery),0);
      if(!isQuery){
        SetFingerprintBits(molPtr,fingerprint,static_cast<int>(nBytes),
                           static_cast<int>(bitFlags),
                           static_cast<int>(0),
                           ACCUMULATE_BITS|USE_DY_AROMATICITY);
      }
      return fingerprint;
    }
    void reaccsToFingerprint(struct reaccs_molecule_t *molPtr,std::vector<boost::uint32_t> &res,
                             unsigned int bitFlags=32767U,bool isQuery=false,bool resetVect=true,
                             unsigned int nBytes=64){
      res.clear();
      char *fingerprint=getFp(molPtr,bitFlags,isQuery,nBytes);
      for(unsigned int i=0;i<nBytes;i+=4){
        boost::uint32_t word;
        word = fingerprint[i] | (fingerprint[i+1]<<8) | (fingerprint[i+2]<<16) | (fingerprint[i+3]<<24);
        res.push_back(word);
      }

      MyFree(fingerprint);
    };
  
    void reaccsToFingerprint(struct reaccs_molecule_t *molPtr,ExplicitBitVect &res,
                             unsigned int bitFlags=32767U,bool isQuery=false,
                             bool resetVect=true,unsigned int nBytes=64){
      PRECONDITION(molPtr,"bad molecule");
      PRECONDITION(res.getNumBits()>=nBytes*8U,"res too small");
      if(resetVect) res.clearBits();

      char *fingerprint=getFp(molPtr,bitFlags,isQuery,nBytes);

      for(unsigned int i=0;i<nBytes;++i){
        char byte = fingerprint[i];
        if(byte){
          char mask=1;
          for (int j=0;j<8;++j){
            if(byte&mask){
              res.setBit(i*8+j);
            }
            mask = mask<<1;
          }
        }
      }
      MyFree(fingerprint);
    };
  
    struct reaccs_molecule_t *reaccsGetCoords(struct reaccs_molecule_t *molPtr){
      PRECONDITION(molPtr,"bad molecule");

      RecolorMolecule(molPtr);
      struct reaccs_molecule_t *res = LayoutMolecule(molPtr);
      POSTCONDITION(res,"could not layout molecule");
      return res;
    };

    struct reaccs_molecule_t *molToReaccs(const ROMol &mol){
      std::string molB=MolToMolBlock(mol,true);
      Utils::LocaleSwitcher ls;
      struct reaccs_molecule_t *res= MolStr2Mol((char *)molB.c_str());
      POSTCONDITION(res,"could not build a molecule");
      return res;
    }

    struct reaccs_molecule_t *stringToReaccs(const std::string &data,bool isSmiles){
      struct reaccs_molecule_t *res;
      if(isSmiles){
        res = SMIToMOL(data.c_str(),DY_AROMATICITY);
      } else {
        Utils::LocaleSwitcher ls;
        res= MolStr2Mol((char *)data.c_str());
      }
      if(!res){
        if(isSmiles){
          BOOST_LOG(rdErrorLog)<<"ERROR could not build molecule from smiles: "<<data<<std::endl;
        } else {
          BOOST_LOG(rdErrorLog)<<"ERROR could not build molecule from molblock: \n"<<data<<std::endl;
        }
      }
      return res;
    }

  }  // end of anonymous namespace

  std::string getCanonSmiles(ROMol &mol,int flags){
    if(flags==-1) flags=DB_STEREO | CENTER_STEREO;
    std::string rdSmi=MolToSmiles(mol,true);
    char *canSmiles = CanSmiles(const_cast<char *>(rdSmi.c_str()),flags);
    std::string res;
    if(canSmiles){
      res=canSmiles;
      MyFree(canSmiles);
    }else {
      BOOST_LOG(rdErrorLog)<<"ERROR: no smiles generated for molecule."<<std::endl;
    }
    return res;
  }

  void getAvalonFP(const ROMol &mol,ExplicitBitVect &res,
                   unsigned int nBits,
                   bool isQuery,
                   bool resetVect,
                   unsigned int bitFlags){
    if(nBits%8) {
      BOOST_LOG(rdWarningLog)<<"Warning: number of bits ("<<nBits<<") is not evenly divisible by 8. Rounding to the nearest byte."<<std::endl;
    }
    unsigned int nBytes = nBits/8;
    struct reaccs_molecule_t *mp=molToReaccs(mol);
    reaccsToFingerprint(mp,res,bitFlags,isQuery,resetVect,nBytes);
    FreeMolecule(mp);
  }
  void getAvalonFP(const ROMol &mol,std::vector<boost::uint32_t> &res,
                   unsigned int nBits,
                   bool isQuery,
                   bool resetVect,
                   unsigned int bitFlags){
    if(nBits%8) {
      BOOST_LOG(rdWarningLog)<<"Warning: number of bits ("<<nBits<<") is not evenly divisible by 8. Rounding to the nearest byte."<<std::endl;
    }
    unsigned int nBytes = nBits/8;
    struct reaccs_molecule_t *mp=molToReaccs(mol);
    reaccsToFingerprint(mp,res,bitFlags,isQuery,resetVect,nBytes);
    FreeMolecule(mp);
  }

  unsigned int set2DCoords(ROMol &mol,bool clearConfs){
    struct reaccs_molecule_t *mp=molToReaccs(mol);
    struct reaccs_molecule_t *mp2=reaccsGetCoords(mp);
    //std::cerr<<"----\n"<<MolToMolStr(mp2)<<"--------\n";

    TEST_ASSERT(mp2->n_atoms==mol.getNumAtoms());

    RDKit::Conformer *conf = new RDKit::Conformer(mol.getNumAtoms());
    conf->set3D(false);
    for(unsigned int i=0;i<mol.getNumAtoms();++i){
      RDGeom::Point3D loc(mp2->atom_array[i].x,mp2->atom_array[i].y,mp2->atom_array[i].z);
      conf->setAtomPos(i,loc);
    }

    unsigned int res;
    if (clearConfs) {
      mol.clearConformers();
      conf->setId(0);
      mol.addConformer(conf);
      res=0;
    } else {
      res=mol.addConformer(conf,true);
    }

    FreeMolecule(mp);
    FreeMolecule(mp2);

    return res;
  }
  std::string set2DCoords(const std::string &data,bool isSmiles){
    struct reaccs_molecule_t *mp=stringToReaccs(data,isSmiles);
    std::string res="";
    if(mp){
      struct reaccs_molecule_t *mp2=reaccsGetCoords(mp);
      FreeMolecule(mp);
      Utils::LocaleSwitcher ls;
      char *molB = MolToMolStr(mp2);
      res=molB;
      FreeMolecule(mp2);
      MyFree(molB);
    } 
    return res;
  }


  std::string getCanonSmiles(const std::string &data,bool isSmiles,int flags){
    if(flags==-1) flags=DB_STEREO | CENTER_STEREO;
    char *smiles=0,*canSmiles=0;
    if(!isSmiles){
      struct reaccs_molecule_t *mp=stringToReaccs(data,isSmiles);
      if(mp){
        smiles = MOLToSMI(mp,ISOMERIC_SMILES);
        FreeMolecule(mp);
        canSmiles = CanSmiles(smiles, flags);
        MyFree(smiles);
      }
    } else {
      canSmiles = CanSmiles((char *)data.c_str(), flags);
    }
    std::string res="";
    if(canSmiles){
      res=canSmiles;
      MyFree(canSmiles);
    } else {
      BOOST_LOG(rdErrorLog)<<"ERROR: no smiles generated for molecule."<<std::endl;
    }
    return res;
  }

  void getAvalonFP(const std::string &data,bool isSmiles,ExplicitBitVect &res,
                   unsigned int nBits,
                   bool isQuery,
                   bool resetVect,
                   unsigned int bitFlags){
    if(nBits%8) {
      BOOST_LOG(rdWarningLog)<<"Warning: number of bits ("<<nBits<<") is not evenly divisible by 8. Rounding to the nearest byte."<<std::endl;
    }
    unsigned int nBytes = nBits/8;
    struct reaccs_molecule_t *mp=stringToReaccs(data,isSmiles);
    if(mp){
      reaccsToFingerprint(mp,res,bitFlags,isQuery,resetVect,nBytes);
      FreeMolecule(mp);
    } else {
      BOOST_LOG(rdErrorLog)<<"ERROR: no fingeprint generated for molecule."<<std::endl;
    }
  }
  void getAvalonFP(const std::string &data,bool isSmiles,std::vector<boost::uint32_t> &res,
                   unsigned int nBits,
                   bool isQuery,
                   bool resetVect,
                   unsigned int bitFlags){
    if(nBits%8) {
      BOOST_LOG(rdWarningLog)<<"Warning: number of bits ("<<nBits<<") is not evenly divisible by 8. Rounding to the nearest byte."<<std::endl;
    }
    unsigned int nBytes = nBits/8;
    struct reaccs_molecule_t *mp=stringToReaccs(data,isSmiles);
    if(mp){
      reaccsToFingerprint(mp,res,bitFlags,isQuery,resetVect,nBytes);
      FreeMolecule(mp);
    } else {
      BOOST_LOG(rdErrorLog)<<"ERROR: no fingeprint generated for molecule."<<std::endl;
    }
  }

  int _checkMolWrapper(struct reaccs_molecule_t **mpp){
    if(!*mpp) return BAD_MOLECULE;
    int res;
    struct reaccs_molecule_t *tmp=*mpp;
    res = RunStruchk(mpp,NULL);
    if(*mpp != tmp) {
      FreeMolecule(tmp);
    }
    return res;
  }
  
  /**
   * Wrapper around struchk.CheckMol
   * The molecule to check is passed in as a string. isSmiles
   * should be set to TRUE if the molecule is encoded as SMILES,
   * to FALSE if the molecule is encoded as sn MDL CTAB.
   * mp is an output parameter - it will point to the checked
   * molecule upon successful checking. In case of errors, mp may be 0.
   **/
  int checkMolString(const std::string &data, const bool isSmiles,
		     struct reaccs_molecule_t **mp) {
    int errs = 0;
    if(isSmiles){
      *mp = SMIToMOL(data.c_str(),DY_AROMATICITY);
    } else {
      Utils::LocaleSwitcher ls;
      *mp= MolStr2Mol((char *)data.c_str());
    }
    if(*mp) {
      errs = _checkMolWrapper(mp);
    } else {
      errs = BAD_MOLECULE;
    }
    return errs;
  }

  int initCheckMol(const std::string &optString) {
    return InitCheckMol((char *) optString.c_str());
  }

  RDKit::ROMOL_SPTR checkMol(int &errs, RDKit::ROMol& inMol) {
    struct reaccs_molecule_t *mp;
    RDKit::ROMol *rMol = 0;
    mp = molToReaccs(inMol);
    errs = _checkMolWrapper(&mp);
    if(mp){
      Utils::LocaleSwitcher ls;
      char *molStr = MolToMolStr(mp);
      FreeMolecule(mp);
      if(molStr){
        rMol = MolBlockToMol(molStr);
        MyFree(molStr);
      }
    }
    return RDKit::ROMOL_SPTR(rMol);
  }

  RDKit::ROMOL_SPTR checkMol(int &errs, const std::string &data, const bool isSmiles) {
    struct reaccs_molecule_t *mp;
    errs = checkMolString(data, isSmiles, &mp);
    if(mp) {
      Utils::LocaleSwitcher ls;
      char *molStr = MolToMolStr(mp);
      RDKit::ROMol *rMol = MolBlockToMol(molStr);
      FreeMolecule(mp);
      MyFree(molStr);
      return RDKit::ROMOL_SPTR(rMol);
    } else {
      return RDKit::ROMOL_SPTR();
    }
  }

  std::pair<std::string,int> checkMolString(const std::string &data, bool isSmiles){
    struct reaccs_molecule_t *mp;
    int errs = checkMolString(data, isSmiles, &mp);
    std::string molStr;
    if(mp) {
      Utils::LocaleSwitcher ls;
      char *tmp=MolToMolStr(mp);
      molStr = std::string(tmp);
      FreeMolecule(mp);
      MyFree(tmp);
    } else {
      molStr="";
    }
    return std::make_pair(molStr,errs);
  }


  void closeCheckMolFiles() {
  	CloseOpenFiles();
  }

}
