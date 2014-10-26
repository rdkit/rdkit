// $Id$
//
//  Copyright (C) 2003-2014 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include "rdmolops.h"
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <string>
#include <math.h>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/MonomerInfo.h> 
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/python_streambuf.h>

namespace python = boost::python;
using boost_adaptbx::python::streambuf;

namespace RDKit{
  python::tuple fragmentOnSomeBondsHelper(const ROMol &mol,python::object pyBondIndices,
                                          unsigned int nToBreak,
                                          bool addDummies,
                                          python::object pyDummyLabels,
                                          python::object pyBondTypes){
    std::vector<unsigned int> *bondIndices=pythonObjectToVect(pyBondIndices,mol.getNumBonds());
    std::vector< std::pair<unsigned int,unsigned int> > *dummyLabels=0;
    if(pyDummyLabels){
      unsigned int nVs=python::extract<unsigned int>(pyDummyLabels.attr("__len__")());
      dummyLabels = new std::vector<std::pair<unsigned int,unsigned int> >(nVs);
      for(unsigned int i=0;i<nVs;++i){
        unsigned int v1=python::extract<unsigned int>(pyDummyLabels[i][0]);
        unsigned int v2=python::extract<unsigned int>(pyDummyLabels[i][1]);
        (*dummyLabels)[i] = std::make_pair(v1,v2);
      }
    }
    std::vector< Bond::BondType > *bondTypes=0;
    if(pyBondTypes){
      unsigned int nVs=python::extract<unsigned int>(pyBondTypes.attr("__len__")());
      if(nVs!=bondIndices->size()) {
	throw_value_error("bondTypes shorter than bondIndices");
      }
      bondTypes = new std::vector< Bond::BondType >(nVs);
      for(unsigned int i=0;i<nVs;++i){
        (*bondTypes)[i] = python::extract< Bond::BondType >(pyBondTypes[i]);
      }
    }
    std::vector<ROMOL_SPTR> frags;
    MolFragmenter::fragmentOnSomeBonds(mol,*bondIndices,frags,nToBreak,addDummies,dummyLabels,bondTypes);
    python::list res;
    for(unsigned int i=0;i<frags.size();++i){
      res.append(frags[i]);
    }
    delete bondIndices;
    delete dummyLabels;
    delete bondTypes;
    return python::tuple(res);
  }

  ROMol *fragmentOnBondsHelper(const ROMol &mol,python::object pyBondIndices,
                               bool addDummies,
                               python::object pyDummyLabels,
                               python::object pyBondTypes){
    std::vector<unsigned int> *bondIndices=pythonObjectToVect(pyBondIndices,mol.getNumBonds());
    std::vector< std::pair<unsigned int,unsigned int> > *dummyLabels=0;
    if(pyDummyLabels){
      unsigned int nVs=python::extract<unsigned int>(pyDummyLabels.attr("__len__")());
      dummyLabels = new std::vector<std::pair<unsigned int,unsigned int> >(nVs);
      for(unsigned int i=0;i<nVs;++i){
        unsigned int v1=python::extract<unsigned int>(pyDummyLabels[i][0]);
        unsigned int v2=python::extract<unsigned int>(pyDummyLabels[i][1]);
        (*dummyLabels)[i] = std::make_pair(v1,v2);
      }
    }
    std::vector< Bond::BondType > *bondTypes=0;
    if(pyBondTypes){
      unsigned int nVs=python::extract<unsigned int>(pyBondTypes.attr("__len__")());
      if(nVs!=bondIndices->size()) {
	throw_value_error("bondTypes shorter than bondIndices");
      }
      bondTypes = new std::vector< Bond::BondType >(nVs);
      for(unsigned int i=0;i<nVs;++i){
        (*bondTypes)[i] = python::extract< Bond::BondType >(pyBondTypes[i]);
      }
    }
    
    ROMol *res=MolFragmenter::fragmentOnBonds(mol,*bondIndices,addDummies,dummyLabels,bondTypes);
    delete bondIndices;
    delete dummyLabels;
    delete bondTypes;
    return res;
  }

  ROMol *renumberAtomsHelper(const ROMol &mol,python::object &pyNewOrder){
    if(python::extract<unsigned int>(pyNewOrder.attr("__len__")())<mol.getNumAtoms()){
      throw_value_error("atomCounts shorter than the number of atoms");
    }
    std::vector<unsigned int> *newOrder=pythonObjectToVect(pyNewOrder,mol.getNumAtoms());
    ROMol *res = MolOps::renumberAtoms(mol,*newOrder);
    delete newOrder;
    return res;
  }

  namespace {
    std::string getResidue(const ROMol &m,const Atom *at){
      if(at->getMonomerInfo()->getMonomerType()!=AtomMonomerInfo::PDBRESIDUE) return "";
      return static_cast<const AtomPDBResidueInfo *>(at->getMonomerInfo())->getResidueName();
    }
    std::string getChainId(const ROMol &m,const Atom *at){
      if(at->getMonomerInfo()->getMonomerType()!=AtomMonomerInfo::PDBRESIDUE) return "";
      return static_cast<const AtomPDBResidueInfo *>(at->getMonomerInfo())->getChainId();
    }
  }
  python::dict splitMolByPDBResidues(const ROMol &mol,
                                     python::object pyWhiteList,
                                     bool negateList){

    std::vector<std::string> *whiteList=NULL;
    if(pyWhiteList){
      unsigned int nVs=python::extract<unsigned int>(pyWhiteList.attr("__len__")());
      whiteList=new std::vector<std::string>(nVs);
      for(unsigned int i=0;i<nVs;++i){
        (*whiteList)[i] = python::extract<std::string>(pyWhiteList[i]);
      }
    }
    std::map<std::string,boost::shared_ptr<ROMol> > res=MolOps::getMolFragsWithQuery(mol,getResidue,false,
                                                                                     whiteList,negateList);
    delete whiteList;

    python::dict pyres;
    for(std::map<std::string,boost::shared_ptr<ROMol> >::const_iterator iter=res.begin();
        iter!=res.end();++iter){
      pyres[iter->first]=iter->second;
    }
    return pyres;
  }
  python::dict splitMolByPDBChainId(const ROMol &mol,
                                     python::object pyWhiteList,
                                     bool negateList){
                                    
    std::vector<std::string> *whiteList=NULL;
    if(pyWhiteList){
      unsigned int nVs=python::extract<unsigned int>(pyWhiteList.attr("__len__")());
      whiteList=new std::vector<std::string>(nVs);
      for(unsigned int i=0;i<nVs;++i){
        (*whiteList)[i] = python::extract<std::string>(pyWhiteList[i]);
      }
    }
    std::map<std::string,boost::shared_ptr<ROMol> > res=MolOps::getMolFragsWithQuery(mol,getChainId,false,
                                                                                     whiteList,negateList);
    delete whiteList;

    python::dict pyres;
    for(std::map<std::string,boost::shared_ptr<ROMol> >::const_iterator iter=res.begin();
        iter!=res.end();++iter){
      pyres[iter->first]=iter->second;
    }
    return pyres;
  }
  
  python::dict parseQueryDefFileHelper(python::object &input,bool standardize,
                                       std::string delimiter,std::string comment,
                                       unsigned int nameColumn,unsigned int smartsColumn){
    python::extract<std::string> get_filename(input);
    std::map<std::string,ROMOL_SPTR> queryDefs;

    if (get_filename.check()) {
        parseQueryDefFile(get_filename(),queryDefs,standardize,delimiter,comment,nameColumn,smartsColumn);
    } else {
        streambuf *sb=new streambuf(input);
        std::istream *istr=new streambuf::istream(*sb);
        parseQueryDefFile(istr,queryDefs,standardize,delimiter,comment,nameColumn,smartsColumn);
        delete istr;
        delete sb;
    }

    python::dict res;
    for(std::map<std::string,ROMOL_SPTR>::const_iterator iter=queryDefs.begin();iter!=queryDefs.end();++iter){
      res[iter->first]=iter->second;
    }
    
    return res;
  }
                                       
  
  void addRecursiveQueriesHelper(ROMol &mol,python::dict replDict,std::string propName){
    std::map<std::string,ROMOL_SPTR> replacements;
    for(unsigned int i=0;i<python::extract<unsigned int>(replDict.keys().attr("__len__")());++i){
      ROMol *m=python::extract<ROMol *>(replDict.values()[i]);
      ROMOL_SPTR nm(new ROMol(*m));
      std::string k=python::extract<std::string>(replDict.keys()[i]);
      replacements[k]=nm;
    }
    addRecursiveQueries(mol,replacements,propName);
    
  }

  ROMol *addHs(const ROMol &orig,bool explicitOnly=false,bool addCoords=false){
    return MolOps::addHs(orig,explicitOnly,addCoords);
  }
  ROMol *removeHs(const ROMol &orig,bool implicitOnly=false){
    return MolOps::removeHs(orig,implicitOnly);
  }
  int getSSSR(ROMol &mol) {
    VECT_INT_VECT rings;
    int nr = MolOps::findSSSR(mol, rings);
    return nr;
  }

  PyObject* replaceSubstructures(const ROMol &orig,
                                 const ROMol &query,
                                 const ROMol &replacement,
                                 bool replaceAll=false,
                                 unsigned int replacementConnectionPoint=0) {
    std::vector<ROMOL_SPTR> v=replaceSubstructs(orig, query,
                                                replacement, replaceAll,
                                                replacementConnectionPoint);
    PyObject *res=PyTuple_New(v.size());
    for(unsigned int i=0;i<v.size();++i){
      PyTuple_SetItem(res,i,
                      python::converter::shared_ptr_to_python(v[i]));
    }
    return res;
  }

  void addRecursiveQuery(ROMol &mol,
                         const ROMol &query,
                         unsigned int atomIdx,
                         bool preserveExistingQuery){
    if(atomIdx>=mol.getNumAtoms()){
      throw_value_error("atom index exceeds mol.GetNumAtoms()");
    }
    RecursiveStructureQuery *q = new RecursiveStructureQuery(new ROMol(query));

    Atom *oAt=mol.getAtomWithIdx(atomIdx);
    if(!oAt->hasQuery()){
      QueryAtom qAt(*oAt);
      static_cast<RWMol &>(mol).replaceAtom(atomIdx,&qAt);
      oAt = mol.getAtomWithIdx(atomIdx);
    }

    
    if(!preserveExistingQuery){
      delete oAt->getQuery();
      oAt->setQuery(q);
    } else {
      oAt->expandQuery(q,Queries::COMPOSITE_AND);
    }
    
  }
#ifdef RDK_32BIT_BUILD
  MolOps::SanitizeFlags sanitizeMol(ROMol &mol,int sanitizeOps,
                                    bool catchErrors) {
#else
  MolOps::SanitizeFlags sanitizeMol(ROMol &mol,unsigned int sanitizeOps,
                                    bool catchErrors) {
#endif    
    RWMol &wmol = static_cast<RWMol &>(mol);
    unsigned int operationThatFailed;
    if(catchErrors){
      try{
        MolOps::sanitizeMol(wmol,operationThatFailed,sanitizeOps);
      } catch (...){
      }
    } else {
      MolOps::sanitizeMol(wmol,operationThatFailed,sanitizeOps);
    }
    return static_cast<MolOps::SanitizeFlags>(operationThatFailed);
  }

  RWMol *getEditable(const ROMol &mol) {
    RWMol *res = static_cast<RWMol *>(new ROMol(mol,false));
    return res;
  }

  ROMol *getNormal(const RWMol &mol) {
    ROMol *res = static_cast<ROMol *>(new RWMol(mol));
    return res;
  }

  void kekulizeMol(ROMol &mol,bool clearAromaticFlags=false) {
    RWMol &wmol = static_cast<RWMol &>(mol);
    MolOps::Kekulize(wmol,clearAromaticFlags);
  }

  void cleanupMol(ROMol &mol){
    RWMol &rwmol = static_cast<RWMol &>(mol);
    MolOps::cleanUp(rwmol);
  }

  void setAromaticityMol(ROMol &mol){
    RWMol &wmol = static_cast<RWMol &>(mol);
    MolOps::setAromaticity(wmol);
  }

  void setConjugationMol(ROMol &mol) {
    RWMol &wmol = static_cast<RWMol &>(mol);
    MolOps::setConjugation(wmol);
  }

  void assignRadicalsMol(ROMol &mol) {
    RWMol &wmol = static_cast<RWMol &>(mol);
    MolOps::assignRadicals(wmol);
  }

  void setHybridizationMol(ROMol &mol) {
    RWMol &wmol = static_cast<RWMol &>(mol);
    MolOps::setHybridization(wmol);
  }
  
  VECT_INT_VECT getSymmSSSR(ROMol &mol) {
    VECT_INT_VECT rings;
    MolOps::symmetrizeSSSR(mol, rings);
    return rings;
  }
  PyObject *getDistanceMatrix(ROMol &mol, bool useBO=false,
                              bool useAtomWts=false,bool force=false,
                              const char *prefix=0) {
    int nats = mol.getNumAtoms();
    npy_intp dims[2];
    dims[0] = nats;
    dims[1] = nats;
    double *distMat;
    
    distMat = MolOps::getDistanceMat(mol, useBO, useAtomWts,force,prefix);
    
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    
    memcpy(static_cast<void *>(res->data),
         static_cast<void *>(distMat),nats*nats*sizeof(double));
    
    return PyArray_Return(res);
  }
  PyObject *get3DDistanceMatrix(ROMol &mol, int confId=-1,
                              bool useAtomWts=false,bool force=false,
                              const char *prefix=0) {
    int nats = mol.getNumAtoms();
    npy_intp dims[2];
    dims[0] = nats;
    dims[1] = nats;
    double *distMat;
    
    distMat = MolOps::get3DDistanceMat(mol, confId, useAtomWts,force,prefix);
    
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    
    memcpy(static_cast<void *>(res->data),
         static_cast<void *>(distMat),nats*nats*sizeof(double));
    
    return PyArray_Return(res);
  }

  PyObject *getAdjacencyMatrix(ROMol &mol, bool useBO=false,
                               int emptyVal=0,bool force=false,
                               const char *prefix=0) {
    int nats = mol.getNumAtoms();
    npy_intp  dims[2];
    dims[0] = nats;
    dims[1] = nats;

    double *tmpMat = MolOps::getAdjacencyMatrix(mol, useBO, emptyVal,force,prefix);
    
    PyArrayObject *res;
    if(useBO){
      // if we're using valence, the results matrix is made up of doubles
      res = (PyArrayObject *)PyArray_SimpleNew(2,dims,
                                              NPY_DOUBLE);
      memcpy(static_cast<void *>(res->data),
             static_cast<void *>(tmpMat),nats*nats*sizeof(double));
    } else {
      res = (PyArrayObject *)PyArray_SimpleNew(2,dims,
                                              NPY_INT);
      int *data = (int *)res->data;
      for(int i=0;i<nats;i++){
        for(int j=0;j<nats;j++){
          data[i*nats+j] = (int)round(tmpMat[i*nats+j]);
        }
      }
    }
    return PyArray_Return(res);
  }

  python::tuple GetMolFrags(const ROMol &mol,bool asMols,bool sanitizeFrags){
    python::list res;

    if(!asMols){
      VECT_INT_VECT frags;
      MolOps::getMolFrags(mol,frags);

      for(unsigned int i=0;i<frags.size();++i){
        python::list tpl;
        for(unsigned int j=0;j<frags[i].size();++j){
          tpl.append(frags[i][j]);
        }
        res.append(python::tuple(tpl));
      }
    } else {
      std::vector<boost::shared_ptr<ROMol> > frags;
      frags=MolOps::getMolFrags(mol,sanitizeFrags);
      for(unsigned int i=0;i<frags.size();++i){
        res.append(frags[i]);
      }
    }
    return python::tuple(res);
  }

  ExplicitBitVect *wrapLayeredFingerprint(const ROMol &mol,unsigned int layerFlags,
                                          unsigned int minPath,unsigned int maxPath,
                                          unsigned int fpSize,
                                          python::list atomCounts,
                                          ExplicitBitVect *includeOnlyBits,
                                          bool branchedPaths,
                                          python::object fromAtoms){
    std::vector<unsigned int> *lFromAtoms=pythonObjectToVect(fromAtoms,mol.getNumAtoms());
    std::vector<unsigned int> *atomCountsV=0;
    if(atomCounts){
      atomCountsV = new std::vector<unsigned int>;
      unsigned int nAts=python::extract<unsigned int>(atomCounts.attr("__len__")());
      if(nAts<mol.getNumAtoms()){
        throw_value_error("atomCounts shorter than the number of atoms");
      }
      atomCountsV->resize(nAts);
      for(unsigned int i=0;i<nAts;++i){
        (*atomCountsV)[i] = python::extract<unsigned int>(atomCounts[i]);
      }
    }

    ExplicitBitVect *res;
    res = RDKit::LayeredFingerprintMol(mol,layerFlags,minPath,maxPath,fpSize,atomCountsV,includeOnlyBits,branchedPaths,
                                       lFromAtoms);

    if(atomCountsV){
      for(unsigned int i=0;i<atomCountsV->size();++i){
        atomCounts[i] = (*atomCountsV)[i];
      }
      delete atomCountsV;
    }
    delete lFromAtoms;

    return res;
  }
  ExplicitBitVect *wrapPatternFingerprint(const ROMol &mol,
                                          unsigned int fpSize,
                                          python::list atomCounts,
                                          ExplicitBitVect *includeOnlyBits){
    std::vector<unsigned int> *atomCountsV=0;
    if(atomCounts){
      atomCountsV = new std::vector<unsigned int>;
      unsigned int nAts=python::extract<unsigned int>(atomCounts.attr("__len__")());
      if(nAts<mol.getNumAtoms()){
        throw_value_error("atomCounts shorter than the number of atoms");
      }
      atomCountsV->resize(nAts);
      for(unsigned int i=0;i<nAts;++i){
        (*atomCountsV)[i] = python::extract<unsigned int>(atomCounts[i]);
      }
    }

    ExplicitBitVect *res;
    res = RDKit::PatternFingerprintMol(mol,fpSize,
                                       atomCountsV,includeOnlyBits);

    if(atomCountsV){
      for(unsigned int i=0;i<atomCountsV->size();++i){
        atomCounts[i] = (*atomCountsV)[i];
      }
      delete atomCountsV;
    }
    
    return res;
  }


  ExplicitBitVect *wrapRDKFingerprintMol(const ROMol &mol,
                                         unsigned int minPath,
                                         unsigned int maxPath,
                                         unsigned int fpSize,
                                         unsigned int nBitsPerHash,
                                         bool useHs,
                                         double tgtDensity,
                                         unsigned int minSize,
                                         bool branchedPaths,
                                         bool useBondOrder,
                                         python::object atomInvariants,
                                         python::object fromAtoms,
                                         python::object atomBits
                                         ){
    std::vector<unsigned int> *lAtomInvariants=pythonObjectToVect<unsigned int>(atomInvariants);
    std::vector<unsigned int> *lFromAtoms=pythonObjectToVect(fromAtoms,mol.getNumAtoms());
    std::vector<std::vector<boost::uint32_t> > *lAtomBits=0;
    //if(!(atomBits.is_none())){
    if(atomBits!=python::object()){
      lAtomBits = new std::vector<std::vector<boost::uint32_t> >(mol.getNumAtoms());
    }
    ExplicitBitVect *res;
    res = RDKit::RDKFingerprintMol(mol,minPath,maxPath,fpSize,nBitsPerHash,
                                   useHs,tgtDensity,minSize,branchedPaths,
                                   useBondOrder,lAtomInvariants,lFromAtoms,lAtomBits);

    delete lAtomInvariants;
    delete lFromAtoms;

    if(lAtomBits){
      python::list &pyl=static_cast<python::list &>(atomBits);
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        python::list tmp;
        BOOST_FOREACH(boost::uint32_t v,(*lAtomBits)[i]){
          tmp.append(v);
        }
        pyl.append(tmp);
      }
      delete lAtomBits;
    }
    
    return res;
  }


  python::object findAllSubgraphsOfLengthsMtoNHelper(const ROMol &mol, unsigned int lowerLen,
                                                     unsigned int upperLen, bool useHs=false,
                                                     int rootedAtAtom=-1){
    if(lowerLen>upperLen){
      throw_value_error("lowerLen > upperLen");
    }
    
    INT_PATH_LIST_MAP oMap=findAllSubgraphsOfLengthsMtoN(mol,lowerLen,upperLen,useHs,rootedAtAtom);
    python::list res;
    for(unsigned int i=lowerLen;i<=upperLen;++i){
      python::list tmp;
      const PATH_LIST &pth=oMap[i];
      for(PATH_LIST_CI pthit=pth.begin();pthit!=pth.end();++pthit){
        tmp.append(python::tuple(*pthit));
      }
      res.append(tmp);
    }
    return python::tuple(res);
  };

  ROMol *pathToSubmolHelper(const ROMol &mol, python::object &path, 
                            bool useQuery,python::object atomMap){
    ROMol *result;
    PATH_TYPE pth;
    for(unsigned int i=0;i<python::extract<unsigned int>(path.attr("__len__")());++i){
      pth.push_back(python::extract<unsigned int>(path[i]));
    }
    std::map<int,int> mapping;
    result = Subgraphs::pathToSubmol(mol,pth,useQuery,mapping);
    if(atomMap!=python::object()){
      // make sure the optional argument actually was a dictionary
      python::dict typecheck=python::extract<python::dict>(atomMap);
      atomMap.attr("clear")();
      for(std::map<int,int>::const_iterator mIt=mapping.begin();
          mIt!=mapping.end();++mIt){
        atomMap[mIt->first]=mIt->second;
      }
    }
    return result;
  }

  struct molops_wrapper {
    static void wrap() {
      std::string docString;
      python::enum_<MolOps::SanitizeFlags>("SanitizeFlags")
        .value("SANITIZE_NONE",MolOps::SANITIZE_NONE)
        .value("SANITIZE_CLEANUP",MolOps::SANITIZE_CLEANUP)
        .value("SANITIZE_PROPERTIES",MolOps::SANITIZE_PROPERTIES)
        .value("SANITIZE_SYMMRINGS",MolOps::SANITIZE_SYMMRINGS)
        .value("SANITIZE_KEKULIZE",MolOps::SANITIZE_KEKULIZE)
        .value("SANITIZE_FINDRADICALS",MolOps::SANITIZE_FINDRADICALS)
        .value("SANITIZE_SETAROMATICITY",MolOps::SANITIZE_SETAROMATICITY)
        .value("SANITIZE_SETCONJUGATION",MolOps::SANITIZE_SETCONJUGATION)
        .value("SANITIZE_SETHYBRIDIZATION",MolOps::SANITIZE_SETHYBRIDIZATION)
        .value("SANITIZE_CLEANUPCHIRALITY",MolOps::SANITIZE_CLEANUPCHIRALITY)
        .value("SANITIZE_ADJUSTHS",MolOps::SANITIZE_ADJUSTHS)
        .value("SANITIZE_ALL",MolOps::SANITIZE_ALL)
        ;

      // ------------------------------------------------------------------------
      docString="Kekulize, check valencies, set aromaticity, conjugation and hybridization\n\
\n\
    - The molecule is modified in place.\n\
\n\
    - If sanitization fails, an exception will be thrown unless catchErrors is set\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
    - sanitizeOps: (optional) sanitization operations to be carried out\n\
                   these should be constructed by or'ing together the\n\
                   operations in rdkit.Chem.SanitizeFlags\n\
    - catchErrors: (optional) if provided, instead of raising an exception\n\
                   when sanitization fails (the default behavior), the \n\
                   first operation that failed (as defined in rdkit.Chem.SanitizeFlags)\n\
                   is returned. Zero is returned on success.\n\
\n";
      python::def("SanitizeMol", sanitizeMol,
                  (python::arg("mol"),
                   python::arg("sanitizeOps")=MolOps::SANITIZE_ALL,
                   python::arg("catchErrors")=false),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Get the smallest set of simple rings for a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
\n\
  RETURNS: a sequence of sequences containing the rings found as atom ids\n\
         The length of this will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.\n\
\n";
      python::def("GetSSSR", getSSSR, 
                  docString.c_str());
      
      // ------------------------------------------------------------------------
      docString="Get a symmetrized SSSR for a molecule.\n\
\n\
  The symmetrized SSSR is at least as large as the SSSR for a molecule.\n\
  In certain highly-symmetric cases (e.g. cubane), the symmetrized SSSR can be\n\
  a bit larger (i.e. the number of symmetrized rings is >= NumBonds-NumAtoms+1).\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
\n\
  RETURNS: a sequence of sequences containing the rings found as atom ids\n\
\n";
      python::def("GetSymmSSSR", getSymmSSSR,
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Does a non-SSSR ring finding for a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
\n\
  RETURNS: Nothing\n\
\n";
      python::def("FastFindRings", MolOps::fastFindRings, 
                  docString.c_str());

      
      // ------------------------------------------------------------------------
      docString="Adds hydrogens to the graph of a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - explicitOnly: (optional) if this toggle is set, only explicit Hs will\n\
      be added to the molecule.  Default value is 0 (add implicit and explicit Hs).\n\
\n\
    - addCoords: (optional) if this toggle is set, The Hs will have 3D coordinates\n\
      set.  Default value is 0 (no 3D coords).\n\
\n\
  RETURNS: a new molecule with added Hs\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
    - Much of the code assumes that Hs are not included in the molecular\n\
      topology, so be *very* careful with the molecule that comes back from\n\
      this function.\n\
\n";
      python::def("AddHs", addHs,
                  (python::arg("mol"),python::arg("explicitOnly")=false,
                   python::arg("addCoords")=false),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Removes any hydrogens from the graph of a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - implicitOnly: (optional) if this toggle is set, only implicit Hs will\n\
      be removed from the graph.  Default value is 0 (remove implicit and explicit Hs).\n\
\n\
  RETURNS: a new molecule with the Hs removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n";
      python::def("RemoveHs", removeHs,
                  (python::arg("mol"),python::arg("implicitOnly")=false),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());

      python::def("MergeQueryHs", (ROMol *(*)(const ROMol &))&MolOps::mergeQueryHs,
                  (python::arg("mol")),
                  "merges hydrogens into their neighboring atoms as queries",
                  python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Removes atoms matching a substructure query from a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as a substructure query\n\
\n\
    - onlyFrags: (optional) if this toggle is set, atoms will only be removed if\n\
      the entire fragment in which they are found is matched by the query.\n\
      See below for examples.\n\
      Default value is 0 (remove the atoms whether or not the entire fragment matches)\n\
\n\
  RETURNS: a new molecule with the substructure removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - DeleteSubstructs('CCOC','OC') -> 'CC'\n\
\n\
    - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'\n\
\n\
    - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'\n\
\n\
    - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'\n\
\n";
      python::def("DeleteSubstructs", deleteSubstructs,
                  (python::arg("mol"),python::arg("query"),
                   python::arg("onlyFrags")=false),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());
      docString="Do a Murcko decomposition and return the scaffold";
      python::def("MurckoDecompose", MurckoDecompose,
                  (python::arg("mol")),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());                  
      docString="Combine the atoms from two molecules to produce a third";
      python::def("CombineMols", combineMols,
                  (python::arg("mol1"),python::arg("mol2"),
                   python::arg("offset")=RDGeom::Point3D(0,0,0)),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());                  

      

      // ------------------------------------------------------------------------
      docString="Replaces atoms matching a substructure query in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as a substructure query\n\
\n\
    - replacement: the molecule to be used as the replacement\n\
\n\
    - replaceAll: (optional) if this toggle is set, all substructures matching\n\
      the query will be replaced in a single result, otherwise each result will\n\
      contain a separate replacement.\n\
      Default value is False (return multiple replacements)\n\
    - replacementConnectionPoint: (optional) index of the atom in the replacement that\n\
      the bond should be made to.\n\
\n\
  RETURNS: a tuple of new molecules with the substructures replaced removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSubstructs('CCOC','OC','NC') -> ('CCNC',)\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','NC') -> ('COCCNC','CNCCOC')\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','NC',True) -> ('CNCCNC',)\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','CN',True,1) -> ('CNCCNC',)\n\
\n";
      python::def("ReplaceSubstructs", replaceSubstructures,
                  (python::arg("mol"),python::arg("query"),
                   python::arg("replacement"),
                   python::arg("replaceAll")=false,
                   python::arg("replacementConnectionPoint")=0),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Adds named recursive queries to atoms\n";
      python::def("MolAddRecursiveQueries",addRecursiveQueriesHelper,
                  (python::arg("mol"),python::arg("queries"),
                   python::arg("propName")),
                  docString.c_str());

      docString="reads query definitions from a simply formatted file\n";
      python::def("ParseMolQueryDefFile",parseQueryDefFileHelper,
                  (python::arg("fileobj"),python::arg("standardize")=true,
                   python::arg("delimiter")="\t",python::arg("comment")="//",
                   python::arg("nameColumn")=0,python::arg("smartsColumn")=1),
                  docString.c_str());

     
      // ------------------------------------------------------------------------
      docString="Returns the molecule's topological distance matrix.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - useBO: (optional) toggles use of bond orders in calculating the distance matrix.\n\
      Default value is 0.\n\
\n\
    - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the\n\
      matrix (to return a \"Balaban\" distance matrix).\n\
      Default value is 0.\n\
\n\
    - force: (optional) forces the calculation to proceed, even if there is a cached value.\n\
      Default value is 0.\n\
\n\
    - prefix: (optional, internal use) sets the prefix used in the property cache\n\
      Default value is "".\n\
\n\
  RETURNS: a Numeric array of floats with the distance matrix\n\
\n";
      python::def("GetDistanceMatrix", getDistanceMatrix,
                  (python::arg("mol"),python::arg("useBO")=false,
                   python::arg("useAtomWts")=false,
                   python::arg("force")=false,
                   python::arg("prefix")=""),
                  docString.c_str());
      // ------------------------------------------------------------------------
      docString="Returns the molecule's 3D distance matrix.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - confId: (optional) chooses the conformer Id to use\n\
      Default value is -1.\n\
\n\
    - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the\n\
      matrix (to return a \"Balaban\" distance matrix).\n\
      Default value is 0.\n\
\n\
    - force: (optional) forces the calculation to proceed, even if there is a cached value.\n\
      Default value is 0.\n\
\n\
    - prefix: (optional, internal use) sets the prefix used in the property cache\n\
      Default value is "".\n\
\n\
  RETURNS: a Numeric array of floats with the distance matrix\n\
\n";
      python::def("Get3DDistanceMatrix", get3DDistanceMatrix,
                  (python::arg("mol"),python::arg("confId")=-1,
                   python::arg("useAtomWts")=false,
                   python::arg("force")=false,
                   python::arg("prefix")=""),
                  docString.c_str());
      // ------------------------------------------------------------------------
      docString="Returns the molecule's adjacency matrix.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - useBO: (optional) toggles use of bond orders in calculating the matrix.\n\
      Default value is 0.\n\
\n\
    - emptyVal: (optional) sets the elements of the matrix between non-adjacent atoms\n\
      Default value is 0.\n\
\n\
    - force: (optional) forces the calculation to proceed, even if there is a cached value.\n\
      Default value is 0.\n\
\n\
    - prefix: (optional, internal use) sets the prefix used in the property cache\n\
      Default value is "".\n\
\n\
  RETURNS: a Numeric array of floats containing the adjacency matrix\n\
\n";
      python::def("GetAdjacencyMatrix", getAdjacencyMatrix, 
                  (python::arg("mol"), python::arg("useBO")=false,
                   python::arg("emptyVal")=0,
                   python::arg("force")=false,
                   python::arg("prefix")=""),
                  docString.c_str());

      
      // ------------------------------------------------------------------------
      docString="Kekulizes the molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the \n\
      molecule will be marked non-aromatic following the kekulization.\n\
      Default value is 0.\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
      python::def("Kekulize", kekulizeMol,
                  (python::arg("mol"),python::arg("clearAromaticFlags")=false),
                  docString.c_str());
      

      // ------------------------------------------------------------------------
      docString="cleans up certain common bad functionalities in the molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
      python::def("Cleanup", cleanupMol,
                  (python::arg("mol")),
                  docString.c_str());
      
      
      // ------------------------------------------------------------------------
      docString="does aromaticity perception\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
      python::def("SetAromaticity", setAromaticityMol,
                  (python::arg("mol")),
                  docString.c_str());
      docString="finds conjugated bonds\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
      python::def("SetConjugation", setConjugationMol,
                  (python::arg("mol")),
                  docString.c_str());
      docString="Assigns hybridization states to atoms\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
      python::def("SetHybridization", setHybridizationMol,
                  (python::arg("mol")),
                  docString.c_str());
      docString="Assigns radical counts to atoms\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
      python::def("AssignRadicals", assignRadicalsMol,
                  (python::arg("mol")),
                  docString.c_str());
      


      // ------------------------------------------------------------------------
      docString="Finds all subgraphs of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target number of bonds for the subgraphs.\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified\n\
      atom will be returned.\n\
\n\
  RETURNS: a tuple of 2-tuples with bond IDs\n\
\n\
  NOTES: \n\
\n\
   - Difference between _subgraphs_ and _paths_ :: \n\
\n\
       Subgraphs are potentially branched, whereas paths (in our \n\
       terminology at least) cannot be.  So, the following graph: \n\
\n\
            C--0--C--1--C--3--C\n\
                  |\n\
                  2\n\
                  |\n\
                  C\n\
  has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)\n\
  but only 2 _paths_ of length 3: (0,1,3),(2,1,3)\n\
\n";
      python::def("FindAllSubgraphsOfLengthN", &findAllSubgraphsOfLengthN,
                  (python::arg("mol"),python::arg("length"),
                   python::arg("useHs")=false,
                   python::arg("rootedAtAtom")=-1),
                  docString.c_str());
      // ------------------------------------------------------------------------
      docString="Finds all subgraphs of a particular length in a molecule\n\
  See documentation for FindAllSubgraphsOfLengthN for definitions\n\
\n";
      python::def("FindAllSubgraphsOfLengthMToN", &findAllSubgraphsOfLengthsMtoNHelper,
                  (python::arg("mol"),python::arg("min"),python::arg("max"),
                   python::arg("useHs")=false,
                   python::arg("rootedAtAtom")=-1),
                  docString.c_str());
      // ------------------------------------------------------------------------
      docString="Finds unique subgraphs of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target number of bonds for the subgraphs.\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - useBO: (optional) Toggles use of bond orders in distinguishing one subgraph from\n\
      another.\n\
      Defaults to 1.\n\
\n\
    - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified\n\
      atom will be returned.\n\
\n\
  RETURNS: a tuple of tuples with bond IDs\n\
\n\
\n";
      python::def("FindUniqueSubgraphsOfLengthN", &findUniqueSubgraphsOfLengthN, 
                  (python::arg("mol"),python::arg("length"),
                   python::arg("useHs")=false,python::arg("useBO")=true,
                   python::arg("rootedAtAtom")=-1),
                  docString.c_str());
                  
      // ------------------------------------------------------------------------
      docString="Finds all paths of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target length for the paths.\n\
\n\
    - useBonds: (optional) toggles the use of bond indices in the paths.\n\
      Otherwise atom indices are used.  *Note* this behavior is different\n\
      from that for subgraphs.\n\
      Defaults to 1.\n\
\n\
    - rootedAtAtom: (optional) if nonzero, only paths from the specified\n\
      atom will be returned.\n\
\n\
  RETURNS: a tuple of tuples with IDs for the bonds.\n\
\n\
  NOTES: \n\
\n\
   - Difference between _subgraphs_ and _paths_ :: \n\
\n\
       Subgraphs are potentially branched, whereas paths (in our \n\
       terminology at least) cannot be.  So, the following graph: \n\
\n\
            C--0--C--1--C--3--C\n\
                  |\n\
                  2\n\
                  |\n\
                  C\n\
\n\
       has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)\n\
       but only 2 _paths_ of length 3: (0,1,3),(2,1,3)\n\
\n";
      python::def("FindAllPathsOfLengthN", &findAllPathsOfLengthN, 
                  (python::arg("mol"),python::arg("length"),
                   python::arg("useBonds")=true,python::arg("useHs")=false,
                   python::arg("rootedAtAtom")=-1),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Finds the bonds within a certain radius of an atom in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - radius: an integer with the target radius for the environment.\n\
\n\
    - rootedAtAtom: the atom to consider\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
  RETURNS: a vector of bond IDs\n\
\n\
\n";
      python::def("FindAtomEnvironmentOfRadiusN", &findAtomEnvironmentOfRadiusN,
                  (python::arg("mol"),python::arg("radius"),
                   python::arg("rootedAtAtom"),
                   python::arg("useHs")=false),
                  docString.c_str());

      python::def("PathToSubmol",pathToSubmolHelper,
                  (python::arg("mol"),python::arg("path"),
                   python::arg("useQuery")=false,
                   python::arg("atomMap")=python::object()),
                  "",
                  python::return_value_policy<python::manage_new_object>());
      
      // ------------------------------------------------------------------------
      docString="Finds the disconnected fragments from a molecule.\n\
\n\
  For example, for the molecule 'CC(=O)[O-].[NH3+]C' GetMolFrags() returns\n\
  ((0, 1, 2, 3), (4, 5))\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - asMols: (optional) if this is provided and true, the fragments\n\
      will be returned as molecules instead of atom ids.\n\
    - sanitizeFrags: (optional) if this is provided and true, the fragments\n\
      molecules will be sanitized before returning them.\n\
\n\
  RETURNS: a tuple of tuples with IDs for the atoms in each fragment\n\
           or a tuple of molecules.\n\
\n";
      python::def("GetMolFrags", &GetMolFrags,
                  (python::arg("mol"),python::arg("asMols")=false,
                   python::arg("sanitizeFrags")=true),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Splits a molecule into pieces based on PDB residue information.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - whiteList: only residues in this list will be returned\n\
    - negateList: if set, negates the white list inclusion logic\n\
\n\
  RETURNS: a dictionary keyed by residue name with molecules as the values\n\
\n";
      python::def("SplitMolByPDBResidues", &splitMolByPDBResidues,
                  (python::arg("mol"),
                   python::arg("whiteList")=python::object(),
                   python::arg("negateList")=false),
                  docString.c_str());
      // ------------------------------------------------------------------------
      docString="Splits a molecule into pieces based on PDB chain information.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - whiteList: only residues in this list will be returned\n\
    - negateList: if set, negates the white list inclusion logic\n\
\n\
  RETURNS: a dictionary keyed by chain id with molecules as the values\n\
\n";
      python::def("SplitMolByPDBChainId", &splitMolByPDBChainId,
                  (python::arg("mol"),
                   python::arg("whiteList")=python::object(),
                   python::arg("negateList")=false),
                  docString.c_str());



      // ------------------------------------------------------------------------
      docString="Returns the formal charge for the molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n";
      python::def("GetFormalCharge", &MolOps::getFormalCharge,docString.c_str());


      // ------------------------------------------------------------------------
      docString="Does the CIP stereochemistry assignment \n\
  for the molecule's atoms (R/S) and double bond (Z/E).\n\
  Chiral atoms will have a property '_CIPCode' indicating\n\
  their chiral code.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - cleanIt: (optional) if provided, atoms with a chiral specifier that aren't\n\
      actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,\n\
      etc.) will have their chiral code set to CHI_UNSPECIFIED\n\
    - force: (optional) causes the calculation to be repeated, even if it has already\n\
      been done\n\
    - flagPossibleStereoCenters (optional)   set the _ChiralityPossible property on\n\
      atoms that are possible stereocenters\n\
\n";
      python::def("AssignStereochemistry", MolOps::assignStereochemistry,
                  (python::arg("mol"),python::arg("cleanIt")=false,python::arg("force")=false,
                   python::arg("flagPossibleStereoCenters")=false),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Removes all stereochemistry info from the molecule.\n\
\n";
      python::def("RemoveStereochemistry", MolOps::removeStereochemistry,
                  (python::arg("mol")),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Sets the chiral tags on a molecule's atoms based on \n\
  a 3D conformation.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - confId: the conformer id to use, -1 for the default \n\
    - replaceExistingTags: if True, existing stereochemistry information will be cleared \n\
                           before running the calculation. \n\
\n";
      python::def("AssignAtomChiralTagsFromStructure", MolOps::assignChiralTypesFrom3D,
                  (python::arg("mol"),python::arg("confId")=-1,python::arg("replaceExistingTags")=true),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Returns an RDKit topological fingerprint for a molecule\n\
\n\
  Explanation of the algorithm below.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - minPath: (optional) minimum number of bonds to include in the subgraphs\n\
      Defaults to 1.\n\
\n\
    - maxPath: (optional) maximum number of bonds to include in the subgraphs\n\
      Defaults to 7.\n\
\n\
    - fpSize: (optional) number of bits in the fingerprint\n\
      Defaults to 2048.\n\
\n\
    - nBitsPerHash: (optional) number of bits to set per path\n\
      Defaults to 2.\n\
\n\
    - useHs: (optional) include paths involving Hs in the fingerprint if the molecule\n\
      has explicit Hs.\n\
      Defaults to True.\n\
\n\
    - tgtDensity: (optional) fold the fingerprint until this minimum density has\n\
      been reached\n\
      Defaults to 0.\n\
\n\
    - minSize: (optional) the minimum size the fingerprint will be folded to when\n\
      trying to reach tgtDensity\n\
      Defaults to 128.\n\
\n\
    - branchedPaths: (optional) if set both branched and unbranched paths will be\n\
      used in the fingerprint.\n\
      Defaults to True.\n\
\n\
    - useBondOrder: (optional) if set both bond orders will be used in the path hashes\n\
      Defaults to True.\n\
\n\
    - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes\n\
      Defaults to empty.\n\
\n\
    - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs \n\
      starting from these atoms will be used.\n\
      Defaults to empty.\n\
\n\
    - atomBits: (optional) an empty list. If provided, the result will contain a list \n\
      containing the bits each atom sets.\n\
      Defaults to empty.\n\
\n\
  RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits\n\
\n\
  ALGORITHM:\n\
\n\
   This algorithm functions by find all subgraphs between minPath and maxPath in\n \
   length.  For each subgraph:\n\
\n\
     1) A hash is calculated.\n\
\n\
     2) The hash is used to seed a random-number generator\n\
\n\
     3) _nBitsPerHash_ random numbers are generated and used to set the corresponding\n\
        bits in the fingerprint\n\
\n\
\n";
      python::def("RDKFingerprint", wrapRDKFingerprintMol,
                  (python::arg("mol"),python::arg("minPath")=1,
                   python::arg("maxPath")=7,python::arg("fpSize")=2048,
                   python::arg("nBitsPerHash")=2,python::arg("useHs")=true,
                   python::arg("tgtDensity")=0.0,python::arg("minSize")=128,
                   python::arg("branchedPaths")=true,
                   python::arg("useBondOrder")=true,
                   python::arg("atomInvariants")=0,
                   python::arg("fromAtoms")=0,
                   python::arg("atomBits")=python::object()
                   ),
                  docString.c_str(),python::return_value_policy<python::manage_new_object>());
      python::scope().attr("_RDKFingerprint_version")=RDKit::RDKFingerprintMolVersion;

      // ------------------------------------------------------------------------
      docString="Returns a layered fingerprint for a molecule\n\
\n\
  NOTE: This function is experimental. The API or results may change from\n\
    release to release.\n\
\n\
  Explanation of the algorithm below.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - layerFlags: (optional) which layers to include in the fingerprint\n\
      See below for definitions. Defaults to all.\n\
\n\
    - minPath: (optional) minimum number of bonds to include in the subgraphs\n\
      Defaults to 1.\n\
\n\
    - maxPath: (optional) maximum number of bonds to include in the subgraphs\n\
      Defaults to 7.\n\
\n\
    - fpSize: (optional) number of bits in the fingerprint\n\
      Defaults to 2048.\n\
\n\
    - atomCounts: (optional) \n\
        if provided, this should be a list at least as long as the number of atoms\n\
        in the molecule. It will be used to provide the count of the number \n                      \
        of paths that set bits each atom is involved in.\n\
        NOTE: the list is not zeroed out here.\n\
\n\
    - setOnlyBits: (optional) \n\
        if provided, only bits that are set in this bit vector will be set\n\
        in the result. This is essentially the same as doing:\n\
           res &= setOnlyBits\n\
        but also has an impact on the atomCounts (if being used)\n\
\n\
    - branchedPaths: (optional) if set both branched and unbranched paths will be\n\
      used in the fingerprint.\n\
      Defaults to True.\n\
\n\
    - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs \n\
      starting from these atoms will be used.\n\
      Defaults to empty.\n\
\n\
  RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits\n\
\n\
  Layer definitions:\n\
     - 0x01: pure topology\n\
     - 0x02: bond order\n\
     - 0x04: atom types\n\
     - 0x08: presence of rings\n\
     - 0x10: ring sizes\n\
     - 0x20: aromaticity\n\
\n\
\n";
      python::def("LayeredFingerprint", wrapLayeredFingerprint,
                  (python::arg("mol"),
                   python::arg("layerFlags")=0xFFFFFFFF,
                   python::arg("minPath")=1,
                   python::arg("maxPath")=7,python::arg("fpSize")=2048,
                   python::arg("atomCounts")=python::list(),
                   python::arg("setOnlyBits")=(ExplicitBitVect *)0,
                   python::arg("branchedPaths")=true,
                   python::arg("fromAtoms")=0
                   ),
                  docString.c_str(),python::return_value_policy<python::manage_new_object>());
      python::scope().attr("_LayeredFingerprint_version")=RDKit::LayeredFingerprintMolVersion;
      python::scope().attr("LayeredFingerprint_substructLayers")=RDKit::substructLayers;

      // ------------------------------------------------------------------------
      docString="A fingerprint using SMARTS patterns \n\
\n\
  NOTE: This function is experimental. The API or results may change from\n\
    release to release.\n";
      python::def("PatternFingerprint", wrapPatternFingerprint,
                  (python::arg("mol"),
                   python::arg("fpSize")=2048,
                   python::arg("atomCounts")=python::list(),
                   python::arg("setOnlyBits")=(ExplicitBitVect *)0),
                  docString.c_str(),python::return_value_policy<python::manage_new_object>());

      docString="Set the wedging on single bonds in a molecule.\n \
   The wedging scheme used is that from Mol files.\n \
\n\
  ARGUMENTS:\n\
\n\
    - molecule: the molecule to update\n \
\n\
\n";
      python::def("WedgeMolBonds", WedgeMolBonds,
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Replaces sidechains in a molecule with dummy atoms for their attachment points.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
  RETURNS: a new molecule with the sidechains removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSidechains('CCC1CCC1','C1CCC1') -> '[Xa]C1CCC1'\n\
\n\
    - ReplaceSidechains('CCC1CC1','C1CCC1') -> ''\n\
\n\
    - ReplaceSidechains('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'\n\
\n";
      python::def("ReplaceSidechains", replaceSidechains,
                  (python::arg("mol"),python::arg("coreQuery")),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Removes the core of a molecule and labels the sidechains with dummy atoms.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
    - replaceDummies: toggles replacement of atoms that match dummies in the query\n\
\n\
    - labelByIndex: toggles labeling the attachment point dummy atoms with \n\
                    the index of the core atom they're attached to.\n\
\n\
    - requireDummyMatch: if the molecule has side chains that attach at points not\n\
                         flagged with a dummy, it will be rejected (None is returned)\n\
\n\
  RETURNS: a new molecule with the core removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceCore('CCC1CCC1','C1CCC1') -> 'CC[1*]'\n\
\n\
    - ReplaceCore('CCC1CC1','C1CCC1') -> ''\n\
\n\
    - ReplaceCore('C1CC2C1CCC2','C1CCC1') -> '[1*]C1CCC1[2*]'\n\
\n\
    - ReplaceCore('C1CNCC1','N') -> '[1*]CCCC[2*]'\n\
\n\
    - ReplaceCore('C1CCC1CN','C1CCC1[*]',False) -> '[1*]CN'\n\
\n";
      python::def("ReplaceCore", replaceCore,
                  (python::arg("mol"),python::arg("coreQuery"),
                   python::arg("replaceDummies")=true,
                   python::arg("labelByIndex")=false,
                   python::arg("requireDummyMatch")=false
                   ),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());


      docString="Return a new molecule with all BRICS bonds broken";
      python::def("FragmentOnBRICSBonds", MolFragmenter::fragmentOnBRICSBonds,
                  (python::arg("mol")),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());

      docString="Return a new molecule with all specified bonds broken\n\
\n\
  ARGUMENTS:\n\
\n\
      - mol            - the molecule to be modified\n\
      - bondIndices    - indices of the bonds to be broken\n\
      - addDummies  - toggles addition of dummy atoms to indicate where \n\
          bonds were broken\n\
      - dummyLabels - used to provide the labels to be used for the dummies.\n\
          the first element in each pair is the label for the dummy\n\
          that replaces the bond's beginAtom, the second is for the \n\
          dummy that replaces the bond's endAtom. If not provided, the\n\
          dummies are labeled with atom indices.\n\
      - bondTypes - used to provide the bond type to use between the\n\
          fragments and the dummy atoms. If not provided, defaults to single. \n\
\n\
  RETURNS:\n\
      a new Mol with the modifications\n\
";
      python::def("FragmentOnBonds", fragmentOnBondsHelper,
                  (python::arg("mol"),
                   python::arg("bondIndices"),
                   python::arg("addDummies")=true,
                   python::arg("dummyLabels")=python::object(),
                   python::arg("bondTypes")=python::object()),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());
      docString="fragment on some bonds";
      python::def("FragmentOnSomeBonds", fragmentOnSomeBondsHelper,
                  (python::arg("mol"),
                   python::arg("bondIndices"),
                   python::arg("numToBreak")=1,
                   python::arg("addDummies")=true,
                   python::arg("dummyLabels")=python::object(),
                   python::arg("bondTypes")=python::object()),
                  docString.c_str());


      // ------------------------------------------------------------------------
      docString="Adds a recursive query to an atom\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as the recursive query (this will be copied)\n\
\n\
    - atomIdx: the atom to modify\n\
\n\
    - preserveExistingQuery: (optional) if this is set, existing query information on the atom will be preserved\n\
\n\
  RETURNS: None\n\
\n";
      python::def("AddRecursiveQuery", addRecursiveQuery,
                  (python::arg("mol"),python::arg("query"),
                   python::arg("atomIdx"),python::arg("preserveExistingQuery")=true),
                  docString.c_str());

      // ------------------------------------------------------------------------
      docString="Returns a copy of a molecule with renumbered atoms\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - newOrder: the new ordering the atoms (should be numAtoms long)\n\
         for example: if newOrder is [3,2,0,1], then atom 3 in the original \n\
         molecule will be atom 0 in the new one\n\
\n\
\n";
      python::def("RenumberAtoms", renumberAtomsHelper,
                  (python::arg("mol"),python::arg("newOrder")),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());


    };
  };
}

void wrap_molops() {
  RDKit::molops_wrapper::wrap();
}

