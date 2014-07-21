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
%include "std_pair.i"
%include "std_vector.i"
%{
#include <RDGeneral/types.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Bond.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolPickler.h>
#include <DistGeom/BoundsMatrix.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolAlign/O3AAlignMolecules.h>
#include <GraphMol/MolDrawing/MolDrawing.h>
#include <GraphMol/MolDrawing/DrawingToSVG.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>

%}

%template(ROMol_Vect) std::vector< boost::shared_ptr<RDKit::ROMol> >;
%template(ROMol_Vect_Vect) std::vector< std::vector< boost::shared_ptr<RDKit::ROMol> > >;
%template(Atom_Vect) std::vector<RDKit::Atom*>;

// These prevent duplicate definitions in Java code 
%ignore RDKit::ROMol::getAtomDegree(const Atom *) const;
%ignore RDKit::ROMol::setAtomBookmark(Atom *,int);
%ignore RDKit::ROMol::clearAtomBookmark(const int, const Atom *);
%ignore RDKit::ROMol::setBondBookmark(Bond *,int);
%ignore RDKit::ROMol::clearBondBookmark(int, const Bond *);
%ignore RDKit::ROMol::replaceAtomBookmark(Atom *,int);
%ignore RDKit::ROMol::hasProp(std::string const) const ;
%ignore RDKit::ROMol::clearProp(std::string const) const ;
%ignore RDKit::ROMol::getAtomWithIdx(unsigned int) const ;
%ignore RDKit::ROMol::getBondWithIdx(unsigned int) const ;
%ignore RDKit::ROMol::getBondBetweenAtoms(unsigned int,unsigned int) const ;
%ignore RDKit::ROMol::getAtomNeighbors(Atom const *at) const;
%ignore RDKit::ROMol::getAtomNeighbors(ATOM_SPTR at) const;
%ignore RDKit::ROMol::getAtomBonds(Atom const *at) const;
%ignore RDKit::ROMol::getVertices() ;
%ignore RDKit::ROMol::getVertices() const ;
%ignore RDKit::ROMol::getEdges() ;
%ignore RDKit::ROMol::getEdges() const ;
%ignore RDKit::ROMol::getTopology() const ;


/*
 * Special handling for Conformer objects which should not be GCed until the molecule is destroyed
 * We want to modify the behavior of the Conformer coming into the addConformer method without 
 * impacting Conformer objects that are arguments to other methods. Therefore we define a pattern
 * that will trigger special handling of the Conformer input (the addConf method  match this pattern).
 * Then add the necessary Java code to modify the Conformer object to no longer be the owner of the
 * underlying C++ object.
 */
%ignore addConformer(Conformer * conf, bool assignId=false);
%rename(addConformer) RDKit::ROMol::addConf;
%typemap(javain) RDKit::Conformer * ownedConf "getCPtrAndReleaseControl($javainput)"
%typemap(javacode) RDKit::ROMol %{
  // Ensure that the GC doesn't collect this item,
  // as the underlying C++ class stores a shallow copy
  private long getCPtrAndReleaseControl(Conformer conf) {
    conf.setSwigCMemOwn(false);
    return Conformer.getCPtr(conf);
  }
%}
%include <GraphMol/ROMol.h>

/* For the time being, assume all properties will be strings */
%template(setProp)  RDKit::ROMol::setProp<std::string>;


%newobject removeHs;
%newobject addHs;
%newobject mergeQueryHs;
%newobject replaceCore;
%newobject replaceSidechains;
%newobject deleteSubstructs;
%newobject getAtoms;


%extend RDKit::ROMol {
  std::string getProp(const std::string key){
    std::string res;
    self->getProp(key, res);
    return res;
  }

  /* Used in the addConformer modifications described above */
  unsigned int RDKit::ROMol::addConf(RDKit::Conformer * ownedConf, bool assignId=false) {
    return self->addConformer(ownedConf, assignId);
  }

  std::string MolToSmiles(bool doIsomericSmiles=false,bool doKekule=false, int rootedAtAtom=-1){
    return RDKit::MolToSmiles(*($self),doIsomericSmiles,doKekule,rootedAtAtom);
  }
  std::string MolToMolBlock(bool includeStereo=true, int confId=-1) {
    return RDKit::MolToMolBlock(*($self),includeStereo,confId);
  }
  void MolToMolFile(std::string fName,bool includeStereo=true, int confId=-1,bool kekulize=true) {
    RDKit::MolToMolFile(*($self), fName, includeStereo, confId, kekulize);
  }
  std::string MolToTPLText(std::string partialChargeProp="_GasteigerCharge", bool writeFirstConfTwice=false) {
    return RDKit::MolToTPLText(*($self), partialChargeProp, writeFirstConfTwice);
  }
  void MolToTPLFile(std::string fName,
    std::string partialChargeProp="_GasteigerCharge",
    bool writeFirstConfTwice=false) {
    RDKit::MolToTPLFile(*($self), fName, partialChargeProp, writeFirstConfTwice);
  }

  std::string MolToPDBBlock(int confId=-1,unsigned int flavor=0) {
    return RDKit::MolToPDBBlock(*($self),confId,flavor);
  }
  void MolToPDBFile(std::string fName,int confId=-1,unsigned int flavor=0) {
    RDKit::MolToPDBFile(*($self), fName, confId, flavor);
  }

  bool hasSubstructMatch(RDKit::ROMol &query,bool useChirality=false){
    RDKit::MatchVectType mv;
    return SubstructMatch(*($self),query,mv,true,useChirality);
  };

  /* From MolOps, Substruct/SubstructMatch */
  std::vector<std::pair<int, int> > getSubstructMatch(RDKit::ROMol &query,bool useChirality=false){
    RDKit::MatchVectType mv;
    SubstructMatch(*($self),query,mv,true,useChirality);
    return mv;
  };

  std::vector< std::vector<std::pair<int, int> > > getSubstructMatches(RDKit::ROMol &query,bool uniquify=true, bool useChirality=false){
    std::vector<RDKit::MatchVectType> mv;
    SubstructMatch(*($self),query,mv,uniquify,true,useChirality);
    return mv;
  };

  RDKit::ROMol *deleteSubstructs(const RDKit::ROMol &query, bool replaceAll) {
    return RDKit::deleteSubstructs(*($self), query, replaceAll);
  };

  std::vector<RDKit::ROMOL_SPTR> replaceSubstructs(const RDKit::ROMol &query,
                                        const RDKit::ROMol &replacement,
                                        bool replaceAll=false) {
    return RDKit::replaceSubstructs(*($self), query, replacement, replaceAll);
  };

  RDKit::ROMol *replaceSidechains(const RDKit::ROMol &coreQuery) {
    return RDKit::replaceSidechains(*($self), coreQuery);
  };

  RDKit::ROMol *replaceCore(const RDKit::ROMol &coreQuery, bool replaceDummies=true,bool labelByIndex=false) {
    return RDKit::replaceCore(*($self), coreQuery, replaceDummies, labelByIndex);
  };

  void AssignAtomCIPRanks(RDKit::INT_VECT &ranks) {
    RDKit::Chirality::assignAtomCIPRanks(*($self), ranks);
  };

  /* Methods from MolFileStereoChem.h */
  void DetectBondStereoChemistry(const RDKit::Conformer *conf) {
    RDKit::DetectBondStereoChemistry(*($self), conf);
  };

  void WedgeMolBonds(const RDKit::Conformer *conf) {
    RDKit::WedgeMolBonds(*($self), conf);
  };

  void pickBondsToWedge() {
    RDKit::pickBondsToWedge(*($self));
  };

  void ClearSingleBondDirFlags() {
    RDKit::ClearSingleBondDirFlags(*($self));
  };

  /* Methods from ConjugHybrid.cpp */
  void setConjugation() {
    RDKit::MolOps::setConjugation(*($self));
  } 

  void setHybridization() {
    RDKit::MolOps::setHybridization(*($self));
  }

  /* From GraphMol/Depictor/RDDepictor.h */
  unsigned int compute2DCoords(const RDGeom::INT_POINT2D_MAP *coordMap=0,
                               bool canonOrient=false,
                               bool clearConfs=true,
                               unsigned int nFlipsPerSample=0,
                               unsigned int nSamples=0,
                               int sampleSeed=0,
                               bool permuteDeg4Nodes=false) {
    return RDDepict::compute2DCoords(*($self),
                               coordMap,
                               canonOrient, 
                               clearConfs, 
                               nFlipsPerSample, 
                               nSamples,
                               sampleSeed,
                               permuteDeg4Nodes);
  }

  unsigned int compute2DCoords(RDKit::ROMol &templ){
    RDKit::MatchVectType matchVect;
    if(templ.getNumConformers() && SubstructMatch(*($self),templ,matchVect)){
      RDGeom::INT_POINT2D_MAP coordMap;
      RDKit::Conformer conf=templ.getConformer();
      for(RDKit::MatchVectType::const_iterator iter=matchVect.begin();
          iter!=matchVect.end();++iter){
        RDGeom::Point2D pt;
        pt.x = conf.getAtomPos(iter->first).x;
        pt.y = conf.getAtomPos(iter->first).y;
        coordMap[iter->second]=pt;
      }
      return RDDepict::compute2DCoords(*($self),&coordMap);
    } else {
      return RDDepict::compute2DCoords(*($self),0);
    }
  }

  
  unsigned int compute2DCoordsMimicDistMat(const RDDepict::DOUBLE_SMART_PTR *dmat=0,
                                           bool canonOrient=true,
                                           bool clearConfs=true,
                                           double weightDistMat=0.5,
                                           unsigned int nFlipsPerSample=3,
                                           unsigned int nSamples=100,
                                           int sampleSeed=25,
                                           bool permuteDeg4Nodes=true) {
    return RDDepict::compute2DCoordsMimicDistMat(*($self),
                                                 dmat,
                                                 canonOrient,
                                                 clearConfs,
                                                 weightDistMat,
                                                 nFlipsPerSample,
                                                 nSamples,
                                                 sampleSeed,
                                                 permuteDeg4Nodes);

  }

  /* From FindRings.cpp, MolOps.h */
  int findSSSR(RDKit::VECT_INT_VECT &res) {
    return RDKit::MolOps::findSSSR(*($self), res);
  };
  int findSSSR() {
    return RDKit::MolOps::findSSSR(*($self));
  };
  int symmetrizeSSSR(RDKit::VECT_INT_VECT &res) {
    return RDKit::MolOps::symmetrizeSSSR(*($self), res);
  }
  int symmetrizeSSSR() {
    return RDKit::MolOps::symmetrizeSSSR(*($self));
  }

  /* From Matrices.cpp, MolOps.h */
  double *getDistanceMat(bool useBO=false,
                         bool useAtomWts=false,
                         bool force=false,
                         const char *propNamePrefix=0) {
    return RDKit::MolOps::getDistanceMat(*($self), useBO, useAtomWts, force, propNamePrefix);
  }

  double *getDistanceMat(const std::vector<int> &activeAtoms,
                         const std::vector<const Bond *> &bonds,
                         bool useBO=false,
                         bool useAtomWts=false) {
    return RDKit::MolOps::getDistanceMat(*($self), activeAtoms, bonds, useBO, useAtomWts);
  }

  double *getAdjacencyMatrix(bool useBO=false,
                             int emptyVal=0,
                             bool force=false,
                             const char *propNamePrefix=0) {
    return RDKit::MolOps::getAdjacencyMatrix(*($self), useBO, emptyVal, force, propNamePrefix);
  }

  RDKit::INT_LIST getShortestPath(int aid1, int aid2) {
    return RDKit::MolOps::getShortestPath(*($self), aid1, aid2);
  }
  /* From MolTransforms.h */
  void transformMolsAtoms(RDGeom::Transform3D &tform) {
    MolTransforms::transformMolsAtoms(($self), tform);
  }
  void canonicalizeMol(bool normalizeCovar=false, bool ignoreHs=true) {
    MolTransforms::canonicalizeMol(*($self), normalizeCovar, ignoreHs);
  }

  /* From Python wrappers -- implied functionality */
  std::vector<RDKit::Atom*> *getAtoms() {
    int c = ($self)->getNumAtoms();
    std::vector<RDKit::Atom*> *atoms = new std::vector<RDKit::Atom*>;
    for (int i = 0; i < c; i++) {
      RDKit::Atom* a = ($self)->getAtomWithIdx(i);
      atoms->push_back(a);
    }
    return atoms;
  }

  /* From MolPickler.h */
  std::vector<int> ToBinary(){
    std::string sres;
    RDKit::MolPickler::pickleMol(*($self),sres);
    std::vector<int> res(sres.length());
    std::copy(sres.begin(),sres.end(),res.begin());
    return res;
  };
  static RDKit::ROMOL_SPTR MolFromBinary(std::vector<int> pkl){
    std::string sres;
    sres.resize(pkl.size());
    std::copy(pkl.begin(),pkl.end(),sres.begin());
    RDKit::ROMol *res=new RDKit::ROMol(sres);
    return RDKit::ROMOL_SPTR(res);
  }

  /* From AddHs.cpp */
  RDKit::ROMol *addHs(bool explicitOnly,bool addCoords=false){
    return RDKit::MolOps::addHs(*($self), explicitOnly, addCoords);
  }

  RDKit::ROMol *removeHs(bool implicitOnly,bool updateExplicitCount=false,bool sanitize=false) {
    return RDKit::MolOps::removeHs(*($self), implicitOnly, updateExplicitCount, sanitize);
  }

  RDKit::ROMol *mergeQueryHs() {
    return RDKit::MolOps::mergeQueryHs(*($self));
  }

  /* From GraphMol/MolAlign/AlignMolecules */
  double alignMol(const RDKit::ROMol &refMol, 
                  int prbCid=-1, int refCid=-1,
                  const std::vector<std::pair<int,int> > *atomMap=0, 
                  const RDNumeric::DoubleVector *weights=0, 
                  bool reflect=false, unsigned int maxIters=50) {
    return RDKit::MolAlign::alignMol(*($self), refMol, prbCid, refCid, atomMap, weights, reflect, maxIters);
  }

  void alignMolConformers(ROMol &mol, const std::vector<unsigned int> *atomIds=0,
                          const std::vector<unsigned int> *confIds=0,
                          const RDNumeric::DoubleVector  *weights=0, 
                          bool reflect=false, unsigned int maxIters=50) {
    RDKit::MolAlign::alignMolConformers(*($self), atomIds, confIds, weights, reflect, maxIters);
  }

  /* From GraphMol/MolAlign/AlignMolecules */
  std::pair<double,double> O3AAlignMol(RDKit::ROMol &refMol, 
                                       int prbCid=-1, int refCid=-1,
                                       bool reflect=false, unsigned int maxIters=50,
                                       unsigned int accuracy=0) {
    RDKit::MMFF::MMFFMolProperties prbMP(*($self));
    RDKit::MMFF::MMFFMolProperties refMP(refMol);
    
    RDKit::MolAlign::O3A o3a(*($self), refMol, &prbMP, &refMP, RDKit::MolAlign::O3A::MMFF94,
                             prbCid, refCid,
                             reflect,maxIters,accuracy);
    double rmsd=o3a.align();
    double score = o3a.score();
    return std::make_pair(rmsd,score);
  }

  void computeGasteigerCharges(const RDKit::ROMol *mol,int nIter=12,bool throwOnParamFailure=false){
    RDKit::computeGasteigerCharges(*mol,nIter,throwOnParamFailure);
  }
  void computeGasteigerCharges(const RDKit::ROMol *mol,
                               std::vector<double> &charges,
                               int nIter=12,bool throwOnParamFailure=false){
    RDKit::computeGasteigerCharges(*mol,charges,nIter,throwOnParamFailure);
  }
}




%extend RDKit::ROMol {
  std::string ToSVG(int lineWidthMult=2,int fontSize=50){
    if(lineWidthMult<0) lineWidthMult *=2;
    if(fontSize<0) fontSize*=2;
    std::vector<int> drawing=RDKit::Drawing::MolToDrawing(*($self),0);
    std::string svg=RDKit::Drawing::DrawingToSVG(drawing,lineWidthMult,fontSize);
    return svg;
  }
  std::string ToSVG(const std::vector<int> &highlightAtoms,
                    int lineWidthMult=2,int fontSize=50){
    if(lineWidthMult<0) lineWidthMult *=2;
    if(fontSize<0) fontSize*=2;
    std::vector<int> drawing=RDKit::Drawing::MolToDrawing(*($self),&highlightAtoms);
    std::string svg=RDKit::Drawing::DrawingToSVG(drawing,lineWidthMult,fontSize);
    return svg;
  }
}
