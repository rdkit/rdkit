/*
*
*  Copyright (c) 2010-2021, Novartis Institutes for BioMedical Research Inc.
*   and other RDKit contributors
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
*       products derived from this software without specific prior written
        permission.
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
%include "std_string.i"
%include "std_vector.i"

#ifdef SWIGCSHARP
%include <std_unique_ptr.i>
%unique_ptr(RDKit::RWMol)
#endif

%{
#include <RDGeneral/types.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <GraphMol/Bond.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/MolPickler.h>
#include <DistGeom/BoundsMatrix.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolAlign/O3AAlignMolecules.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/Chirality.h>
#include <sstream>
%}

%template(ROMol_Vect) std::vector< boost::shared_ptr<RDKit::ROMol> >;
%template(ROMol_Vect_Vect) std::vector< std::vector< boost::shared_ptr<RDKit::ROMol> > >;
%template(Atom_Vect) std::vector<RDKit::Atom*>;
%template(StereoGroup_Vect) std::vector<RDKit::StereoGroup>;
%template(UChar_Vect) std::vector<unsigned char>;

// These prevent duplicate definitions in Java code
%ignore RDKit::ROMol::hasProp(std::string const) const ;
%ignore RDKit::ROMol::clearProp(std::string const) const ;
%ignore RDKit::ROMol::getAtomWithIdx(unsigned int) const ;
%ignore RDKit::ROMol::getBondWithIdx(unsigned int) const ;
%ignore RDKit::ROMol::getBondBetweenAtoms(unsigned int,unsigned int) const ;
%ignore RDKit::ROMol::getAtomNeighbors(Atom const *at) const;
%ignore RDKit::ROMol::getAtomBonds(Atom const *at) const;
%ignore RDKit::ROMol::atomNeighbors(Atom const *at) const;
%ignore RDKit::ROMol::atomBonds(Atom const *at) const;
%ignore RDKit::ROMol::atomNeighbors(Atom const *at);
%ignore RDKit::ROMol::atomBonds(Atom const *at);
%ignore RDKit::ROMol::atoms();
%ignore RDKit::ROMol::atoms() const;
%ignore RDKit::ROMol::bonds();
%ignore RDKit::ROMol::bonds() const;

%ignore RDKit::ROMol::getVertices() ;
%ignore RDKit::ROMol::getVertices() const ;
%ignore RDKit::ROMol::getEdges() ;
%ignore RDKit::ROMol::getEdges() const ;
%ignore RDKit::ROMol::getTopology() const ;


#ifdef SWIGJAVA
%typemap(jni) std::string RDKit::ROMol::toByteArray "jbyteArray"
%typemap(jtype) std::string RDKit::ROMol::toByteArray "byte[]"
%typemap(jstype) std::string RDKit::ROMol::toByteArray "byte[]"
%typemap(javaout) std::string RDKit::ROMol::toByteArray {
  return $jnicall;
}
%typemap(out) std::string RDKit::ROMol::toByteArray {
  $result = JCALL1(NewByteArray, jenv, $1.size());
  JCALL4(SetByteArrayRegion, jenv, $result, 0, $1.size(), (const jbyte*)$1.c_str());
}
#endif

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
  public static ROMol fromByteArray(byte[] pkl) {
    UChar_Vect vec = null;
    try {
      vec = new UChar_Vect();
      vec.reserve(pkl.length);
      for (int i = 0; i < pkl.length; ++i) {
        vec.add((byte)pkl[i]);
      }
      return ROMol.fromUCharVect(vec);
    } finally {
      if (vec != null) {
        vec.delete();
      }
    }
  }
%}
%typemap(cscode) RDKit::ROMol %{
  public static ROMol FromByteArray(byte[] pkl) {
    UChar_Vect vec = null;
    try {
      vec = new UChar_Vect();
      vec.Capacity = pkl.Length;
      for (int i = 0; i < pkl.Length; ++i) {
        vec.Add((byte)pkl[i]);
      }
      return ROMol.fromUCharVect(vec);
    } finally {
      if (vec != null) {
        vec.Dispose();
      }
    }
  }
  public
   byte[] ToByteArray(int propertyFlags = -1) {
     UChar_Vect vec = null;
     try {
       vec = toUCharVect(propertyFlags);
       byte[] res = new byte[vec.Count];
       vec.CopyTo(res);
       return res;
     } finally {
       if (vec != null) {
         vec.Dispose();
       }
     }
   }
%}
%include <GraphMol/ROMol.h>

%ignore SubstructMatch;
%include <GraphMol/Substruct/SubstructMatch.h>

%ignore RDKit::MolPickler;
#ifdef SWIGJAVA
%include "enums.swg"
%javaconst(1);
#endif
%include <GraphMol/MolPickler.h>
#ifdef SWIGJAVA
%javaconst(0);
#endif



%newobject removeHs;
%newobject addHs;
%newobject mergeQueryHs;
%newobject replaceCore;
%newobject replaceSidechains;
%newobject deleteSubstructs;
%newobject getAtoms;
%newobject getBonds;
%newobject getAtomNeighbors;
%newobject getAtomBonds;

%{
#ifdef RDK_BUILD_COORDGEN_SUPPORT
bool getPreferCoordGen() {
  return RDDepict::preferCoordGen;
}
void setPreferCoordGen(bool val) {
  RDDepict::preferCoordGen = val;
}
#else
bool getPreferCoordGen() {
  return false;
}
void setPreferCoordGen(bool val) {
}
#endif
%}

bool getPreferCoordGen();
void setPreferCoordGen(bool);

%{
bool getUseLegacyStereoPerception() {
  return RDKit::Chirality::getUseLegacyStereoPerception();
}
void setUseLegacyStereoPerception(bool val) {
  RDKit::Chirality::setUseLegacyStereoPerception(val);
}
bool getAllowNontetrahedralChirality() {
  return RDKit::Chirality::getAllowNontetrahedralChirality();
}
void setAllowNontetrahedralChirality(bool val) {
  RDKit::Chirality::setAllowNontetrahedralChirality(val);
}
%}

bool getUseLegacyStereoPerception();
void setUseLegacyStereoPerception(bool);
bool getAllowNontetrahedralChirality();
void setAllowNontetrahedralChirality(bool);

#ifdef SWIGJAVA
%javamethodmodifiers RDKit::ROMol::fromUCharVect "private";
#endif
#ifdef SWIGCSHARP
%csmethodmodifiers RDKit::ROMol::fromUCharVect "private";
%csmethodmodifiers RDKit::ROMol::toUCharVect "private";
#endif

%{
  /* From MolPickler.h */
void setDefaultPickleProperties(unsigned int propertyFlags) {
  RDKit::MolPickler::setDefaultPickleProperties(propertyFlags);
}
unsigned int getDefaultPickleProperties() {
  return RDKit::MolPickler::getDefaultPickleProperties();
}
%}

void setDefaultPickleProperties(unsigned int propertyFlags);
unsigned int getDefaultPickleProperties();

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

  std::string MolToSmiles(bool doIsomericSmiles=true, bool doKekule=false, int rootedAtAtom=-1, bool canonical=true,
                          bool allBondsExplicit=false, bool allHsExplicit=false, bool doRandom=false) {
    return RDKit::MolToSmiles(*($self), doIsomericSmiles, doKekule, rootedAtAtom, canonical, allBondsExplicit, allHsExplicit, doRandom);
  }
  std::string MolToSmiles(const RDKit::SmilesWriteParams &params) {
    return RDKit::MolToSmiles(*($self), params);
  }
  std::string MolToMolBlock(bool includeStereo=true, int confId=-1, bool kekulize=true, bool forceV3000=false) {
    return RDKit::MolToMolBlock(*($self), includeStereo, confId, kekulize, forceV3000);
  }
  void MolToMolFile(std::string fName,bool includeStereo=true, int confId=-1, bool kekulize=true, bool forceV3000=false) {
    RDKit::MolToMolFile(*($self), fName, includeStereo, confId, kekulize, forceV3000);
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
  std::string MolToSequence() {
    return RDKit::MolToSequence(*($self));
  }
  std::string MolToFASTA() {
    return RDKit::MolToFASTA(*($self));
  }
  std::string MolToHELM() {
    return RDKit::MolToHELM(*($self));
  }

  std::string MolToCMLBlock(int confId=-1, bool kekulize=true) {
    return RDKit::MolToCMLBlock(*($self), confId, kekulize);
  }
  void MolToCMLFile(std::string fName, int confId=-1, bool kekulize=true) {
    RDKit::MolToCMLFile(*($self), fName, confId, kekulize);
  }

  std::string MolToXYZBlock(int confId=-1) {
    return RDKit::MolToXYZBlock(*($self), confId);
  }
  void MolToXYZFile(std::string fName, int confId=-1) {
    RDKit::MolToXYZFile(*($self), fName, confId);
  }

  bool hasSubstructMatch(RDKit::ROMol &query,bool useChirality=false){
    RDKit::MatchVectType mv;
    return SubstructMatch(*($self),query,mv,true,useChirality);
  };

  bool hasSubstructMatch(RDKit::ROMol &query,RDKit::SubstructMatchParameters ps){
    ps.maxMatches = 1;
    std::vector<RDKit::MatchVectType> mv = SubstructMatch(*($self),query,ps);
    return mv.size()>0;
  };

  std::vector<std::pair<int, int> > getSubstructMatch(RDKit::ROMol &query,RDKit::SubstructMatchParameters ps){
    std::vector<RDKit::MatchVectType> mvs = SubstructMatch(*($self),query,ps);
    RDKit::MatchVectType mv;
    if(mvs.size()) mv = mvs[0];
    return mv;
  };

  std::vector< std::vector<std::pair<int, int> > > getSubstructMatches(RDKit::ROMol &query,RDKit::SubstructMatchParameters ps){
    std::vector<RDKit::MatchVectType> mvs = SubstructMatch(*($self),query,ps);
    return mvs;
  };

  bool hasSubstructMatch(RDKit::MolBundle & query,
                         RDKit::SubstructMatchParameters ps) {
    ps.maxMatches = 1;
    std::vector<RDKit::MatchVectType> mv = SubstructMatch(*($self), query, ps);
    return mv.size() > 0;
  };

  std::vector<std::pair<int, int>> getSubstructMatch(
      RDKit::MolBundle & query, RDKit::SubstructMatchParameters ps) {
    std::vector<RDKit::MatchVectType> mvs = SubstructMatch(*($self), query, ps);
    RDKit::MatchVectType mv;
    if (mvs.size()) mv = mvs[0];
    return mv;
  };

  std::vector<std::vector<std::pair<int, int>>> getSubstructMatches(
      RDKit::MolBundle & query, RDKit::SubstructMatchParameters ps) {
    std::vector<RDKit::MatchVectType> mvs = SubstructMatch(*($self), query, ps);
    return mvs;
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

  /* Methods from SubstructUtils.h */
  const RDKit::MatchVectType &getMostSubstitutedCoreMatch(const RDKit::ROMol& core,
    const std::vector<RDKit::MatchVectType>& matches) {
    return RDKit::getMostSubstitutedCoreMatch(*($self), core, matches);
  };

  std::vector<RDKit::MatchVectType> sortMatchesByDegreeOfCoreSubstitution(
    const RDKit::ROMol& core, const std::vector<RDKit::MatchVectType>& matches) {
    return RDKit::sortMatchesByDegreeOfCoreSubstitution(*($self), core, matches);
  };

  /* Methods from MolFileStereoChem.h */
  void DetectBondStereoChemistry(const RDKit::Conformer *conf) {
    RDKit::DetectBondStereoChemistry(*($self), conf);
  };

  void WedgeMolBonds(const RDKit::Conformer *conf) {
    RDKit::WedgeMolBonds(*($self), conf);
  };

  void pickBondsToWedge() {
    RDKit::Chirality::pickBondsToWedge(*($self));
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
                               bool permuteDeg4Nodes=false,
			       bool forceRDKit=false) {
    return RDDepict::compute2DCoords(*($self),
                               coordMap,
                               canonOrient,
                               clearConfs,
                               nFlipsPerSample,
                               nSamples,
                               sampleSeed,
				     permuteDeg4Nodes, forceRDKit);
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
                                           bool permuteDeg4Nodes=true, bool forceRDKit=false) {
    return RDDepict::compute2DCoordsMimicDistMat(*($self),
                                                 dmat,
                                                 canonOrient,
                                                 clearConfs,
                                                 weightDistMat,
                                                 nFlipsPerSample,
                                                 nSamples,
                                                 sampleSeed,
                                                 permuteDeg4Nodes, forceRDKit);

  }

  void generateDepictionMatching2DStructure(const RDKit::ROMol &reference,
                                          const RDKit::MatchVectType &refMatchVect,
                                          int confId=-1,
                                          bool forceRDKit=false) {
    RDDepict::generateDepictionMatching2DStructure(*($self),reference,refMatchVect,confId,forceRDKit);
  }
  RDKit::MatchVectType generateDepictionMatching2DStructure(const RDKit::ROMol &reference,
                                          int confId=-1,
                                          bool acceptFailure=false, bool forceRDKit=false,
                                          bool allowOptionalAttachments=false) {
    return RDDepict::generateDepictionMatching2DStructure(*($self),reference,confId,nullptr,
            acceptFailure,forceRDKit,allowOptionalAttachments);
  }
  RDKit::MatchVectType generateDepictionMatching2DStructure(const RDKit::ROMol &reference,
                                          int confId,
                                          const RDKit::ROMol &referencePattern,
                                          bool acceptFailure=false, bool forceRDKit=false,
                                          bool allowOptionalAttachments=false) {
    return RDDepict::generateDepictionMatching2DStructure(*($self),reference,confId,
           &referencePattern,acceptFailure,forceRDKit,allowOptionalAttachments);
  }

  double normalizeDepiction(int confId=-1, int canonicalize=1,
                            double scaleFactor=-1.0) {
    return RDDepict::normalizeDepiction(*($self), confId, canonicalize, scaleFactor);
  }

  void straightenDepiction(int confId=-1, bool minimizeRotation=false) {
    RDDepict::straightenDepiction(*($self), confId, minimizeRotation);
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
    auto atoms = ($self)->atoms();
    return new std::vector<RDKit::Atom*>(atoms.begin(), atoms.end());
  }

  std::vector<RDKit::Bond*> *getBonds() {
    auto bonds = ($self)->bonds();
    return new std::vector<RDKit::Bond*>(bonds.begin(), bonds.end());
  }

  std::vector<RDKit::Atom*> *getAtomNeighbors(RDKit::Atom *at) {
    auto atomNbrs = ($self)->atomNeighbors(at);
    return new std::vector<RDKit::Atom*>(atomNbrs.begin(), atomNbrs.end());
  }

  std::vector<RDKit::Bond*> *getAtomBonds(RDKit::Atom *at) {
    auto bondNbrs = ($self)->atomBonds(at);
    return new std::vector<RDKit::Bond*>(bondNbrs.begin(), bondNbrs.end());
  }

  std::vector<int> ToBinary(int propertyFlags=-1){
    std::string sres;
    if(propertyFlags>=0) {
      RDKit::MolPickler::pickleMol(*($self),sres,propertyFlags);
    } else {
      RDKit::MolPickler::pickleMol(*($self),sres);
    }
    std::vector<int> res(sres.length());
    std::copy(sres.begin(),sres.end(),res.begin());
    return res;
  };
  static RDKit::ROMOL_SPTR MolFromBinary(const std::vector<int> &pkl){
    std::string sres;
    sres.resize(pkl.size());
    std::copy(pkl.begin(),pkl.end(),sres.begin());
    RDKit::ROMol *res;
    try {
      res = new RDKit::ROMol(sres);
    } catch (const RDKit::MolPicklerException &e) {
      res = nullptr;
      throw;
    }
    return RDKit::ROMOL_SPTR(res);
  }
#ifdef SWIGJAVA
  const std::string toByteArray(int propertyFlags=-1) {
    std::string sres;
    if(propertyFlags>=0) {
      RDKit::MolPickler::pickleMol(*($self),sres,propertyFlags);
    } else {
      RDKit::MolPickler::pickleMol(*($self),sres);
    }
    return sres;
  }
#endif
#ifdef SWIGCSHARP
  const std::vector<unsigned char> toUCharVect(int propertyFlags=-1) {
    std::string sres;
    if(propertyFlags>=0) {
      RDKit::MolPickler::pickleMol(*($self),sres,propertyFlags);
    } else {
      RDKit::MolPickler::pickleMol(*($self),sres);
    }
    const std::vector<unsigned char> vec(sres.begin(), sres.end());
    return vec;
  }
#endif
  static RDKit::ROMOL_SPTR fromUCharVect(const std::vector<unsigned char> &pkl) {
    std::string sres(pkl.begin(), pkl.end());
    RDKit::ROMol *res;
    try {
      res = new RDKit::ROMol(sres);
    } catch (const RDKit::MolPicklerException &e) {
      res = nullptr;
      throw;
    }
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

  void alignMolConformers(const std::vector<unsigned int> *atomIds=0,
                          const std::vector<unsigned int> *confIds=0,
                          const RDNumeric::DoubleVector  *weights=0,
                          bool reflect=false, unsigned int maxIters=50) {
    RDKit::MolAlign::alignMolConformers(*($self), atomIds, confIds, weights, reflect, maxIters);
  }

  /* From GraphMol/MolAlign/AlignMolecules */
  double getAlignmentTransform(const RDKit::ROMol &refMol,
                             RDGeom::Transform3D &trans, int prbCid = -1,
                             int refCid = -1, const std::vector<std::pair<int,int> > *atomMap = 0,
                             const RDNumeric::DoubleVector *weights = 0,
                             bool reflect = false, unsigned int maxIters = 50){
    return RDKit::MolAlign::getAlignmentTransform(*($self), refMol, trans, prbCid, refCid, atomMap, weights, reflect, maxIters);
  }

  double getBestAlignmentTransform(
    const RDKit::ROMol &refMol, RDGeom::Transform3D &bestTrans,
    std::vector<std::pair<int,int> > &bestMatch, int prbCid = -1, int refCid = -1,
    const std::vector<std::vector<std::pair<int,int> > > &map = std::vector<std::vector<std::pair<int,int> > >(),
    int maxMatches = 1e6, bool symmetrizeConjugatedTerminalGroups = true,
    const RDNumeric::DoubleVector *weights = nullptr, bool reflect = false,
    unsigned int maxIters = 50) {
    return RDKit::MolAlign::getBestAlignmentTransform(*($self), refMol, bestTrans, bestMatch,
      prbCid, refCid, map, maxMatches, symmetrizeConjugatedTerminalGroups, weights, reflect, maxIters);
  }

  double calcRMS(
    const ROMol &refMol, int prbCid = -1, int refCid = -1,
    const std::vector<std::vector<std::pair<int,int> > > &map = std::vector<std::vector<std::pair<int,int> > >(),
    int maxMatches = 1e6, bool symmetrizeConjugatedTerminalGroups = true,
    const RDNumeric::DoubleVector *weights = nullptr) {
    return RDKit::MolAlign::CalcRMS(*($self), refMol, prbCid, refCid, map, maxMatches, symmetrizeConjugatedTerminalGroups, weights);
  }

  double getBestRMS(
    const ROMol &refMol, int prbCid = -1, int refCid = -1,
    const std::vector<std::vector<std::pair<int,int> > > &map = std::vector<std::vector<std::pair<int,int> > >(),
    int maxMatches = 1e6, bool symmetrizeConjugatedTerminalGroups = true,
    const RDNumeric::DoubleVector *weights = nullptr) {
    return RDKit::MolAlign::getBestRMS(*($self), refMol, prbCid, refCid, map, maxMatches, symmetrizeConjugatedTerminalGroups, weights);
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

  void computeGasteigerCharges(int nIter=12,bool throwOnParamFailure=false){
    RDKit::computeGasteigerCharges(*($self),nIter,throwOnParamFailure);
  }
  void computeGasteigerCharges(std::vector<double> &charges,
                               int nIter=12,bool throwOnParamFailure=false){
    RDKit::computeGasteigerCharges(*($self),charges,nIter,throwOnParamFailure);
  }

  /* From new_canon.h*/
  void rankMolAtoms(UINT_VECT &ranks,
                  bool breakTies = true, bool includeChirality = true,
                  bool includeIsotopes = true){
	RDKit::Canon::rankMolAtoms(*($self), ranks, breakTies, includeChirality, includeIsotopes);
  }
}




%extend RDKit::ROMol {
  std::string ToSVG(const std::vector<int> &highlightAtoms,
                    int lineWidthMult=2,int fontSize=50){
    // FIX: not sure any more what these are for
    if(fontSize<0) fontSize*=2;
    if(lineWidthMult<0) lineWidthMult *=2;
    std::stringstream outs;
    RDKit::MolDraw2DSVG drawer(300,300,outs);
    //drawer.setFontSize(static_cast<float>(fontSize)*drawer.fontSize()/50);
    drawer.drawMolecule(*($self),&highlightAtoms);
    drawer.finishDrawing();
    outs.flush();

    return outs.str();
  }
  std::string ToSVG(int lineWidthMult=2,int fontSize=50){
    // FIX: not sure any more what these are for
    if(fontSize<0) fontSize*=2;
    if(lineWidthMult<0) lineWidthMult *=2;
    std::stringstream outs;
    RDKit::MolDraw2DSVG drawer(300,300,outs);
    //drawer.setFontSize(static_cast<float>(fontSize)*drawer.fontSize()/50);
    drawer.drawMolecule(*($self));
    drawer.finishDrawing();
    outs.flush();

    return outs.str();

  }
}
