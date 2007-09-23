// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/types.h>
#include <map>

namespace RDKit{

#ifndef OLD_PICKLE
  //const int MolPickler::versionMajor=2;
  const int MolPickler::versionMajor=4;
#else
  const int MolPickler::versionMajor=1;
#endif
  const int MolPickler::versionMinor=0;
  const int MolPickler::versionPatch=0;
  const int MolPickler::endianId=0xDEADBEEF;

  template <typename T>
  void streamWrite(std::ostream &ss,MolPickler::Tags tag,const T &what){
    streamWrite(ss,tag);
    streamWrite(ss,what);
  };  



  void MolPickler::pickleMol(const ROMol *mol,std::ostream &ss){
    PRECONDITION(mol,"empty molecule");
    streamWrite(ss,endianId);
    streamWrite(ss,VERSION);
    streamWrite(ss,versionMajor);
    streamWrite(ss,versionMinor);
    streamWrite(ss,versionPatch);
#ifndef OLD_PICKLE
    if(mol->getNumAtoms()>255){
      _pickle<int>(mol,ss);
    } else {
      _pickle<unsigned char>(mol,ss);
    }
#else
    _pickleV1(mol,ss);
#endif    
  }
  void MolPickler::pickleMol(const ROMol *mol,std::string &res){
    PRECONDITION(mol,"empty molecule");
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    MolPickler::pickleMol(mol,ss);
    res = ss.str();
  }

  // NOTE: if the mol passed in here already has atoms and bonds, they will
  // be left intact.  The side effect is that ALL atom and bond bookmarks
  // will be blown out by the end of this process.
  void MolPickler::molFromPickle(std::istream &ss,ROMol *mol){
    PRECONDITION(mol,"empty molecule");
    Tags tag;
    int tmpInt;

    mol->clearAllAtomBookmarks();
    mol->clearAllBondBookmarks();
  
    streamRead(ss,tmpInt);
    if(tmpInt!=endianId){
      throw MolPicklerException("Bad pickle format: bad endian ID or invalid file format");
    }
  
    streamRead(ss,tag);
    if(tag!=VERSION){
      throw MolPicklerException("Bad pickle format: no version tag");
    }
    int majorVersion,minorVersion,patchVersion;
    streamRead(ss,majorVersion);
    streamRead(ss,minorVersion);
    streamRead(ss,patchVersion);
    if(majorVersion==1){
      _depickleV1(ss,mol);
    } else {
      // check to see how many atoms there are:
      //int pos=ss.tellg();
      int numAtoms;
      streamRead(ss,numAtoms);
      // move back in the stream:
      //ss.seekg(pos);
      if(numAtoms>255){
	_depickle<int>(ss,mol,majorVersion,numAtoms);
      } else {
	_depickle<unsigned char>(ss,mol,majorVersion,numAtoms);
      }
    }
    mol->clearAllAtomBookmarks();
    mol->clearAllBondBookmarks();
    if(majorVersion<4){
      // FIX for issue 220 - probably better to change the pickle format later
      MolOps::assignBondStereoCodes(*mol,true);
    }
  }
  void MolPickler::molFromPickle(const std::string &pickle,ROMol *mol){
    PRECONDITION(mol,"empty molecule");
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    ss.write(pickle.c_str(),pickle.length());
    MolPickler::molFromPickle(ss,mol);
  }


  //--------------------------------------
  //
  //            Molecules
  //
  //--------------------------------------
  template <typename T>
  void MolPickler::_pickle(const ROMol *mol,std::ostream &ss){
    PRECONDITION(mol,"empty molecule");
    int tmpInt;
    bool includeAtomCoords=true;
    std::map<int,int> atomIdxMap;

    tmpInt = static_cast<int>(mol->getNumAtoms());
    streamWrite(ss,tmpInt);
    tmpInt = static_cast<int>(mol->getNumBonds());
    streamWrite(ss,tmpInt);

    char flag = 0;
    if(includeAtomCoords) flag |= 0x1<<7;
    streamWrite(ss,flag);

    // -------------------
    //
    // Write Atoms
    //
    // -------------------
    streamWrite(ss,BEGINATOM);
    ROMol::ConstAtomIterator atIt;
    int nWritten=0;
    for(atIt=mol->beginAtoms();atIt!=mol->endAtoms();++atIt){
      _pickleAtom<T>(ss,*atIt);
      atomIdxMap[(*atIt)->getIdx()] = nWritten;
      nWritten++;
    }
  
    // -------------------
    //
    // Write Bonds
    //
    // -------------------
    streamWrite(ss,BEGINBOND);
    for(unsigned int i=0;i<mol->getNumBonds();i++){
      _pickleBond<T>(ss,mol->getBondWithIdx(i),atomIdxMap);
    }

    // -------------------
    //
    // Write Rings (if present)
    //
    // -------------------
    const RingInfo *ringInfo=mol->getRingInfo();
    if(ringInfo && ringInfo->isInitialized()){
      streamWrite(ss,BEGINSSSR);
      _pickleSSSR<T>(ss,ringInfo,atomIdxMap);
    }
    
    // pickle the conformations if necessary
    
    if (includeAtomCoords) {
      streamWrite(ss,BEGINCONFS);
      tmpInt = static_cast<int>(mol->getNumConformers());
      streamWrite(ss,tmpInt);
      
      ROMol::ConstConformerIterator ci;
      for (ci = mol->beginConformers(); ci != mol->endConformers(); ++ci) {
        const Conformer *conf = ci->get();
        _pickleConformer<T>(ss, conf);
      }
    }
    streamWrite(ss,ENDMOL);
  }

  template <typename T>
  void MolPickler::_depickle(std::istream &ss,ROMol *mol, int version,int numAtoms){
    PRECONDITION(mol,"empty molecule");
    bool directMap= mol->getNumAtoms()==0;
    Tags tag;
    int tmpInt;
    //int numAtoms,numBonds;
    int numBonds;
    
    streamRead(ss,tmpInt);
    numBonds = tmpInt;

    // did we include coordinates
    bool includeCoords=false;
    if (version >= 3) {
      char flag;
      streamRead(ss,flag);
      if (flag & 0x1<<7) includeCoords = true;
    }
    // -------------------
    //
    // Read Atoms
    //
    // -------------------
    streamRead(ss,tag);
    if(tag != BEGINATOM){
      throw MolPicklerException("Bad pickle format: BEGINATOM tag not found.");
    }
    Conformer *conf = 0;
    if ((version == 2) && includeCoords) {
      // there can only one conformation - since the poositions were stored on
      // the atoms themselves in this version
      conf = new Conformer(numAtoms);
      mol->addConformer(conf, true);
    }
    for(int i=0;i<numAtoms;i++){
      RDGeom::Point3D pos;
      Atom *atom=_addAtomFromPickle<T>(ss,mol,pos,version,directMap);
      if ((version == 2) && includeCoords) {
        // this is a older pickle so we go the pos
        conf->setAtomPos(i, pos);
      }
      if(!directMap){
	mol->setAtomBookmark(atom,i);
      }
    }

    // -------------------
    //
    // Read Bonds
    //
    // -------------------
    streamRead(ss,tag);
    if(tag != BEGINBOND){
      throw MolPicklerException("Bad pickle format: BEGINBOND tag not found.");
    }
    for(int i=0;i<numBonds;i++){
      Bond *bond=_addBondFromPickle<T>(ss,mol,version,directMap);
      if(!directMap){
	mol->setBondBookmark(bond,i);
      }
    }

    // -------------------
    //
    // Read Rings (if needed)
    //
    // -------------------
    streamRead(ss,tag);
    if(tag == BEGINSSSR){
      _addRingInfoFromPickle<T>(ss,mol,version,directMap);
      streamRead(ss,tag);
    }
    
    if (tag == BEGINCONFS) {
      //if ((version == 3) && includeCoords) {
      // read in the conformation
      streamRead(ss, tmpInt);
      int i;
      for (i = 0; i < tmpInt; i++) {
        Conformer *conf = _conformerFromPickle<T>(ss,version);
        mol->addConformer(conf);
      }
      streamRead(ss,tag);
    }

    if(tag != ENDMOL){
      throw MolPicklerException("Bad pickle format: ENDMOL tag not found.");
    }

  }



  //--------------------------------------
  //
  //            Atoms
  //
  //--------------------------------------

  // T refers to the type of the atom indices written
  template <typename T>
  void MolPickler::_pickleAtom(std::ostream &ss,const Atom *atom) {
    PRECONDITION(atom,"empty atom");
    // bool includeAtomCoords){ - coordinates are no longer included on the atom - we have conformations for that
    char tmpChar;
    signed char tmpSchar;
    char flags;

    tmpChar = atom->getAtomicNum()%128;
    streamWrite(ss,tmpChar);


    flags = 0;
    //if(includeAtomCoords) flags |= 0x1<<7;
    if(atom->getIsAromatic()) flags |= 0x1<<6;
    if(atom->getNoImplicit()) flags |= 0x1<<5;
    streamWrite(ss,flags);
    
    tmpSchar=static_cast<signed char>(atom->getMass() -
				      PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
    streamWrite(ss,tmpSchar);
    tmpSchar=static_cast<signed char>(atom->getFormalCharge());
    streamWrite(ss,tmpSchar);
    tmpChar = static_cast<char>(atom->getChiralTag());
    streamWrite(ss,tmpChar);
    tmpChar = static_cast<char>(atom->getHybridization());
    streamWrite(ss,tmpChar);
    tmpChar = static_cast<char>(atom->getNumExplicitHs());
    streamWrite(ss,tmpChar);
    tmpChar = static_cast<char>(atom->getExplicitValence());
    streamWrite(ss,tmpChar);
    tmpChar = static_cast<char>(atom->getImplicitValence());
    streamWrite(ss,tmpChar);
  }

  template <typename T>
  void MolPickler::_pickleConformer(std::ostream &ss,const Conformer *conf) {
    PRECONDITION(conf,"empty conformer");
    float tmpFloat;
    int tmpInt = static_cast<int>(conf->getId());
    streamWrite(ss,tmpInt);
    T tmpT = static_cast<T>(conf->getNumAtoms());
    streamWrite(ss,tmpT);
    const RDGeom::POINT3D_VECT &pts = conf->getPositions();
    for (RDGeom::POINT3D_VECT_CI pti = pts.begin(); pti != pts.end(); pti++) {
      tmpFloat = static_cast<float>(pti->x);
      streamWrite(ss,tmpFloat);
      tmpFloat = static_cast<float>(pti->y);
      streamWrite(ss,tmpFloat);
      tmpFloat = static_cast<float>(pti->z);
      streamWrite(ss,tmpFloat);
    } 
  }
    
  template <typename T> 
  Conformer *MolPickler::_conformerFromPickle(std::istream &ss,int version) {
    float tmpFloat;
    int tmpInt;
    streamRead(ss, tmpInt);
    unsigned int cid = static_cast<unsigned int>(tmpInt);
    T tmpT;
    streamRead(ss, tmpT);
    unsigned int numAtoms = static_cast<unsigned int>(tmpT);
    Conformer *conf = new Conformer(numAtoms);
    conf->setId(cid);
    for (unsigned int i = 0; i < numAtoms; i++) {
      streamRead(ss, tmpFloat);
      conf->getAtomPos(i).x = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat);
      conf->getAtomPos(i).y = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat);
      conf->getAtomPos(i).z = static_cast<double>(tmpFloat);
    }
    return conf;
  }
    
  template <typename T>
  Atom *MolPickler::_addAtomFromPickle(std::istream &ss,ROMol *mol,
                                       RDGeom::Point3D &pos,
				       int version,bool directMap){
    PRECONDITION(mol,"empty molecule");
    float x,y,z;
    char tmpChar;
    signed char tmpSchar;
    char flags;


    Atom *atom = new Atom();

    streamRead(ss,tmpChar);
    atom->setAtomicNum(static_cast<int>(tmpChar));

    streamRead(ss,flags);
    atom->setIsAromatic(flags & 0x1<<6);
    atom->setNoImplicit(flags & 0x1<<5);
    
    // are coordinates present?
    if(flags & 0x1<<7){
      streamRead(ss,x);
      pos.x = static_cast<double>(x);
      streamRead(ss,y);
      pos.y = static_cast<double>(y);
      streamRead(ss,z);
      pos.z = static_cast<double>(z);
      //atom->setPos(x,y,z);
    }
    
    streamRead(ss,tmpSchar);
    atom->setMass(PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum())+
		  static_cast<int>(tmpSchar));

    streamRead(ss,tmpSchar);
    atom->setFormalCharge(static_cast<int>(tmpSchar));

    streamRead(ss,tmpChar);
    atom->setChiralTag(static_cast<Atom::ChiralType>(tmpChar));
    streamRead(ss,tmpChar);
    atom->setHybridization(static_cast<Atom::HybridizationType>(tmpChar));
    streamRead(ss,tmpChar);
    atom->setNumExplicitHs(static_cast<int>(tmpChar));
    streamRead(ss,tmpChar);
    atom->d_explicitValence  = tmpChar;
    streamRead(ss,tmpChar);
    atom->d_implicitValence  = tmpChar;
    
    mol->addAtom(atom,false,true);
    return atom;
  }

  //--------------------------------------
  //
  //            Bonds
  //
  //--------------------------------------

  template <typename T>
  void MolPickler::_pickleBond(std::ostream &ss,const Bond *bond,
			       std::map<int,int> &atomIdxMap){
    PRECONDITION(bond,"empty bond");
    T tmpT;
    char tmpChar;
    char flags;

    tmpT = static_cast<T>(atomIdxMap[bond->getBeginAtomIdx()]);
    streamWrite(ss,tmpT);
    tmpT = static_cast<T>(atomIdxMap[bond->getEndAtomIdx()]);
    streamWrite(ss,tmpT);

    flags = 0;
    if(bond->getIsAromatic())	flags |= 0x1<<6;
    if(bond->getIsConjugated()) flags |= 0x1<<5;
    streamWrite(ss,flags);
    tmpChar = static_cast<char>(bond->getBondType());
    streamWrite(ss,tmpChar);
    tmpChar = static_cast<char>(bond->getBondDir());
    streamWrite(ss,tmpChar);

    // write info about the stereochemistry:
    tmpChar = static_cast<char>(bond->getStereo());
    streamWrite(ss,tmpChar);
    if(bond->getStereo()!=Bond::STEREONONE){
      const INT_VECT &stereoAts=bond->getStereoAtoms();
      tmpChar = stereoAts.size();
      streamWrite(ss,tmpChar);
      for(INT_VECT_CI idxIt=stereoAts.begin();idxIt!=stereoAts.end();++idxIt){
	tmpT = static_cast<T>(*idxIt);
	streamWrite(ss,tmpT);
      }
    }    
  }

  template <typename T>
  Bond *MolPickler::_addBondFromPickle(std::istream &ss,ROMol *mol,int version,
				       bool directMap){
    PRECONDITION(mol,"empty molecule");
    char tmpChar;
    char flags;
    T tmpT;

    Bond *bond = new Bond();
    streamRead(ss,tmpT);
    if(directMap){
      bond->setBeginAtomIdx(static_cast<int>(tmpT));
    } else {
      bond->setBeginAtomIdx(mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx());
    }
    streamRead(ss,tmpT);
    if(directMap){
      bond->setEndAtomIdx(static_cast<int>(tmpT));
    } else {
      bond->setEndAtomIdx(mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx());

    }
    streamRead(ss,flags);
    bond->setIsAromatic(flags & 0x1<<6);
    bond->setIsConjugated(flags & 0x1<<5);

    streamRead(ss,tmpChar);
    bond->setBondType(static_cast<Bond::BondType>(tmpChar));
    streamRead(ss,tmpChar);
    bond->setBondDir(static_cast<Bond::BondDir>(tmpChar));

    if(version>3){
      streamRead(ss,tmpChar);
      Bond::BondStereo stereo=static_cast<Bond::BondStereo>(tmpChar);
      bond->setStereo(stereo);
      if(stereo!=Bond::STEREONONE){
	streamRead(ss,tmpChar);
	for(char i=0;i<tmpChar;++i){
	  streamRead(ss,tmpT);
	  bond->getStereoAtoms().push_back(static_cast<int>(tmpT));
	}
      }
    }
    mol->addBond(bond,true);
    return bond;
  }

  //--------------------------------------
  //
  //            Rings
  //
  //--------------------------------------
  template <typename T>
  void MolPickler::_pickleSSSR(std::ostream &ss,const RingInfo *ringInfo,
			       std::map<int,int> &atomIdxMap){
    PRECONDITION(ringInfo,"missing ring info");
    T tmpT;
    tmpT = ringInfo->numRings();
    streamWrite(ss,tmpT);
    for(unsigned int i=0;i<ringInfo->numRings();i++){
      INT_VECT ring = ringInfo->atomRings()[i];
      tmpT = static_cast<T>(ring.size());
      streamWrite(ss,tmpT);

      for(unsigned int j=0;j<ring.size();j++){
	tmpT = static_cast<T>(atomIdxMap[ring[j]]);
	streamWrite(ss,tmpT);
      }
      ring = ringInfo->bondRings()[i];
      for(unsigned int j=0;j<ring.size();j++){
	tmpT = static_cast<T>(ring[j]);
	streamWrite(ss,tmpT);
      }
    }
  }

  template <typename T>
  void MolPickler::_addRingInfoFromPickle(std::istream &ss,ROMol *mol,
					  int version,
					  bool directMap){
    PRECONDITION(mol,"empty molecule");
    RingInfo *ringInfo=mol->getRingInfo();
    if(!ringInfo->isInitialized()) ringInfo->initialize();

    T numRings;
    streamRead(ss,numRings);

    if(numRings>0){
      ringInfo->preallocate(mol->getNumAtoms(),mol->getNumBonds());
    }
    
    for(int i=0;i<numRings;i++){
      T tmpT;
      T ringSize;
      streamRead(ss,ringSize);

      INT_VECT atoms(static_cast<int>(ringSize));
      INT_VECT bonds(static_cast<int>(ringSize));
      for(int j=0; j<ringSize; j++){
	streamRead(ss,tmpT);
	if(directMap){
	  atoms[j] = static_cast<int>(tmpT);
	} else {
	  atoms[j] = mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx();
	}
      }
      for(int j=0; j<ringSize; j++){
	streamRead(ss,tmpT);
	if(directMap){
	  bonds[j] = static_cast<int>(tmpT);
	} else {
	  bonds[j] = mol->getBondWithBookmark(static_cast<int>(tmpT))->getIdx();
	}
      }
      ringInfo->addRing(atoms,bonds);
    }
  }


  //--------------------------------------
  //
  //            Version 1 Pickler:
  //
  //--------------------------------------

  void MolPickler::_pickleV1(const ROMol *mol,std::ostream &ss){
    PRECONDITION(mol,"empty molecule");
    ROMol::ConstAtomIterator atIt;
    const Conformer *conf = 0;
    if (mol->getNumConformers() > 0) {
      conf = &(mol->getConformer());
    }
    for(atIt=mol->beginAtoms();atIt!=mol->endAtoms();atIt++){
      const Atom *atom = *atIt;

      streamWrite(ss,BEGINATOM);
      streamWrite(ss,ATOM_NUMBER,atom->getAtomicNum());

      streamWrite(ss,ATOM_INDEX,atom->getIdx());

      streamWrite(ss,ATOM_POS);
      RDGeom::Point3D p;
      if (conf) {
        p = conf->getAtomPos(atom->getIdx());
      } 
      streamWrite(ss,p.x);
      streamWrite(ss,p.y);
      streamWrite(ss,p.z);
    

      if(atom->getFormalCharge() != 0){
	streamWrite(ss,ATOM_CHARGE,atom->getFormalCharge());
      }
      if(atom->getNumExplicitHs() != 0){
	streamWrite(ss,ATOM_NEXPLICIT,atom->getNumExplicitHs());
      }
      if(atom->getChiralTag() != 0){
	streamWrite(ss,ATOM_CHIRALTAG,atom->getChiralTag());
      }
      if(atom->getMass() != 0.0){
	streamWrite(ss,ATOM_MASS,atom->getMass());
      }
      if(atom->getIsAromatic()){
	streamWrite(ss,ATOM_ISAROMATIC,static_cast<char>(atom->getIsAromatic()));
      }
      streamWrite(ss,ENDATOM);
    }
  
    ROMol::ConstBondIterator bondIt;
    for(bondIt=mol->beginBonds();bondIt!=mol->endBonds();bondIt++){
      const Bond *bond = *bondIt;
      streamWrite(ss,BEGINBOND);
      streamWrite(ss,BOND_INDEX,bond->getIdx());
      streamWrite(ss,BOND_BEGATOMIDX,bond->getBeginAtomIdx());
      streamWrite(ss,BOND_ENDATOMIDX,bond->getEndAtomIdx());
      streamWrite(ss,BOND_TYPE,bond->getBondType());
      if(bond->getBondDir()){
	streamWrite(ss,BOND_DIR,bond->getBondDir());
      }
      streamWrite(ss,ENDBOND);
    }
    streamWrite(ss,ENDMOL);
  }

  void MolPickler::_depickleV1(std::istream &ss,ROMol *mol){
    PRECONDITION(mol,"empty molecule");
    Tags tag;

    Conformer *conf = new Conformer();
    mol->addConformer(conf);
    streamRead(ss,tag);
    while(tag != ENDMOL){
      switch(tag){
      case BEGINATOM:
	_addAtomFromPickleV1(ss,mol);
	break;
      case BEGINBOND:
	_addBondFromPickleV1(ss,mol);
	break;
      default:
	UNDER_CONSTRUCTION("bad tag in pickle");
      }
      streamRead(ss,tag);
    }
    mol->clearAllAtomBookmarks();
    mol->clearAllBondBookmarks();
  }
  
  void MolPickler::_addAtomFromPickleV1(std::istream &ss,ROMol *mol) {
    PRECONDITION(mol,"empty molecule");
    Tags tag;
    int intVar;
    double dblVar;
    char charVar;
    streamRead(ss,tag);
    Atom *atom = new Atom();
    Conformer &conf = mol->getConformer();
    RDGeom::Point3D pos;
    while(tag != ENDATOM){
      switch(tag){
      case ATOM_INDEX:
	streamRead(ss,intVar);
	mol->setAtomBookmark(atom,intVar);
	break;
      case ATOM_NUMBER:
	streamRead(ss,intVar);
	atom->setAtomicNum(intVar);
	break;
      case ATOM_POS:
	streamRead(ss,pos.x);
	streamRead(ss,pos.y);
	streamRead(ss,pos.z);
        break;
      case ATOM_CHARGE:
	streamRead(ss,intVar);
	atom->setFormalCharge(intVar);
	break;
      case ATOM_NEXPLICIT:
	streamRead(ss,intVar);
	atom->setNumExplicitHs(intVar);
	break;
      case ATOM_CHIRALTAG:
	streamRead(ss,intVar);
	atom->setChiralTag(static_cast<Atom::ChiralType>(intVar));
	break;
      case ATOM_MASS:
	streamRead(ss,dblVar);
	atom->setMass(dblVar);
	break;
      case ATOM_ISAROMATIC:
	streamRead(ss,charVar);
	atom->setIsAromatic(charVar);
	break;
      default:
	ASSERT_INVARIANT(0,"bad tag in atom block of pickle");
      }
      streamRead(ss,tag);
    }
    unsigned int id = mol->addAtom(atom,false,true);
    conf.setAtomPos(id, pos);
  }
  void MolPickler::_addBondFromPickleV1(std::istream &ss,ROMol *mol){
    PRECONDITION(mol,"empty molecule");
    Tags tag;
    int intVar,idx=-1;
    Bond::BondType bt;
    Bond::BondDir bd;
    streamRead(ss,tag);
    Bond *bond = new Bond();
    while(tag != ENDBOND){
      switch(tag){
      case BOND_INDEX:
	streamRead(ss,idx);
	break;
      case BOND_BEGATOMIDX:
	streamRead(ss,intVar);
	bond->setBeginAtomIdx(mol->getAtomWithBookmark(intVar)->getIdx());
	break;
      case BOND_ENDATOMIDX:
	streamRead(ss,intVar);
	bond->setEndAtomIdx(mol->getAtomWithBookmark(intVar)->getIdx());
	break;
      case BOND_TYPE:
	streamRead(ss,bt);
	bond->setBondType(bt);
	break;
      case BOND_DIR:
	streamRead(ss,bd);
	bond->setBondDir(bd);
	break;
      default:
	ASSERT_INVARIANT(0,"bad tag in bond block of pickle");
      }
      streamRead(ss,tag);
    }
    mol->addBond(bond,true);
  }


};
