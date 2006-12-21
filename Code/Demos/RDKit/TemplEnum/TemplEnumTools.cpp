//  $Id$
//
//  Copyright (C) 2003,2004 Rational Discovery LLC
//   All Rights Reserved
//
#include "TemplEnum.h"
#include <Geometry/Transform3D.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/FileParsers.h>
#define FEQ(_a_,_b_) (fabs((_a_)-(_b_))<1e-4)

namespace TemplateEnum {
  using namespace RDKit;

  // ------------------------------------------------------------------
  //
  // transforms a sidechain so that it is oriented better for attachment
  // to a molecule
  //
  //  Arguments:
  //   mol: core molecule 
  //   sidechain: sidechain to attach.
  //   molConnectorIdx: the index of the attachment point atom in the
  //     molecule.
  //   sidechainConnectorIdx: the index of the attachment point atom
  //     in the sidechain.
  // 
  // ------------------------------------------------------------------
  void orientSidechain(RWMol *mol,RWMol *sidechain,
		       int molAttachIdx,int sidechainAttachIdx){
    PRECONDITION(mol,"bad molecule");
    PRECONDITION(sidechain,"bad molecule");

    // ---------
    // start by getting our 4 atoms
    // ---------
    Conformer &molConf=mol->getConformer();
    Conformer &sidechainConf=sidechain->getConformer();
    Atom *molConnAtom,*molAttachAtom;
    Atom *chainConnAtom,*chainAttachAtom;
    molAttachAtom = mol->getAtomWithIdx(molAttachIdx);
    chainAttachAtom = sidechain->getAtomWithIdx(sidechainAttachIdx);
    PRECONDITION(molAttachAtom->getDegree()==1,"attachment points must be degree 1");
    PRECONDITION(chainAttachAtom->getDegree()==1,"attachment points must be degree 1");
    RWMol::ADJ_ITER nbrIdx,endNbrs;
    boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(molAttachAtom);
    molConnAtom=mol->getAtomWithIdx(*nbrIdx);
    boost::tie(nbrIdx,endNbrs) = sidechain->getAtomNeighbors(chainAttachAtom);
    chainConnAtom=sidechain->getAtomWithIdx(*nbrIdx);

    //-----------------------------------------
    //  Notation:
    //    Pmc: molecule connection point (the atom that will be
    //     removed from the molecule).
    //    Pma: molecule attachment point (the atom to which we'll form
    //     the bond).
    //    Psc: sidechain connection point
    //    Psa: sidechain attachment point
    //    Vm: Pmc-Pma (molecular attachment vector)
    //    Vs: Psc-Psa (sidechain attachment vector)
    //
    //-----------------------------------------
    RDGeom::Transform3D sidechainTform,templateTform,tmpTform;

    RDGeom::Point3D Vm,Um,Pmc,Pma;
    RDGeom::Point3D Vs,Us,Psc,Psa;
    Pmc = molConf.getAtomPos(molConnAtom->getIdx());
    Pma = molConf.getAtomPos(molAttachAtom->getIdx());
    std::cerr << "p=array(["<<Pma.x<<","<<Pma.y<<","<<Pma.z<<"])" << std::endl;
    Psc = sidechainConf.getAtomPos(chainConnAtom->getIdx());
    Psa = sidechainConf.getAtomPos(chainAttachAtom->getIdx());

    templateTform.setToIdentity();
    Vm = Pmc - Pma;
    // note the opposite direction here:
    Vs = Psa - Psc;
    Um = Vm;
    Um.normalize();
    std::cerr << "Um=array(["<<Um.x<<","<<Um.y<<","<<Um.z<<"])" << std::endl;
    Us = Vs;
    Us.normalize();
    std::cerr << "Us=array(["<<Us.x<<","<<Us.y<<","<<Us.z<<"])" << std::endl;
    // translate Psc -> Pma
    //RDGeom::Point3D headTrans = Pma-Psc;
    
    templateTform.setToIdentity();

    tmpTform.setToIdentity();
    tmpTform.SetTranslation(Pma);
    templateTform *= tmpTform;

    double sinT,cosT;
    cosT = Us.dotProduct(Um);
    if(cosT>1.0) cosT = 1.0;
    if(fabs(cosT)<1.0){
      tmpTform.setToIdentity();
      sinT = sqrt(1.0-cosT*cosT);
      RDGeom::Point3D rotnAxis=Us.crossProduct(Um);
      rotnAxis.normalize();
      std::cerr << "ax=array(["<<rotnAxis.x<<","<<rotnAxis.y<<","<<rotnAxis.z<<"])" << std::endl;
      tmpTform.SetRotation(cosT,sinT,rotnAxis);
      templateTform *= tmpTform;
    } else if(cosT==1.0){
      RDGeom::Point3D normal(1,0,0);
      if(fabs(Us.dotProduct(normal))==1.0){
	normal = RDGeom::Point3D(0,1,0);
      }
      RDGeom::Point3D rotnAxis=Us.crossProduct(normal);
      templateTform.SetRotation(-1,0,rotnAxis);
    }
    
    tmpTform.setToIdentity();
    tmpTform.SetTranslation(Psc*-1.0);
    templateTform *= tmpTform;

    // ---------
    // transform the atomic positions in the sidechain:
    // ---------
    MolTransforms::transformMolsAtoms(sidechain,templateTform);


    // that's it!
  }
  // ------------------------------------------------------------------
  //
  // attaches a sidechain fragment to a molecule.
  //
  //  Arguments:
  //   mol: molecule to be modified
  //   sidechain: sidechain to attach.  The sidechain is copied in.
  //   molConnectorIdx: the index of the attachment point atom in the
  //     molecule.
  //   sidechainConnectorIdx: the index of the attachment point atom
  //     in the sidechain.
  //   bondType: type of the bond to form between the atoms
  //
  //  The connector atoms are *NOT* part of the final molecule, they
  //  merely serve to establish where things connect.
  // 
  // ------------------------------------------------------------------
  void molAddSidechain(RWMol *mol,RWMol *sidechain,
		       int molConnectorIdx,int sidechainConnectorIdx,
		       Bond::BondType bondType){
    PRECONDITION(mol,"bad molecule provided");
    PRECONDITION(sidechain,"bad sidechain provided");
    int origNumAtoms=mol->getNumAtoms();
    mol->insertMol(*sidechain);

    // get pointers to the two connectors (these are the atoms
    // we'll end up removing)
    Atom *molConnAtom,*sidechainConnAtom;
    molConnAtom = mol->getAtomWithIdx(molConnectorIdx);
    sidechainConnAtom = mol->getAtomWithIdx(sidechainConnectorIdx+origNumAtoms);

    // now use those pointers to get the atoms which will remain
    // (we're going to connect these) and remove the original
    // connection points from the molecule
    Atom *tmpAtom;
    RWMol::ADJ_ITER nbrIdx,endNbrs;
    boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(molConnAtom);

    // we are assuming that the first neighbor is the correct one.
    // Really we should be able to assume that there is only a single
    // attachment point.
    tmpAtom = molConnAtom;
    molConnAtom = mol->getAtomWithIdx(*nbrIdx);
    mol->removeAtom(tmpAtom);
    // repeat that process for the sidechain:
    boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(sidechainConnAtom);
    tmpAtom = sidechainConnAtom;
    sidechainConnAtom = mol->getAtomWithIdx(*nbrIdx);
    mol->removeAtom(tmpAtom);

    // finally connect the remaining atoms:
    mol->addBond(molConnAtom,sidechainConnAtom,bondType);
  }


  // used as a starting point for connection bookmarks
  const int CONNECT_BOOKMARK_START=0x23424223;

  // ------------------------------------------------------------------
  //
  // Loop through all the atoms and bookmark the attachment points
  //
  //  attachment points are assumed to have one of these properties:
  //   - "molFileAlias" (set by the mol file parser)
  //   - "dummyLabel"   (set by the SMILES parser)
  //  frontMarker is the "recognition character" to be used to pick
  //  valid labels.  e.g. if the frontMarker is 'X', a label beginning
  //  with 'Y' will not be marked.
  //
  //  In addition to bookmarking the attachment points, we also set
  //  the "maxAttachIdx" property, which holds an integer with the
  //  maximum attachment bookmark index. This is used in the
  //  enumeration to prevent us from having to scan through all the
  //  molecule's bookmarks
  //
  //
  // ------------------------------------------------------------------
  void markAttachmentPoints(RWMOL_SPTR mol,char frontMarker){
    markAttachmentPoints(mol.get(),frontMarker);
  }
  void markAttachmentPoints(RWMol *mol,char frontMarker){
    PRECONDITION(mol,"bad molecule");

    RWMol::AtomIterator atomIt;
    int maxAttachIdx=0;

    // scan through the atoms and mark those that have aliases
    //  (these might be attachment points)
    for(atomIt=mol->beginAtoms();atomIt!=mol->endAtoms();atomIt++){
      // start by finding possible attachment point properties:
      std::string attachLabel="";
      if((*atomIt)->hasProp("molFileAlias")){
	(*atomIt)->getProp("molFileAlias",attachLabel);
      } else if((*atomIt)->hasProp("dummyLabel")){
	(*atomIt)->getProp("dummyLabel",attachLabel);
      }

      // if we got one and it starts with the appropriate front
      // marker, proceed:
      if(attachLabel!="" && attachLabel[0]==frontMarker ){
	// to avoid trouble later, we guarantee that the attachment
	// point has degree 1 (one bond to it).
	if((*atomIt)->getDegree()>1)
	  throw EnumException("More than one bond to an attachment point.");

	int offset = CONNECT_BOOKMARK_START;
	if(attachLabel.length()>1){
	  if(attachLabel[1]>='a' && attachLabel[1]<='z'){
	    offset += (int)attachLabel[1] - (int)'a';
	  }
	}
	mol->setAtomBookmark(*atomIt,offset);
	if(offset>maxAttachIdx) maxAttachIdx=offset;
      }
    }
    if(maxAttachIdx){
      mol->setProp("maxAttachIdx",maxAttachIdx);
    }
  }

  // ------------------------------------------------------------------
  //
  // loops through the sidechain molecules and calls
  // _markAttachmentPoints()_ on each.
  //
  // ------------------------------------------------------------------
  void prepareSidechains(RWMOL_SPTR_VECT *sidechains,char frontMarker){
    PRECONDITION(sidechains,"bad sidechain list");
    RWMOL_SPTR_VECT::iterator mpvI;
    for(mpvI=sidechains->begin();mpvI!=sidechains->end();mpvI++){
      markAttachmentPoints(*mpvI,frontMarker);
    }
  }

  // ------------------------------------------------------------------
  //
  // Enumerates the library around a template and returns the result.
  //
  //
  // ------------------------------------------------------------------
  RWMOL_SPTR_VECT enumerateLibrary(RWMol *templateMol,
				   VECT_RWMOL_SPTR_VECT &sidechains,
				   bool orientSidechains){
    PRECONDITION(templateMol,"bad molecule");
    RWMOL_SPTR_VECT res,tmp;
    res.push_back(RWMOL_SPTR(new RWMol(*templateMol)));

    // if there's no attachment point on the molecule or no
    // sidechains, return now: 
    if(!templateMol->hasProp("maxAttachIdx") || sidechains.size()==0 )
      return res;

    int maxIdx;
    templateMol->getProp("maxAttachIdx",maxIdx);

    tmp.clear();
    // loop over the sidechains and attach them
    for(unsigned int i=0;i<sidechains.size();i++){
      int tgtMark=CONNECT_BOOKMARK_START+i;
      // here's another boundary condition
      if(tgtMark>maxIdx) break;

      /// loop over all atoms with the appropriate mark
      //  This means that if a mol has two attachment points with the
      //  same name (e.g. two Xa's) they'll always have the same
      //  sidechain attached to them.  This is a feature.
      RWMOL_SPTR_VECT::iterator sidechainIt;
      for(sidechainIt=sidechains[i].begin();
	  sidechainIt!=sidechains[i].end();
	  sidechainIt++){
	// we've got our sidechain, find the atom it attaches from
	if( (*sidechainIt)->hasAtomBookmark(CONNECT_BOOKMARK_START) ){
	  //
	  // NOTE: If there's more than one marked atom in the sidechain,
	  ///      we'll only use the first for the moment.
	  //
	  int sidechainAtomIdx = (*sidechainIt)->getAtomWithBookmark(CONNECT_BOOKMARK_START)->getIdx();
	
	  // now add the sidechain to each molecule
	  RWMOL_SPTR_VECT::iterator templMolIt;
	  // loop over all the mols we've generated to this point
	  for(templMolIt=res.begin();templMolIt!=res.end();templMolIt++){
	    RWMol *templ =  new RWMol(**templMolIt);
	    std::string name,tmpStr;
	    if(templ->hasProp("_Name")){
	      templ->getProp("_Name",tmpStr);
	      name = name + " " + tmpStr;
	    }
	    while(templ->hasAtomBookmark(tgtMark)){
	      // this is the atom we'll be replacing in the template
	      Atom *at = templ->getAtomWithBookmark(tgtMark);

	      // copy and transform the sidechain:
	      RWMol *sidechain;
	      if(orientSidechains){
		sidechain = new RWMol(*(sidechainIt->get()));
		orientSidechain(templ,sidechain,
				at->getIdx(),sidechainAtomIdx);
	      } else {
		sidechain = sidechainIt->get();
	      }
	      // FIX: need to use the actual bond order here:
	      molAddSidechain(templ,sidechain,
			      at->getIdx(),sidechainAtomIdx,
			      Bond::SINGLE);
	      if(sidechain->hasProp("_Name")){
		sidechain->getProp("_Name",tmpStr);
		name = name + " " + tmpStr;
	      }
	      templ->clearAtomBookmark(tgtMark,at);
	      if(orientSidechains){
		delete sidechain;
	      }
	    }
	    //std::cout << templ << "> " << MolToSmiles(*templ) << std::endl;
	    if(name != "") templ->setProp("_Name",name);
	    tmp.push_back(RWMOL_SPTR(templ));
	  }
	}
      }

      //
      // if we just made any molecules, free up the memory used by the
      // existing result set and move the molecules we just generated
      // over
      if(tmp.size()){
#if 0      
	RWMOL_SPTR_VECT::iterator tmpMolIt;
	for(tmpMolIt=res.begin();tmpMolIt!=res.end();tmpMolIt++){
	  delete *tmpMolIt;
	}
#endif
	res = tmp;
	tmp.clear();
      }
    }
    return res;
  }

  // ------------------------------------------------------------------
  //
  //  Reads a template and library of sidechains from input files.
  //   the template file should be a mol file and the sidechain files
  //   SD files
  //
  // ------------------------------------------------------------------
  RWMOL_SPTR_VECT enumFromFiles(const char *templateName,
			    std::vector<const char *> &sidechainNames){
    PRECONDITION(templateName,"bad template file name passed in");

    // build and mark the template molecule
    RWMol *templ = MolFileToMol(templateName,false);
    if(!templ) throw EnumException("could not construct template molecule"); 
    markAttachmentPoints(templ,'X');

    // now build and mark each set of sidechains:
    RWMOL_SPTR_VECT sidechains;
    VECT_RWMOL_SPTR_VECT allSidechains;
    for(std::vector<const char*>::const_iterator i=sidechainNames.begin();
	i!=sidechainNames.end();i++){
      sidechains = SDFileToMols(*i,false);
      if(!sidechains.size()){
	std::string err="no sidechains read from file: ";
	err += *i;
	throw EnumException(err.c_str());
      }
      prepareSidechains(&sidechains,'X');
      allSidechains.push_back(sidechains);
    }

    // enumerate the library:
    RWMOL_SPTR_VECT library=enumerateLibrary(templ,allSidechains);


    //--------------------------
    //
    // Clean up the molecules and sidechains we constructed along the
    // way. 
    //
    //--------------------------
    delete templ;
#if 0
    VECT_RWMOL_SPTR_VECT::iterator vmpvI;
    for(vmpvI=allSidechains.begin();vmpvI!=allSidechains.end();vmpvI++){
      RWMOL_SPTR_VECT::iterator mpvI;
      for(mpvI=vmpvI->begin();mpvI!=vmpvI->end();mpvI++){
	delete *mpvI;
      }
      vmpvI->clear();
    }
#endif
    allSidechains.clear();

    return library;
  }  

} // end of TemplateEnum namespace

