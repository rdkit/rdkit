// $Id: MolPickler.cpp 1123 2009-06-01 13:04:33Z glandrum $
//
//  Copyright (C) 2009 Greg Landrum
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <sstream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/types.h>
#include <boost/cstdint.hpp>
using boost::int32_t;
using boost::uint32_t;

namespace RDKit{
  const int32_t ReactionPickler::versionMajor=1;
  const int32_t ReactionPickler::versionMinor=1;
  const int32_t ReactionPickler::versionPatch=0;
  const int32_t ReactionPickler::endianId=0xDEADBEEF;

  void streamWrite(std::ostream &ss,ReactionPickler::Tags tag){
    int32_t tmp=tag;
    streamWrite(ss,tmp);
  }
  template <typename T>
  void streamWrite(std::ostream &ss,ReactionPickler::Tags tag,const T &what){
    streamWrite(ss,tag);
    streamWrite(ss,what);
  };  

  void streamRead(std::istream &ss,ReactionPickler::Tags &tag){
    int32_t tmp;
    streamRead(ss,tmp);
    tag=static_cast<ReactionPickler::Tags>(tmp);
  }

  void ReactionPickler::pickleReaction(const ChemicalReaction *rxn,std::ostream &ss){
    PRECONDITION(rxn,"empty reaction");
    streamWrite(ss,endianId);
    streamWrite(ss,VERSION);
    streamWrite(ss,versionMajor);
    streamWrite(ss,versionMinor);
    streamWrite(ss,versionPatch);
    _pickle(rxn,ss);
  }
  void ReactionPickler::pickleReaction(const ChemicalReaction *rxn,std::string &res){
    PRECONDITION(rxn,"empty reaction");
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    ReactionPickler::pickleReaction(rxn,ss);
    res = ss.str();
  }

  // NOTE: if the reaction passed in here already has reactants and products
  // be left intact. 
  void ReactionPickler::reactionFromPickle(std::istream &ss,ChemicalReaction *rxn){
    PRECONDITION(rxn,"empty reaction");
    Tags tag;
    int32_t tmpInt;

    streamRead(ss,tmpInt);
    if(tmpInt!=endianId){
      throw ReactionPicklerException("Bad pickle format: bad endian ID or invalid file format");
    }
  
    streamRead(ss,tag);
    if(tag!=VERSION){
      throw ReactionPicklerException("Bad pickle format: no version tag");
    }
    int32_t majorVersion,minorVersion,patchVersion;
    streamRead(ss,majorVersion);
    streamRead(ss,minorVersion);
    streamRead(ss,patchVersion);
    if(majorVersion>versionMajor||(majorVersion==versionMajor&&minorVersion>versionMinor)){
      BOOST_LOG(rdWarningLog)<<"Depickling from a version number ("<<majorVersion<<"." << minorVersion<<")" << "that is higher than our version ("<<versionMajor<<"."<<versionMinor<<").\nThis probably won't work."<<std::endl;
    }
    majorVersion=1000*majorVersion+minorVersion*10+patchVersion;

    _depickle(ss,rxn,majorVersion);
  }
  void ReactionPickler::reactionFromPickle(const std::string &pickle,ChemicalReaction *rxn){
    PRECONDITION(rxn,"empty reaction");
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    ss.write(pickle.c_str(),pickle.length());
    ReactionPickler::reactionFromPickle(ss,rxn);
  }

  void ReactionPickler::_pickle(const ChemicalReaction *rxn,std::ostream &ss){
    PRECONDITION(rxn,"empty reaction");
    uint32_t tmpInt;

    tmpInt = static_cast<int32_t>(rxn->getNumReactantTemplates());
    streamWrite(ss,tmpInt);
    tmpInt = static_cast<int32_t>(rxn->getNumProductTemplates());
    streamWrite(ss,tmpInt);
    tmpInt = static_cast<int32_t>(rxn->getNumAgentTemplates());
    streamWrite(ss,tmpInt);

    uint32_t flag = 0;
    if(rxn->getImplicitPropertiesFlag()) flag |= 0x1;
    if(rxn->df_needsInit) flag |= 0x2;
    streamWrite(ss,flag);

    // -------------------
    //
    // Write Reactants
    //
    // -------------------
    streamWrite(ss,BEGINREACTANTS);
    for(MOL_SPTR_VECT::const_iterator tmpl=rxn->beginReactantTemplates();
        tmpl!=rxn->endReactantTemplates();++tmpl){
      MolPickler::pickleMol(tmpl->get(),ss);
    }
    streamWrite(ss,ENDREACTANTS);
  
    streamWrite(ss,BEGINPRODUCTS);
    for(MOL_SPTR_VECT::const_iterator tmpl=rxn->beginProductTemplates();
        tmpl!=rxn->endProductTemplates();++tmpl){
      MolPickler::pickleMol(tmpl->get(),ss);
    }
    streamWrite(ss,ENDPRODUCTS);

    if(rxn->getNumAgentTemplates()){
      streamWrite(ss,BEGINAGENTS);
      for(MOL_SPTR_VECT::const_iterator tmpl=rxn->beginAgentTemplates();
          tmpl!=rxn->endAgentTemplates();++tmpl){
        MolPickler::pickleMol(tmpl->get(),ss);
      }
      streamWrite(ss,ENDAGENTS);
    }
    streamWrite(ss,ENDREACTION);
  } // end of _pickle

  void ReactionPickler::_depickle(std::istream &ss,ChemicalReaction *rxn, int version){
    PRECONDITION(rxn,"empty reaction");

    Tags tag;
    uint32_t numReactants,numProducts,numAgents=0;
    
    streamRead(ss,numReactants);
    streamRead(ss,numProducts);
    if(version>1000){
      streamRead(ss,numAgents);
    }
    // we use this here and below to set df_needsInit, so don't re-use the variable
    uint32_t flag = 0;
    streamRead(ss,flag);
    rxn->setImplicitPropertiesFlag(flag & 0x1);

    // -------------------
    //
    // Read Reactants
    //
    // -------------------
    streamRead(ss,tag);
    if(tag != BEGINREACTANTS){
      throw ReactionPicklerException("Bad pickle format: BEGINREACTANTS tag not found.");
    }
    for(unsigned int i=0;i<numReactants;++i){
      ROMol *mol=new ROMol();
      MolPickler::molFromPickle(ss,mol);
      rxn->addReactantTemplate(ROMOL_SPTR(mol));
    }
    streamRead(ss,tag);
    if(tag != ENDREACTANTS){
      throw ReactionPicklerException("Bad pickle format: ENDREACTANTS tag not found.");
    }
    streamRead(ss,tag);
    if(tag != BEGINPRODUCTS){
      throw ReactionPicklerException("Bad pickle format: BEGINPRODUCTS tag not found.");
    }
    for(unsigned int i=0;i<numProducts;++i){
      ROMol *mol=new ROMol();
      MolPickler::molFromPickle(ss,mol);
      rxn->addProductTemplate(ROMOL_SPTR(mol));
    }
    streamRead(ss,tag);
    if(tag != ENDPRODUCTS){
      throw ReactionPicklerException("Bad pickle format: ENDPRODUCTS tag not found.");
    }
    if(numAgents!=0){
      streamRead(ss,tag);
      if(tag != BEGINAGENTS){
        throw ReactionPicklerException("Bad pickle format: BEGINAGENTS tag not found.");
      }
      for(unsigned int i=0;i<numAgents;++i){
        ROMol *mol=new ROMol();
        MolPickler::molFromPickle(ss,mol);
        rxn->addAgentTemplate(ROMOL_SPTR(mol));
      }
      streamRead(ss,tag);
      if(tag != ENDAGENTS){
        throw ReactionPicklerException("Bad pickle format: ENDAGENTS tag not found.");
      }
    }

    // need to do this after we add reactants and products
    rxn->df_needsInit = flag & 0x2;

  } // end of _depickle
}; // end of RDKit namespace
