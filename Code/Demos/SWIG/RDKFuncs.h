// $Id$
//
// Copyright (C) 2008-2010 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <RDGeneral/versions.h>

#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>

#include <Demos/RDKit/Draw/MolDrawing.h>


RDKit::ROMOL_SPTR MolFromSmiles(std::string smi){
  RDKit::ROMol *mol=0;
  try{
    mol=static_cast<RDKit::ROMol *>(RDKit::SmilesToMol(smi));
  } catch (...){
    mol=0;
  }
  return RDKit::ROMOL_SPTR(mol);
};
RDKit::ROMOL_SPTR MolFromSmarts(std::string sma){
  RDKit::ROMol *mol=0;
  try{
    mol=static_cast<RDKit::ROMol *>(RDKit::SmartsToMol(sma));
  } catch (...){
    mol=0;
  }
  return RDKit::ROMOL_SPTR(mol);
};
RDKit::ROMOL_SPTR MolFromMolBlock(std::string molB,
                                  bool sanitize=true,bool removeHs=true){
  RDKit::ROMol *mol=0;
  try{
    mol=static_cast<RDKit::ROMol *>(RDKit::MolBlockToMol(molB,sanitize,removeHs));
  } catch (...){
    mol=0;
  }
  return RDKit::ROMOL_SPTR(mol);
};
RDKit::ROMOL_SPTR MolFromMolFile(std::string filename,
                                 bool sanitize=true,bool removeHs=true){
  RDKit::ROMol *mol=0;
  try{
    mol=static_cast<RDKit::ROMol *>(RDKit::MolFileToMol(filename,sanitize,removeHs));
  } catch (...){
    mol=0;
  }
  return RDKit::ROMOL_SPTR(mol);
};
RDKit::ChemicalReaction *ReactionFromSmarts(std::string sma){
  RDKit::ChemicalReaction *res=0;
  try {
    res=RDKit::RxnSmartsToChemicalReaction(sma);
    if(res) res->initReactantMatchers();
  } catch (...){
    res=0;
  }
  return res;
};
RDKit::ChemicalReaction *ReactionFromRxnBlock(std::string block){
  RDKit::ChemicalReaction *res=0;
  try {
    res=RDKit::RxnBlockToChemicalReaction(block);
    if(res) res->initReactantMatchers();
  } catch (...){
    res=0;
  }
  return res;
};
RDKit::ChemicalReaction *ReactionFromRxnFile(std::string filename){
  RDKit::ChemicalReaction *res=0;
  try {
    res=RDKit::RxnFileToChemicalReaction(filename);
    if(res) res->initReactantMatchers();
  } catch (...){
    res=0;
  }
  return res;
};

std::string MolToSmiles(RDKit::ROMOL_SPTR mol,bool doIsomericSmiles=false,
                        bool doKekule=false, int rootedAtAtom=-1){
  return RDKit::MolToSmiles(*mol,doIsomericSmiles,doKekule,rootedAtAtom);
};
std::string MolToMolBlock(RDKit::ROMOL_SPTR mol, bool includeStereo=true, 
                          int confId=-1) {
  return RDKit::MolToMolBlock(*mol,includeStereo,confId);
}

std::vector<int> MolToBinary(RDKit::ROMOL_SPTR mol){
  std::string sres;
  RDKit::MolPickler::pickleMol(*mol,sres);
  std::vector<int> res(sres.length());
  std::copy(sres.begin(),sres.end(),res.begin());
  return res;
};
RDKit::ROMOL_SPTR MolFromBinary(std::vector<int> pkl){
  std::string sres;
  sres.resize(pkl.size());
  std::copy(pkl.begin(),pkl.end(),sres.begin());
  RDKit::ROMol *res=new RDKit::ROMol(sres);
  return RDKit::ROMOL_SPTR(res);
};


std::vector<int> RxnToBinary(RDKit::ChemicalReaction *rxn){
  std::string sres;
  RDKit::ReactionPickler::pickleReaction(rxn,sres);
  std::vector<int> res(sres.length());
  std::copy(sres.begin(),sres.end(),res.begin());
  return res;
};
RDKit::ChemicalReaction *RxnFromBinary(std::vector<int> pkl){
  std::string sres;
  sres.resize(pkl.size());
  std::copy(pkl.begin(),pkl.end(),sres.begin());
  RDKit::ChemicalReaction *res=new RDKit::ChemicalReaction(sres);
  return res;
};

std::string ReactionToSmarts(RDKit::ChemicalReaction *rxn){
  return RDKit::ChemicalReactionToRxnSmarts(*rxn);
};

std::string rdkitVersion(){
  return RDKit::rdkitVersion;
}

unsigned int compute2DCoords(RDKit::ROMol &mol,bool canonOrient=false,
                             bool clearConfs=true){
  return RDDepict::compute2DCoords(mol,0,canonOrient,clearConfs);
}

unsigned int compute2DCoords(RDKit::ROMol &mol,
                             RDKit::ROMol &templ){
  RDKit::MatchVectType matchVect;
  if(templ.getNumConformers() && SubstructMatch(mol,templ,matchVect)){
    RDGeom::INT_POINT2D_MAP coordMap;
    RDKit::Conformer conf=templ.getConformer();
    for(RDKit::MatchVectType::const_iterator iter=matchVect.begin();
        iter!=matchVect.end();++iter){
      RDGeom::Point2D pt;
      pt.x = conf.getAtomPos(iter->first).x;
      pt.y = conf.getAtomPos(iter->first).y;
      coordMap[iter->second]=pt;
    }
    return RDDepict::compute2DCoords(mol,&coordMap);
  } else {
    return RDDepict::compute2DCoords(mol,0);
  }
}


unsigned int compute3DCoords(RDKit::ROMol &mol,int seed=23,
                             bool clearConfs=true,bool minimize=true){
  int res;
  res= RDKit::DGeomHelpers::EmbedMolecule(mol,0,seed,clearConfs);
  if(res<0){
    res= RDKit::DGeomHelpers::EmbedMolecule(mol,0,seed,clearConfs,true);
  }
  if(res>=0){
    if(minimize){
      ForceFields::ForceField *ff=RDKit::UFF::constructForceField(mol);
      ff->initialize();
      ff->minimize(200);
      delete ff;
    }
  } else {
    BOOST_LOG(rdWarningLog)<<"3D coordination generation failed."<<std::endl;
  }
  return static_cast<unsigned int>(res);
}

std::vector<int> MolToDrawing(const RDKit::ROMol &mol){
  RDKit::RWMol cp(mol);
  RDKit::MolOps::Kekulize(cp);
  if(!mol.getNumConformers()) RDDepict::compute2DCoords(cp);
  std::vector<int> drawing=RDKit::Drawing::DrawMol(cp);
  return drawing;
}

namespace {
  std::string getColor(int atNum){
    static std::map<int,std::string> colors;
    if(colors.empty()){
      colors[7]="#0000FF";
      colors[8]="#FF0000";
      colors[9]="#33CCCC";
      colors[15]="#FF7F00";
      colors[16]="#CCCC00";
      colors[17]="#00CC00";
      colors[35]="#7F4C19";
      colors[0]="#7F7F7F";
    }
    std::string res="#000000";
    if(colors.find(atNum)!=colors.end()) res= colors[atNum];
    return res;
  }
  void drawLine(std::vector<int>::const_iterator &pos,std::ostringstream &sstr){
    int width=*pos++;
    int dashed=*pos++;
    int an1=*pos++;
    int an2=*pos++;
    std::string c1=getColor(an1);
    std::string c2=getColor(an2);
    if(c1==c2){
      sstr<<"<svg:path ";
      sstr<<"d='M "<<*pos<<","<<*(pos+1)<<" "<<*(pos+2)<<","<<*(pos+3)<<"' ";
      pos+=4;
      sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c1<<";stroke-width:4px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1'";
      sstr<<" />\n";
    } else {
      int xp1 = *pos++;
      int yp1 = *pos++;
      int xp2 = *pos++;
      int yp2 = *pos++;
      int mx = xp1+(xp2-xp1)/2;
      int my = yp1+(yp2-yp1)/2;
      sstr<<"<svg:path ";
      sstr<<"d='M "<<xp1<<","<<yp1<<" "<<mx<<","<<my<<"' ";
      sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c1<<";stroke-width:4px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1'";
      sstr<<" />\n";
      sstr<<"<svg:path ";
      sstr<<"d='M "<<mx<<","<<my<<" "<<xp2<<","<<yp2<<"' ";
      sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c2<<";stroke-width:4px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1'";
      sstr<<" />\n";
    }
  }

  void drawAtom(std::vector<int>::const_iterator &pos,std::ostringstream &sstr){
    int fontSz=50;
    int atNum=*pos++;
    int xp=*pos++;
    int yp=*pos++;
    int slen=*pos++;
    std::string label="";
    for(unsigned int i=0;i<slen;++i){
      label+=(char)*pos++;
    }
    int width=fontSz*label.length();
    int height=fontSz;
    sstr<<"<svg:g transform='translate("<<xp<<","<<yp<<")'><svg:rect ";
    sstr<<"style='opacity:1.0;fill:#FFFFFF;stroke:none'";
    sstr<<" width='"<<width<<"' height='"<<height<<"'";
    sstr<<" x='-"<<width/2<<"' y='-"<<height/2<<"'";
    sstr<<"> </svg:rect>\n";
    sstr<<"<svg:text";
    sstr<<" style='font-size:"<<fontSz<<"px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill-opacity:1;stroke:none;font-family:Sans;text-anchor:middle"<<";fill:"<<getColor(atNum)<<"'";
    sstr<<" y='"<<.75*fontSz/2<<"'>";
    sstr<<"<svg:tspan>";
    sstr<<label<<"</svg:tspan>";
    sstr<<"</svg:text>";
    sstr<<"</svg:g>\n";
  }

  std::string ToSVG(const std::vector<int> &drawing){
    std::vector<int>::const_iterator pos=drawing.begin()+2;
    if(*pos!= RDKit::Drawing::BOUNDS){
      std::cerr<<"no bounds token found"<<std::endl;
      return "";
    }
    pos+=3;
    int width=300,height=300;
    width = *pos++;
    height = *pos++;
    std::ostringstream sstr;
    sstr<<"<?xml version='1.0' encoding='iso-8859-1'?>\n";
    sstr << "<svg:svg version='1.1' baseProfile='full'\n      \
        xmlns:svg='http://www.w3.org/2000/svg'\n                \
        xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
        xml:space='preserve'\n";
    sstr<<"width='"<<width<<"px' height='"<<height<<"px' >\n";
    sstr<<"<svg:g transform='translate("<<width*.05<<","<<height*.05<<") scale(.9,.9)'>";
    while(pos!=drawing.end()){
      int token=*pos++;
      switch(token){
      case  RDKit::Drawing::LINE:
        drawLine(pos,sstr);
        break;
      case  RDKit::Drawing::ATOM:
        drawAtom(pos,sstr);
        break;
      default:
        std::cerr<<"unrecognized token: "<<token<<std::endl;
      }
    }

    sstr<<"</svg:g></svg:svg>";
    return sstr.str();
  }
}
std::string MolToSVG(const RDKit::ROMol &mol){
  std::vector<int> drawing=MolToDrawing(mol);
  std::string svg=ToSVG(drawing);
  return svg;
}

