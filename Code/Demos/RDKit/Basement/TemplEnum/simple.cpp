//  $Id$
//
//  Copyright (C) 2004 Rational Discovery LLC
//   All Rights Reserved
//
#include "TemplEnum.h"
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;
using namespace TemplateEnum;

#include <math.h>

bool feq(double v1,double v2,double tol=1e-4){
  return fabs(v1-v2)<=tol;
}
bool feq(RDGeom::Point3D p1,RDGeom::Point3D p2,double tol=1e-4){
  return feq(p1.x,p2.x,tol)&&feq(p1.y,p2.y,tol)&&feq(p1.z,p2.z,tol);
}


int main(int argc,const char *argv[]){
  if(argc<3){
    std::cerr << " Usage: simple.exe <template.mol> [attach1.sdf, attach2.sdf, ...]" << std::endl;
    exit(-1);
  }
  std::vector<const char *>fileNames;
  for(int i=2;i<argc;i++){
    fileNames.push_back( argv[i] );
  }
  RWMOL_SPTR_VECT library=enumFromFiles(argv[1],fileNames);

  std::cerr << "Created: " << library.size() << " compounds." << std::endl;
  for(int i=0;i<library.size();i++){
    RWMOL_SPTR mol = library[i];
    std::cout << MolToMolBlock(mol.get(),false);
    std::cout << "$$$$" << std::endl;
  }
  
  exit(0);

}
