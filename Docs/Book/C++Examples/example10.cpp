//
// Generating depictions - example10.cpp

#include <iostream>
#include <Geometry/point.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

int main( int argc , char **argv ) {


  auto mol = "c1nccc2n1ccc2"_smiles;
  RDDepict::compute2DCoords( *mol );

 #ifdef RDK_BUILD_COORDGEN_SUPPORT 
  RDDepict::preferCoordGen = true;
#else
  std::cout << "CoordGen support not available" << std::endl;
#endif
  
  RDDepict::compute2DCoords( *mol , nullptr , true );
  
  std::shared_ptr<RDKit::ROMol> templ( RDKit::SmilesToMol( "c1nccc2n1ccc2" ) );
  RDDepict::compute2DCoords( *templ );
  std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "c1cccc2ncn3cccc3c21" ) );
  
  MatchVectType matchVect;
  if(SubstructMatch( *mol1 , *templ , matchVect )) {
    RDKit::Conformer &conf = templ->getConformer();
    RDGeom::INT_POINT2D_MAP coordMap;
    for(auto mv: matchVect) {
      RDGeom::Point3D pt3 = conf.getAtomPos( mv.first );
      RDGeom::Point2D pt2( pt3.x , pt3.y );
      coordMap[mv.second] = pt2;
    }
    RDDepict::compute2DCoords( *mol1 , &coordMap );
  }  

}
