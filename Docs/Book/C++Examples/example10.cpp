//
// Generating depictions - example10.cpp

#include <iostream>
#include <Geometry/point.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Substruct/SubstructMatch.h>

int main( int argc , char **argv ) {

  RDKit::RWMOL_SPTR mol( new RDKit::RWMol( *RDKit::SmilesToMol( "c1nccc2n1ccc2" ) ) );
  RDDepict::compute2DCoords( *mol );

  RDDepict::compute2DCoords( *mol , static_cast<RDGeom::INT_POINT2D_MAP *>( 0 ) ,
			     true );
  
  RDKit::ROMOL_SPTR templ( RDKit::SmilesToMol( "c1nccc2n1ccc2" ) );
  RDDepict::compute2DCoords( *templ );
  RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "c1cccc2ncn3cccc3c21" ) );
  
  RDKit::MatchVectType matchVect;
  if( RDKit::SubstructMatch( *mol1 , *templ , matchVect ) ) {
    RDKit::Conformer &conf = templ->getConformer();
    RDGeom::INT_POINT2D_MAP coordMap;
    for( RDKit::MatchVectType::const_iterator mv = matchVect.begin() ;
	 mv != matchVect.end() ; ++mv ) {
      RDGeom::Point3D pt3 = conf.getAtomPos( mv->first );
      RDGeom::Point2D pt2( pt3.x , pt3.y );
      coordMap[mv->second] = pt2;
    }
    RDDepict::compute2DCoords( *mol1 , &coordMap );
  }  

}
