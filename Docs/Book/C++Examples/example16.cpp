//
// Substructure searching with SMARTS atom map indices - example16.cpp

#include <iostream>
#include <map>
#include <vector>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

int main( int argc , char **argv ) {

  std::shared_ptr<RDKit::RWMol> patt1( RDKit::SmartsToMol( "[cH0:1][c:2]([cH0])!@[CX3!r:3]=[NX2!r:4]" ) );
  std::map<int,unsigned int> ind_map;
  for(auto atom: patt1->atoms()) {
    int map_num = atom->getAtomMapNum();
    if( map_num ) {
      ind_map[map_num-1] = atom->getIdx();
    }
  }
  std::vector<unsigned int> map_list;
  for(auto im: ind_map) {
    map_list.push_back(im.second);
  }
  for( size_t i = 0 , is = map_list.size() ; i < is ; ++i ) {
    std::cout << map_list[i] << " ";
  }
  std::cout << std::endl;

  std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "Cc1cccc(C)c1C(C)=NC" ) );
  std::vector<RDKit::MatchVectType> hits_vect;
  if( RDKit::SubstructMatch( *mol1 , *patt1 , hits_vect ) ) {
    for( size_t i = 0 ; i < hits_vect.size() ; ++i ) {
      std::cout << "Match " << i + 1 << " : ";
      for( size_t j = 0 ; j < map_list.size() ; ++j ) {
	std::cout << hits_vect[i][map_list[j]].second << " ";
      }
      std::cout << std::endl;
    }
  }
  
}
