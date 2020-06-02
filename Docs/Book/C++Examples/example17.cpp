//
// Depictions with highlights - example17.cpp

#include <fstream>
#include <string>
#include <vector>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>

using namespace RDKit;

std::vector<int> get_all_hit_atoms(ROMol &mol, const std::string &smt) {
  std::vector<int> hit_atoms;
  RWMol *query = SmartsToMol(smt);
  std::vector<MatchVectType> hits_vect;
  SubstructMatch(mol, *query, hits_vect);
  for( size_t i = 0 ; i < hits_vect.size() ; ++i ) {
    for( size_t j = 0 ; j < hits_vect[i].size() ; ++j ) {
      hit_atoms.emplace_back(hits_vect[i][j].second);
    }
  }
  delete query;
  return hit_atoms;
}

std::vector<int> get_all_hit_bonds(ROMol &mol, const std::vector<int> &hit_atoms) {
  std::vector<int> hit_bonds;
  for(int i: hit_atoms) {
    for(int j: hit_atoms) {
      if(i > j) {
	Bond *bnd = mol.getBondBetweenAtoms(i, j);
	if(bnd) {
	  hit_bonds.emplace_back(bnd->getIdx());
	}
      }
    }
  }
  return hit_bonds;
}
    
void update_colour_map(const std::vector<int> &ats, DrawColour col,
		       std::map<int, std::vector<DrawColour> > &ha_map) {
  for(auto h: ats) {
    auto ex = ha_map.find(h);
    if(ex == ha_map.end()) {
      std::vector<DrawColour> cvec(1, col);
      ha_map.insert(make_pair(h, cvec));
    } else {
      if(ex->second.end() == find(ex->second.begin(), ex->second.end(), col)) {
	ex->second.emplace_back(col);
      }
    }
  }
}

int main() {

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book/images/";

  auto mol = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"_smiles;
  std::vector<int> hit1_atoms = get_all_hit_atoms(*mol, "CONN");
  std::vector<int> hit1_bonds = get_all_hit_bonds(*mol, hit1_atoms);
  {
    // all atoms and bonds highlighted the same colour
    std::ofstream outs(file_root + "example_17_1.svg");
    MolDraw2DSVG drawer(500, 500, outs);
    drawer.drawMolecule(*mol, &hit1_atoms, &hit1_bonds);
    drawer.finishDrawing();
    outs.flush();
  }
  {
    // individual colours for each highlight
    std::ofstream outs(file_root + "example_17_2.svg");
    MolDraw2DSVG drawer(500, 500, outs);
    std::map<int, DrawColour> atom_cols;
    atom_cols.insert(std::make_pair(13, DrawColour(1.0, 0.0, 0.0)));
    atom_cols.insert(std::make_pair(14, DrawColour(0.0, 1.0, 0.0)));
    atom_cols.insert(std::make_pair(15, DrawColour(0.0, 0.0, 1.0)));
    atom_cols.insert(std::make_pair(16, DrawColour(1.0, 0.55, 0.0)));
    std::map<int, DrawColour> bond_cols;
    bond_cols.insert(std::make_pair(13, DrawColour(0.8, 0.8, 0.0)));
    bond_cols.insert(std::make_pair(14, DrawColour(0.0, 1.0, 1.0)));
    bond_cols.insert(std::make_pair(15, DrawColour(1.0, 0.0, 1.0)));
    drawer.drawMolecule(*mol, &hit1_atoms, &hit1_bonds, &atom_cols, &bond_cols);
    drawer.finishDrawing();
    outs.flush();
  }

  std::vector<int> hit2_atoms = get_all_hit_atoms(*mol, "N#CC~CO");
  std::vector<int> hit2_bonds = get_all_hit_bonds(*mol, hit2_atoms);
  std::vector<int> hit3_atoms = get_all_hit_atoms(*mol, "C=CON");
  std::vector<int> hit3_bonds = get_all_hit_bonds(*mol, hit3_atoms);
  std::vector<int> hit4_atoms = get_all_hit_atoms(*mol, "CONNCN");
  std::vector<int> hit4_bonds = get_all_hit_bonds(*mol, hit4_atoms);
  
  std::map<int, std::vector<DrawColour> > ha_map;
  update_colour_map(hit1_atoms, DrawColour(1.0, 0.0, 0.0), ha_map);
  update_colour_map(hit2_atoms, DrawColour(0.0, 1.0, 0.0), ha_map);
  update_colour_map(hit3_atoms, DrawColour(0.0, 0.0, 1.0), ha_map);
  update_colour_map(hit4_atoms, DrawColour(1.0, 0.55, 0.0), ha_map);
  std::map<int, std::vector<DrawColour> > hb_map;
  update_colour_map(hit1_bonds, DrawColour(1.0, 0.0, 0.0), hb_map);
  update_colour_map(hit2_bonds, DrawColour(0.0, 1.0, 0.0), hb_map);
  update_colour_map(hit3_bonds, DrawColour(0.0, 0.0, 1.0), hb_map);
  update_colour_map(hit4_bonds, DrawColour(1.0, 0.55, 0.0), hb_map);
  std::map<int, double> h_rads;
  std::map<int, int> h_lw_mult;

  {
    // atoms and bonds can have multiple highlights
    std::ofstream outs(file_root + "example_17_3.svg");
    MolDraw2DSVG drawer(500, 500, outs);
    drawer.drawOptions().fillHighlights = false;
    drawer.drawOptions().continuousHighlight = true;
    drawer.drawOptions().atomHighlightsAreCircles = true;
    drawer.drawMoleculeWithHighlights(*mol, "Test 1", ha_map,
    				      hb_map, h_rads, h_lw_mult);
    drawer.finishDrawing();
    outs.flush();
  }
  
}
