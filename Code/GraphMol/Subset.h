#ifndef RD_SUBSET_H
#define RD_SUBSET_H

#include <RDGeneral/export.h>
#include <memory>
#include <vector>

namespace RDKit {
class RWMol;

//! Subsetting Methods for copyMolSubset
/*
 * These control what the path structures mean in copyMolSubset
 *     \param BONDS_BETWEEN_ATOMS - extract the specified atoms and the bonds between them
 *     \param BONDS - extract the specified bonds and their attached atoms
 *
 *  NOTE: when the atoms and bonds are fully specified, the subset method is ignored.
 */ 

enum class SubsetMethod {
   BONDS_BETWEEN_ATOMS,
     BONDS
};  
  
struct RDKIT_GRAPHMOL_EXPORT SubsetOptions {
  bool sanitize = false;
  bool clearComputedProps = false;
  bool copyAsQuery = false;
  bool copyCoordinates = true;
  unsigned int conformerIdx = -1;
  
  SubsetMethod method = SubsetMethod::BONDS_BETWEEN_ATOMS;

  SubsetOptions() = default;
  
  SubsetOptions(bool sanitize, bool copyAsQuery,
		bool copyCoordinates, unsigned int conformerIdx) :
    sanitize(sanitize), copyAsQuery(copyAsQuery), copyCoordinates(copyCoordinates),
    conformerIdx(conformerIdx)
  {
  }
};

struct RDKIT_GRAPHMOL_EXPORT SubsetInfo {
  std::vector<bool> selected_atoms;
  std::vector<bool> selected_bonds;
  std::map<unsigned int, unsigned int> atom_mapping;
  std::map<unsigned int, unsigned int> bond_mapping;

  void debug(std::ostream &out) {
    out << "Selected Atoms Size: " << selected_atoms.size() << std::endl;
    for(size_t i=0; i<selected_atoms.size(); ++i) {
      if(selected_atoms[i]) out << i << " ";
    }
    out << std::endl;
    out << "Selected Bonds Size: " << selected_bonds.size() << std::endl;
    for(size_t i=0; i<selected_bonds.size(); ++i) {
      if(selected_bonds[i]) out << i << " ";
    }
    out << std::endl;
    out << "Atom Mapping" << std::endl;
    for(auto &v: atom_mapping) {
      std::cerr << " " << v.first << " -> " << v.second << std::endl;
    }
    out << "Bond Mapping" << std::endl;
    for(auto &v: bond_mapping) {
      std::cerr << " " << v.first << " -> " << v.second << std::endl;
    }
  }
};

//!
/*
 * Extract a subgraph from an ROMol. Bonds, atoms, substance groups and
 * stereo groups are only extracted to the subgraph if all participant entities
 * are contained within the given atoms and bonds.
 *
 * \param mol - starting mol
 * \param path - the indices of atoms or bonds to extract. If an index falls
 *             outside of the acceptable indices, it is ignored.yes

 * \param method - the method by which to extract this subgraph.
 * \param sanitize - whether to sanitize the extracted mol.
 *
 * NOTE: Bookmarks are currently copied, StereoGroups that are not entirely
 *        included in the subset are not copied.
 *
 * NOTE: when using this method the SubsetOptions.method is ignored
 *   for resolving the atoms and bonds to use to subset
 */

RDKIT_GRAPHMOL_EXPORT
std::unique_ptr<RDKit::RWMol> copyMolSubset(const RDKit::ROMol& mol,
					    const std::vector<unsigned int> &atoms,
					    const std::vector<unsigned int> &bonds,
					    const SubsetOptions &options, SubsetInfo *selectionInfo=nullptr);

//!
/*
 * Helper api to extract a subgraph from an ROMol. Bonds, substance groups and
 * stereo groups are only extracted to the subgraph if all participant entities
 * are selected by the `path` parameter.
 *
 * \param mol - starting mol
 * \param path - the indices of atoms or bonds to extract. If an index falls
 *             outside of the acceptable indices, it is ignored.yes

 * \param method - the method by which to extract this subgraph.
 * \param sanitize - whether to sanitize the extracted mol.
 *
 * NOTE: Bookmarks are currently copied, StereoGroups that are not entirely
 *        included in the subset are not copied.
 *
 */
RDKIT_GRAPHMOL_EXPORT std::unique_ptr<RDKit::RWMol>
copyMolSubset(const RDKit::ROMol& mol,
	      const std::vector<unsigned int>& path,
	      const SubsetOptions &options = SubsetOptions(),
	      SubsetInfo *mappings = nullptr);


}  // namespace RDKit
#endif
