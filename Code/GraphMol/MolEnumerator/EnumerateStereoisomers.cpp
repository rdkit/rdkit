#include <GraphMol/MolEnumerator/EnumerateStereoisomers.h>

#include <GraphMol/MolBundle.h>
#include <GraphMol/ROMol.h>

namespace RDKit {

[[nodiscard]] unsigned int get_stereoisomer_count(
    const ROMol& mol, const StereoEnumerationOptions options) {
  return 0;
}

[[nodiscard]] MolBundle enumerate_stereoisomers(
    const ROMol& mol, const StereoEnumerationOptions options, bool verbose) {
  return {};
}

}  // namespace RDKit
