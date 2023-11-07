#include <GraphMol/MolEnumerator/EnumerateStereoisomers.h>

#include <GraphMol/MolBundle.h>
#include <GraphMol/ROMol.h>

namespace RDKit {

[[nodiscard]] MolBundle enumerate_stereoisomers(
    const ROMol& mol, const StereoEnumerationOptions options, bool verbose) {
  return {};
}

}  // namespace RDKit
