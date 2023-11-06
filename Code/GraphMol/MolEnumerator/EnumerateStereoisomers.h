
#include <RDGeneral/export.h>

namespace RDKit {

class MolBundle;
class ROMol;

struct RDKIT_MOLENUMERATOR_EXPORT StereoEnumerationOptions {
  bool try_embedding;
  bool only_unassigned;
  bool only_stereo_groups;
  bool unique;
  unsigned int max_isomers;
  unsigned int rand;
};

[[nodiscard]] RDKIT_MOLENUMERATOR_EXPORT unsigned int get_stereoisomer_count(
    const ROMol& mol, const StereoEnumerationOptions options = {});

[[nodiscard]] RDKIT_MOLENUMERATOR_EXPORT MolBundle enumerate_stereoisomers(
    const ROMol& mol, const StereoEnumerationOptions options = {},
    bool verbose = false);

}  // namespace RDKit
