#include <GraphMol/ROMol.h>
#include <map>
#include <vector>

using ShapeInput = struct ShapeInputStruct {
  std::vector<float> coord;
  std::vector<double> alpha_vector;
  std::vector<unsigned int> atom_type_vector;
  std::vector<unsigned int> volumeAtomIndexVector;
  std::map<unsigned int, std::vector<unsigned int>>
      colorAtomType2IndexVectorMap;
  std::vector<double> shift;
  double sov{0.0};
  double sof{0.0};
};

ShapeInput PrepareConformer(const RDKit::ROMol &mol, int confId = -1,
                            bool useColors = true);

std::pair<double, double>
AlignMolecule(const ShapeInput &refShape, RDKit::ROMol &fit,
              std::vector<float> &matrix, int fitConfId = -1,
              bool useColors = true, double opt_param = 0.5,
              unsigned int max_preiters = 3u, unsigned int max_postiters = 16u);

std::pair<double, double>
AlignMolecule(const RDKit::ROMol &ref, RDKit::ROMol &fit,
              std::vector<float> &matrix, int refConfId = -1,
              int fitConfId = -1, bool useColors = true, double opt_param = 0.5,
              unsigned int max_preiters = 3u, unsigned int max_postiters = 16u);
