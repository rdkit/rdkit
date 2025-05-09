#ifndef RDKIT_PUBCHEMSHAPE_GUARD
#define RDKIT_PUBCHEMSHAPE_GUARD

#include <GraphMol/ROMol.h>
#include <map>
#include <vector>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

//! The input for the pubchem shape alignment code
struct RDKIT_PUBCHEMSHAPE_EXPORT ShapeInput {
  ShapeInput() = default;
  ShapeInput(const std::string &str) {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    std::stringstream ss(str);
    boost::archive::text_iarchive ia(ss);
    ia &*this;
#endif
  }
  ShapeInput(const ShapeInput &other) = default;
  ShapeInput(ShapeInput &&other) = default;
  ShapeInput &operator=(const ShapeInput &other) = default;
  ShapeInput &operator=(ShapeInput &&other) = default;
  ~ShapeInput() = default;

  std::string toString() const {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa &*this;
    return ss.str();
#endif
  }

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int) {
    ar & coord;
    ar & alpha_vector;
    ar & atom_type_vector;
    ar & volumeAtomIndexVector;
    ar & colorAtomType2IndexVectorMap;
    ar & shift;
    ar & sov;
    ar & sof;
  }
#endif

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

struct RDKIT_PUBCHEMSHAPE_EXPORT ShapeInputOptions {
  bool useColors{true};  // Whether to use colors (pharmacophore features) in
                         // the score.
  bool includeDummies{false};  // Hydrogen and dummy atoms are normally skipped.
                               // This forces the inclusion of dummy atoms.
  double dummyRadius{2.16};    // This is the radius used for Xe.
  std::vector<unsigned int> atomSubset;  // If not empty, use just these atoms
                                         // in the molecule to form the
                                         // ShapeInput object.
  std::vector<unsigned int> notColorAtoms;  // Any atoms mentioned here by index
                                            // should not be used in a color
                                            // feature.
  std::vector<std::pair<unsigned int, double>>
      atomRadii;  // Use these non-standard radii for these atoms.
                  // The int is for the atom index in the molecule, not
                  // the atomic number.
};

//! Prepare the input for the shape comparison
/*!
  \param mol        the molecule to prepare
  \param confId     (optional) the conformer to use
  \param shapeOpts  (optional) Change the default behaviour.

  \return a ShapeInput object, translated to the origin
*/
RDKIT_PUBCHEMSHAPE_EXPORT ShapeInput
PrepareConformer(const RDKit::ROMol &mol, int confId = -1,
                 const ShapeInputOptions &shapeOpts = ShapeInputOptions());

//! Align a shape onto a reference shape.  Assumes the shapes are both
//! centred on the origin.
/*!
  \param refShape      the reference shape
  \param fit           the shape to align
  \param matrix        the transformation matrix (populated on return)
  \param opt_param     (optional) the optimization parameter \param
  max_preiters  (optional) the max number of pre-optimization iterations \param
  max_postiters (optional) the max number of post-optimization iterationsa

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
  if useColors is false)
*/
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> AlignShape(
    const ShapeInput &refShape, ShapeInput &fitShape,
    std::vector<float> &matrix, double opt_param = 1.0,
    unsigned int max_preiters = 10u, unsigned int max_postiters = 30u);

//! Assuming that fitShape has been overlaid onto a reference shape to give
//! the ! transformation matrix, apply the same transformation to the given
//! conformer.
/*!
  \param finalTrans    the final translation to apply to the fitConf coords.
  \param matrix        the transformation matrix produced from alignment
  \param fitShape      the shape that was aligned
  \param fitConf       the conformation to be transformed
*/
RDKIT_PUBCHEMSHAPE_EXPORT void TransformConformer(
    const std::vector<double> &finalTrans, const std::vector<float> &matrix,
    ShapeInput &fitShape, RDKit::Conformer &fitConf);

//! Align a molecule to a reference shape
/*!
  \param refShape      the reference shape
  \param fit           the molecule to align
  \param matrix        the transformation matrix (populated on return)
  \param fitConfId     (optional) the conformer to use for the fit molecule
  \param useColors     (optional) whether or not to use colors in the scoring
  \param opt_param     (optional) the optimization parameter
  \param max_preiters  (optional) the max number of pre-optimization iterations
  \param max_postiters (optional) the max number of post-optimization iterations
  \param applyRefShift (optional) if true, apply the reference shape's shift
                                  translation to the final coordinates

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
  if useColors is false)
*/
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> AlignMolecule(
    const ShapeInput &refShape, RDKit::ROMol &fit, std::vector<float> &matrix,
    int fitConfId = -1, bool useColors = true, double opt_param = 1.0,
    unsigned int max_preiters = 10u, unsigned int max_postiters = 30u,
    bool applyRefShift = false);

//! Align a molecule to a reference molecule
/*!
  \param ref           the reference molecule
  \param fit           the molecule to align
  \param matrix        the transformation matrix (populated on return)
  \param refConfId     (optional) the conformer to use for the reference
                       molecule
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule
  \param useColors     (optional) whether or not to use colors in the
                       scoring
  \param opt_param     (optional) the optimization parameter
  \param max_preiters  (optional) the max number of pre-optimization iterations
  \param max_postiters (optional) the max number of post-optimization
                       iterations

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
  if useColors is false)
*/
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> AlignMolecule(
    const RDKit::ROMol &ref, RDKit::ROMol &fit, std::vector<float> &matrix,
    int refConfId = -1, int fitConfId = -1, bool useColors = true,
    double opt_param = 1.0, unsigned int max_preiters = 10u,
    unsigned int max_postiters = 30u);

#endif
