#include "MolEnumerator.h"
#include <RDGeneral/Exceptions.h>

#include <boost/range/irange.hpp>
#include <unordered_set>
#include <GraphMol/MolEnumerator/MolEnumerator.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
namespace MolEnumerator {

inline namespace detail {
class StereoFlipper {
 public:
  virtual void flip(RWMol& mol, bool flag) = 0;
};
}  // namespace detail

namespace {

class AtomStereoFlipper : public StereoFlipper {
 public:
  AtomStereoFlipper(Atom* atom) : d_atom(atom->getIdx()) {}

  void flip(RWMol& mol, bool flag) override {
    auto next_atom_stereo =
        flag ? Atom::CHI_TETRAHEDRAL_CW : Atom::CHI_TETRAHEDRAL_CCW;
    mol.getAtomWithIdx(d_atom)->setChiralTag(next_atom_stereo);
  }
  static void add_flippers_from_mol(const ROMol& mol,
                                    const StereoEnumerationOptions& options,
                                    std::vector<stereo_flipper_t>& flippers) {
    if (options.only_stereo_groups) {
      return;
    }

    for (auto* atom : mol.atoms()) {
      if (!atom->hasProp("_ChiralityPossible")) {
        continue;
      }

      if (!options.only_unassigned ||
          atom->getChiralTag() == Atom::CHI_UNSPECIFIED) {
        flippers.push_back(std::make_shared<AtomStereoFlipper>(atom));
      }
    }
  }

 private:
  int d_atom;
};

class BondStereoFlipper : public StereoFlipper {
 public:
  BondStereoFlipper(Bond* bond) : d_bond(bond->getIdx()) {}

  void flip(RWMol& mol, bool flag) override {
    auto next_bond_stereo = flag ? Bond::STEREOCIS : Bond::STEREOTRANS;
    mol.getBondWithIdx(d_bond)->setStereo(next_bond_stereo);
  }
  static void add_flippers_from_mol(const ROMol& mol,
                                    const StereoEnumerationOptions& options,
                                    std::vector<stereo_flipper_t>& flippers) {
    if (options.only_stereo_groups) {
      return;
    }

    for (auto bond : mol.bonds()) {
      auto bond_stereo = bond->getStereo();
      if (bond_stereo == Bond::STEREONONE) {
        continue;
      }

      if (!options.only_unassigned || bond_stereo == Bond::STEREOANY) {
        flippers.push_back(std::make_shared<BondStereoFlipper>(bond));
      }
    }
  }

 private:
  int d_bond;
};

class StereoGroupFlipper : public StereoFlipper {
 public:
  StereoGroupFlipper(StereoGroup* group) {
    auto& stereo_group_atoms = group->getAtoms();
    d_input_atom_parities.reserve(stereo_group_atoms.size());

    std::transform(stereo_group_atoms.begin(), stereo_group_atoms.end(),
                   std::back_inserter(d_input_atom_parities), [](auto* atom) {
                     return std::make_pair(atom->getIdx(),
                                           atom->getChiralTag());
                   });
  }

  void flip(RWMol& mol, bool flag) override {
    auto flip_atom_parity = [](auto& atom_parity) {
      if (atom_parity != Atom::CHI_TETRAHEDRAL_CW &&
          atom_parity != Atom::CHI_TETRAHEDRAL_CCW) {
        return atom_parity;
      }

      return atom_parity == Atom::CHI_TETRAHEDRAL_CW ? Atom::CHI_TETRAHEDRAL_CCW
                                                     : Atom::CHI_TETRAHEDRAL_CW;
    };

    for (auto& [atom, parity] : d_input_atom_parities) {
      mol.getAtomWithIdx(atom)->setChiralTag(flag ? parity
                                                  : flip_atom_parity(parity));
    }
  }

  static void add_flippers_from_mol(const ROMol& mol,
                                    const StereoEnumerationOptions& options,
                                    std::vector<stereo_flipper_t>& flippers) {
    if (!options.only_unassigned) {
      return;
    }

    for (auto group : mol.getStereoGroups()) {
      if (group.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
        flippers.push_back(std::make_shared<StereoGroupFlipper>(&group));
      }
    }
  }

 private:
  std::vector<std::pair<int, Atom::ChiralType>> d_input_atom_parities;
};

}  // namespace

namespace {

void preprocess_mol_for_stereoisomer_enumeration(RWMol& rwmol) {
  for (auto* atom : rwmol.atoms()) {
    atom->clearProp("_CIPCode");
  }

  for (auto bond : rwmol.bonds()) {
    if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
      bond->setBondDir(Bond::BondDir::NONE);
    }
  }
}

std::vector<stereo_flipper_t> get_flippers(
    const ROMol& input_mol, const StereoEnumerationOptions& options = {}) {
  std::vector<stereo_flipper_t> flippers;

  RWMol mol(input_mol);
  preprocess_mol_for_stereoisomer_enumeration(mol);

  static_cast<void>(RDKit::Chirality::findPotentialStereo(mol));

  AtomStereoFlipper::add_flippers_from_mol(mol, options, flippers);
  BondStereoFlipper::add_flippers_from_mol(mol, options, flippers);
  StereoGroupFlipper::add_flippers_from_mol(mol, options, flippers);

  return flippers;
};

}  // namespace

[[nodiscard]] unsigned int get_stereoisomer_count(
    const ROMol& mol, const StereoEnumerationOptions options) {
  auto flippers = get_flippers(mol, options);
  return std::pow(2, flippers.size());
}

[[nodiscard]] MolBundle enumerate_stereoisomers(
    const ROMol& mol, const StereoEnumerationOptions options) {
  std::vector<MolEnumeratorParams> paramsList;

  MolEnumerator::MolEnumeratorParams stereoParams;
  stereoParams.dp_operation = MolEnumerator::StereoIsomerOp::createOp();
  if (options.max_isomers > 0) {
    stereoParams.maxToEnumerate = options.max_isomers;
  }
  paramsList.push_back(stereoParams);

  // NOTE: we don't need the values, so we don't have to worry about underlying
  // values being invalid
  MolBundle stereoisomers;
  std::unordered_set<std::string_view> seen_isomers;
  auto all_stereoisomers = MolEnumerator::enumerate(mol, paramsList);
  for (auto& stereoisomer : all_stereoisomers.getMols()) {
    if (options.unique) {
      std::string canon_smiles = MolToSmiles(*stereoisomer, true);
      if (!seen_isomers.insert(canon_smiles).second) {
        continue;
      }
    }

    if (options.try_embedding) {
      MolOps::addHs(*stereoisomer);
      DGeomHelpers::EmbedMolecule(*stereoisomer);
    }

    stereoisomers.addMol(std::move(stereoisomer));
  }
  return stereoisomers;
}

StereoIsomerOp::StereoIsomerOp(const std::shared_ptr<ROMol> mol) : dp_mol(mol) {
  PRECONDITION(mol, "bad molecule");
  initFromMol();
}

StereoIsomerOp::StereoIsomerOp(const ROMol& mol) : dp_mol(new ROMol(mol)) {
  initFromMol();
}

std::shared_ptr<StereoIsomerOp> StereoIsomerOp::createOp() {
  return std::shared_ptr<StereoIsomerOp>(new StereoIsomerOp);
}

std::unique_ptr<MolEnumeratorOp> StereoIsomerOp::copy() const {
  return std::unique_ptr<MolEnumeratorOp>(new StereoIsomerOp(*this));
}

StereoIsomerOp::StereoIsomerOp(const StereoIsomerOp& other)
    : dp_mol(other.dp_mol), d_variationPoints(other.d_variationPoints) {}

StereoIsomerOp& StereoIsomerOp::operator=(const StereoIsomerOp& other) {
  if (&other == this) {
    return *this;
  }
  dp_mol = other.dp_mol;
  d_variationPoints = other.d_variationPoints;
  return *this;
}

void StereoIsomerOp::initFromMol(const ROMol& mol) {
  dp_mol.reset(new ROMol(mol));
  initFromMol();
}
void StereoIsomerOp::initFromMol() {
  d_variationPoints.clear();
  if (!dp_mol) {
    return;
  }
  if (!dp_mol->hasProp(detail::idxPropName)) {
    detail::preserveOrigIndices(*dp_mol);
  }

  d_variationPoints = get_flippers(*dp_mol);
}

std::vector<size_t> StereoIsomerOp::getVariationCounts() const {
  return std::vector<size_t>(d_variationPoints.size(), 2);
}

std::unique_ptr<ROMol> StereoIsomerOp::operator()(
    const std::vector<size_t>& which) const {
  PRECONDITION(dp_mol, "no molecule");

  if (which.size() != d_variationPoints.size()) {
    throw ValueErrorException("bad element choice in enumeration");
  }

  auto isomer = std::make_unique<RWMol>(*dp_mol);
  for (size_t i = 0; i < which.size(); ++i) {
    d_variationPoints[i]->flip(*isomer, which[i]);
  }

  if (isomer->getStereoGroups().size()) {
    isomer->setStereoGroups({});  // set to empty
  }

  MolOps::setDoubleBondNeighborDirections(*isomer);
  isomer->clearComputedProps(false);

  MolOps::assignStereochemistry(*isomer, true, true, true);

  return isomer;
}

}  // namespace MolEnumerator

}  // namespace RDKit
