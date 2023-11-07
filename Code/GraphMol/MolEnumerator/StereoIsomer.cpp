#include "MolEnumerator.h"
#include <RDGeneral/Exceptions.h>

#include <boost/range/irange.hpp>
#include <unordered_set>
#include <variant>
#include <GraphMol/MolEnumerator/MolEnumerator.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
namespace MolEnumerator {

namespace {
class StereoFlipper {
 public:
  virtual void flip(bool flag) = 0;
};

class BondStereoFlipper;
class AtomStereoFlipper;
class StereoGroupFlipper;

using stereo_flipper_t =
    std::variant<AtomStereoFlipper, BondStereoFlipper, StereoGroupFlipper>;

class BondStereoFlipper : public StereoFlipper {
 public:
  BondStereoFlipper(Bond* bond) : m_bond(bond) {}

  void flip(bool flag) override;
  static void add_flippers_from_mol(const ROMol& mol,
                                    const StereoEnumerationOptions& options,
                                    std::vector<stereo_flipper_t>& flippers);

 private:
  Bond* m_bond;
};

class AtomStereoFlipper : public StereoFlipper {
 public:
  AtomStereoFlipper(Atom* atom) : m_atom(atom) {}

  void flip(bool flag) override;
  static void add_flippers_from_mol(const ROMol& mol,
                                    const StereoEnumerationOptions& options,
                                    std::vector<stereo_flipper_t>& flippers);

 private:
  Atom* m_atom;
};

class StereoGroupFlipper : public StereoFlipper {
 public:
  StereoGroupFlipper(StereoGroup* group) {
    auto& stereo_group_atoms = group->getAtoms();
    m_input_atom_parities.reserve(stereo_group_atoms.size());

    std::transform(stereo_group_atoms.begin(), stereo_group_atoms.end(),
                   std::back_inserter(m_input_atom_parities), [](auto* atom) {
                     return std::make_pair(atom, atom->getChiralTag());
                   });
  }

  StereoGroupFlipper(StereoGroupFlipper&& other) noexcept = default;

  void flip(bool flag) override;

  static void add_flippers_from_mol(const ROMol& mol,
                                    const StereoEnumerationOptions& options,
                                    std::vector<stereo_flipper_t>& flippers);

 private:
  std::vector<std::pair<Atom*, Atom::ChiralType>> m_input_atom_parities;
};

void BondStereoFlipper::flip(bool flag) {
  auto next_bond_stereo = flag ? Bond::STEREOCIS : Bond::STEREOTRANS;
  m_bond->setStereo(next_bond_stereo);
}

void BondStereoFlipper::add_flippers_from_mol(
    const ROMol& mol, const StereoEnumerationOptions& options,
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
      flippers.push_back(BondStereoFlipper(bond));
    }
  }
}

void AtomStereoFlipper::flip(bool flag) {
  auto next_atom_stereo =
      flag ? Atom::CHI_TETRAHEDRAL_CW : Atom::CHI_TETRAHEDRAL_CCW;
  m_atom->setChiralTag(next_atom_stereo);
}

void AtomStereoFlipper::add_flippers_from_mol(
    const ROMol& mol, const StereoEnumerationOptions& options,
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
      flippers.push_back(AtomStereoFlipper(atom));
    }
  }
}

void StereoGroupFlipper::flip(bool flag) {
  auto flip_atom_parity = [](auto& atom_parity) {
    if (atom_parity != Atom::CHI_TETRAHEDRAL_CW &&
        atom_parity != Atom::CHI_TETRAHEDRAL_CCW) {
      return atom_parity;
    }

    return atom_parity == Atom::CHI_TETRAHEDRAL_CW ? Atom::CHI_TETRAHEDRAL_CCW
                                                   : Atom::CHI_TETRAHEDRAL_CW;
  };

  for (auto& [atom, parity] : m_input_atom_parities) {
    atom->setChiralTag(flag ? parity : flip_atom_parity(parity));
  }
}

void StereoGroupFlipper::add_flippers_from_mol(
    const ROMol& mol, const StereoEnumerationOptions& options,
    std::vector<stereo_flipper_t>& flippers) {
  if (!options.only_unassigned) {
    return;
  }

  for (auto group : mol.getStereoGroups()) {
    if (group.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
      flippers.push_back(StereoGroupFlipper(&group));
    }
  }
}

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
    const ROMol& mol, const StereoEnumerationOptions options, bool verbose) {
  using namespace MolEnumerator;

  MolBundle stereo_isomers;

  std::vector<MolEnumeratorParams> paramsList;
  MolEnumerator::MolEnumeratorParams stereoParams;
  auto sOp = new MolEnumerator::StereoIsomerOp;
  stereoParams.dp_operation =
      std::shared_ptr<MolEnumerator::MolEnumeratorOp>(sOp);
  if (options.max_isomers > 0) {
    stereoParams.maxToEnumerate = options.max_isomers;
  }
  paramsList.push_back(stereoParams);

  // we don't need the values, so we can use temporary values
  std::unordered_set<std::string_view> seen_isomers;
  auto all_stereoisomers = MolEnumerator::enumerate(mol, paramsList);
  for (auto mm : all_stereoisomers.getMols()) {
    if (options.unique) {
      std::string canon_smiles = MolToSmiles(*mm, true);
      if (!seen_isomers.insert(canon_smiles).second) {
        continue;
      }
    }

    if (options.try_embedding) {
      MolOps::addHs(*mm);
      DGeomHelpers::EmbedMolecule(*mm);
    }

    stereo_isomers.addMol(std::move(mm));
  }
  return stereo_isomers;
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
}

std::vector<size_t> StereoIsomerOp::getVariationCounts() const {
  auto flippers = get_flippers(*dp_mol);
  return std::vector<size_t>(flippers.size(), 2);
}

std::unique_ptr<ROMol> StereoIsomerOp::operator()(
    const std::vector<size_t>& which) const {
  PRECONDITION(dp_mol, "no molecule");

  RWMol* mol = new RWMol(*dp_mol);
  auto flippers = get_flippers(*mol);
  if (which.size() != flippers.size()) {
    throw ValueErrorException("bad element choice in enumeration");
  }

  for (size_t i = 0; i < which.size(); ++i) {
    std::visit([&](auto&& flipper) { flipper.flip(which[i]); }, flippers[i]);
  }

  RWMol* isomer;
  if (mol->getStereoGroups().size()) {
    std::vector<StereoGroup> empty_group;
    isomer = new RWMol(*mol);
    isomer->setStereoGroups(std::move(empty_group));

  } else {
    isomer = new RWMol(*mol);
  }
  MolOps::setDoubleBondNeighborDirections(*isomer);
  isomer->clearComputedProps(false);

  MolOps::assignStereochemistry(*isomer, true, true, true);

  return std::unique_ptr<ROMol>(static_cast<ROMol*>(isomer));
}

}  // namespace MolEnumerator

}  // namespace RDKit
