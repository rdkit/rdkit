#include "MolEnumerator.h"
#include <RDGeneral/Exceptions.h>

#include <variant>
#include <GraphMol/MolEnumerator/EnumerateStereoisomers.h>
#include <GraphMol/Chirality.h>

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

}  // namespace

std::vector<stereo_flipper_t> get_flippers(
    const ROMol& mol, const StereoEnumerationOptions& options) {
  auto potential_stereo = RDKit::Chirality::findPotentialStereo(mol);

  std::vector<stereo_flipper_t> flippers;

  AtomStereoFlipper::add_flippers_from_mol(mol, options, flippers);
  BondStereoFlipper::add_flippers_from_mol(mol, options, flippers);
  StereoGroupFlipper::add_flippers_from_mol(mol, options, flippers);

  return flippers;
};

unsigned int get_stereoisomer_count(
    ROMol* mol,
    const StereoEnumerationOptions options = StereoEnumerationOptions()) {
  auto flippers = get_flippers(*mol, options);
  return std::pow(2, flippers.size());
}

/*
std::vector<ROMol*> enumerate_stereoisomers(
    const ROMol& input_mol,
    const StereoEnumerationOptions options = StereoEnumerationOptions(),
    bool verbose = false) {

  RWMol rwmol(input_mol);
  for (auto* atom : rwmol.atoms()) {
    atom->clearProp("_CIPCode");
  }

  for (auto bond : rwmol.bonds()) {
    if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
      bond->setBondDir(Bond::BondDir::NONE);
    }
  }
  auto flippers = get_flippers(mol, options);
  const unsigned int n_centers = flippers.size();

  if (!n_centers) {
    return std::vector<ROMol*>{mol};
  }

  if (options.max_isomers == 0 ||
      std::pow(2, n_centers) <= options.max_isomers) {
    std::vector<unsigned int> bitsource = _RangeBitsGenerator(n_centers);
  } else {
    if (!options.rand) {
      // deterministic random seed invariant to input atom order
      std::vector<std::pair<int, int>> ordered_atoms;
      for (auto atom : mol->atoms()) {
        ordered_atoms.push_back(
            std::make_pair(atom->getDegree(), atom->getAtomicNum()));
      }
      std::sort(ordered_atoms.begin(), ordered_atoms.end());
      std::size_t seed =
          std::hash<std::vector<std::pair<int, int>>>(ordered_atoms);
      rand = srand(seed);
    } else {
      rand = srand(options.rand);
    }
  }

  bitsource = _UniqueRandomBitsGenerator(n_centers, options.max_isomers, rand);

  std::set<std::string> seen_isomers;
  int num_isomers = 0;
  for (auto bitflag : bitsource) {
    for (int i = 0; i < n_centers; ++i) {
      flippers[i]->flip(bitflag & (1 << i));
    }

    RWMol* isomer;
    if (mol->getStereoGroups()) {
      std::vector<StereoGroup> empty_group;
      isomer = new RWMol(*mol);
      isomer->setStereoGroups(std::move(empty_group));

    } else {
      isomer = new RWMol(*mol);
    }
    MolOps::setDoubleBondNeighborDirections(*isomer);
    isomer->clearComputedProps(false);

    MolOps::assignStereochemistry(*isomer, true, true, true);
    if (options.unique) {
      std::string cansmi = MolToSmiles(*isomer, true);
      if (seen_isomers.find(cansmi) != seen_isomers.end()) {
        continue;
      }

      seen_isomers.insert(cansmi);
    }

    if (options.try_embedding) {
      MolOps::addHs(*isomer);
      DGeomHelpers::EmbedMolecule(*isomer, bitflag & 0x7fffffff);
    }
  }
}
*/

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
  std::cout << "hereeeeeeeeeeee" << std::endl;
  /*
  for (const auto bond : dp_mol->bonds()) {
    std::string endpts;
    std::string attach;
    if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, endpts) &&
        bond->getPropIfPresent(common_properties::_MolFileBondAttach, attach) &&
        attach == "ANY") {
      const Atom *atom = bond->getBeginAtom();
      if (atom->getAtomicNum() == 0) {
        atom = bond->getEndAtom();
        if (atom->getAtomicNum() == 0) {
          // marvin sketch seems to place the position-variation dummy at the
          // beginning of the bond, so we're going to favor taking the end atom.
          // In case other tools construct this differently, we have an
          // exception to that if the end atom is an AtomNull query and the
          // beginning atom is not one.
          if (atom->hasQuery() &&
              atom->getQuery()->getDescription() == "AtomNull" &&
              bond->getBeginAtom()->hasQuery() &&
              bond->getBeginAtom()->getQuery()->getDescription() !=
                  "AtomNull") {
            atom = bond->getBeginAtom();
          }
        }
      }
      d_dummiesAtEachPoint.push_back(bond->getOtherAtomIdx(atom->getIdx()));
      std::vector<unsigned int> oats =
          RDKit::SGroupParsing::ParseV3000Array<unsigned int>(
              endpts, dp_mol->getNumAtoms(), false);
      // decrement the indices and do error checking and whatever additional
      // cleanup is required:
      for (auto &oat : oats) {
        if (oat == 0 || oat > dp_mol->getNumAtoms()) {
          throw ValueErrorException("Bad variation point index");
        }
        --oat;
        // github #4381: if we're connecting to an aromatic heteroatom which
        // has implicit Hs, we should remove those
        auto attachAtom = dp_mol->getAtomWithIdx(oat);
        if (attachAtom->getIsAromatic() && attachAtom->getAtomicNum() != 6) {
          attachAtom->setNumExplicitHs(0);
        }
      }
      d_variationPoints.push_back(std::make_pair(atom->getIdx(), oats));
    }
  }
  */
}

std::vector<size_t> StereoIsomerOp::getVariationCounts() const {
  std::vector<size_t> res(d_variationPoints.size());
  std::transform(
      d_variationPoints.begin(), d_variationPoints.end(), res.begin(),
      [](std::pair<unsigned int, std::vector<unsigned int>> pr) -> size_t {
        return pr.second.size();
      });
  return res;
}

std::unique_ptr<ROMol> StereoIsomerOp::operator()(
    const std::vector<size_t>& which) const {
  PRECONDITION(dp_mol, "no molecule");
  if (which.size() != d_variationPoints.size()) {
    throw ValueErrorException("bad element choice in enumeration");
  }
  // a bit of quick error checking before starting the real work
  for (unsigned int i = 0; i < d_variationPoints.size(); ++i) {
    if (which[i] >= d_variationPoints[i].second.size()) {
      throw ValueErrorException("bad element value in enumeration");
    }
  }
  RWMol* res = new RWMol(*dp_mol);
  for (unsigned int i = 0; i < d_variationPoints.size(); ++i) {
    const auto tpl = d_variationPoints[i];
    auto begAtomIdx = tpl.first;
    auto endAtomIdx = tpl.second[which[i]];
    // do we already have a bond?
    if (res->getBondBetweenAtoms(begAtomIdx, endAtomIdx)) {
      continue;
    }
    res->addBond(begAtomIdx, endAtomIdx, Bond::BondType::SINGLE);
  }
  // now remove the dummies:
  res->beginBatchEdit();
  for (auto idx : d_dummiesAtEachPoint) {
    res->removeAtom(idx);
  }
  res->commitBatchEdit();
  return std::unique_ptr<ROMol>(static_cast<ROMol*>(res));
}

}  // namespace MolEnumerator

}  // namespace RDKit
