//  Copyright (C) 2020 Eisuke Kawashima
//
//  Chemical Markup Language, CML, schema 3
//  http://www.xml-cml.org/
//  See
//  http://www.xml-cml.org/convention/molecular
//  http://www.xml-cml.org/schema/schema3/schema.xsd
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <RDGeneral/Invariant.h>

#if __cplusplus >= 201402L && __has_cpp_attribute(fallthrough)
#define FALLTHROUGH [[fallthrough]]
#else
#define FALLTHROUGH
#endif

namespace RDKit {
namespace {
boost::property_tree::ptree molToPTree(const ROMol& mol, int confId,
                                       bool kekulize) {
  RWMol rwmol{mol};
  if (kekulize) {
    MolOps::Kekulize(rwmol);
  }

  boost::property_tree::ptree pt;
  // prefix namespaces
  auto& root = pt.add("cml", "");
  root.put("<xmlattr>.xmlns", "http://www.xml-cml.org/schema");
  root.put("<xmlattr>.xmlns:convention", "http://www.xml-cml.org/convention/");
  root.put("<xmlattr>.convention", "convention:molecular");

  auto& molecule = root.add("molecule", "");

  // molecule/@id MUST start with an alphabetical character
  // http://www.xml-cml.org/convention/molecular#molecule-id
  const auto molecule_id_default_prefix = "m";
  molecule.put("<xmlattr>.id",
               boost::format{"%1%%2%"} % molecule_id_default_prefix % confId);

  std::string name;
  rwmol.getPropIfPresent(common_properties::_Name, name);
  if (!name.empty()) {
    molecule.put("name", name);
  }

  int mol_formal_charge = 0;
  unsigned mol_num_radical_electrons = 0u;

  // atom/@id MUST start with an alphabetical character
  // http://www.xml-cml.org/convention/molecular#atom-id
  const auto atom_id_prefix = "a";

  const Conformer* conf = nullptr;
  if (rwmol.getNumConformers()) {
    conf = &rwmol.getConformer(confId);
    // wedge bonds so that we can use that info:
    WedgeMolBonds(rwmol, conf);
  }
  auto& atomArray = molecule.put("atomArray", "");
  for (unsigned i = 0u, nAtoms = rwmol.getNumAtoms(); i < nAtoms; i++) {
    auto& atom = atomArray.add("atom", "");
    const auto& a = rwmol.getAtomWithIdx(i);

    atom.put("<xmlattr>.id", boost::format{"%1%%2%"} % atom_id_prefix % i);
    if (a->getAtomicNum()) {
      atom.put("<xmlattr>.elementType", a->getSymbol());
    } else {
      atom.put("<xmlattr>.elementType", "Du");  // dummy
    }

    const auto charge = a->getFormalCharge();
    mol_formal_charge += charge;
    atom.put("<xmlattr>.formalCharge", charge);

    atom.put("<xmlattr>.hydrogenCount", a->getTotalNumHs(true));

    const auto isotope = a->getIsotope();
    if (isotope) {
      atom.put("<xmlattr>.isotopeNumber", isotope);
    }

    const auto n_rad_es = a->getNumRadicalElectrons();
    mol_num_radical_electrons += n_rad_es;

    if (conf != nullptr) {
      const auto& pos = conf->getAtomPos(i);
      boost::format xyz_fmt{"%.6f"};

      if (!conf->is3D()) {
        atom.put("<xmlattr>.x2", xyz_fmt % pos.x);
        atom.put("<xmlattr>.y2", xyz_fmt % pos.y);
      } else {
        atom.put("<xmlattr>.x3", xyz_fmt % pos.x);
        atom.put("<xmlattr>.y3", xyz_fmt % pos.y);
        atom.put("<xmlattr>.z3", xyz_fmt % pos.z);
      }
    }
    // atom/@atomParity if chiral
    // http://www.xml-cml.org/convention/molecular#atom-atomParity
    // the parity is the sign of the chiral volume. We can determine that from
    // the ChiralTag:
    int parity = 0;
    switch (a->getChiralTag()) {
      case Atom::CHI_TETRAHEDRAL_CCW:
        parity = 1;
        break;
      case Atom::CHI_TETRAHEDRAL_CW:
        parity = -1;
        break;
      default:
        parity = 0;
    }
    if (parity) {
      std::vector<unsigned> neighbors;
      for (auto nbri : boost::make_iterator_range(rwmol.getAtomNeighbors(a))) {
        const auto at = rwmol[nbri];
        neighbors.push_back(at->getIdx());
      }

      auto& atomParity = atom.add("atomParity", parity);
      atomParity.put("<xmlattr>.atomRefs4",
                     boost::format{"%1%%2% %1%%3% %1%%4% %1%%5%"} %
                         atom_id_prefix % neighbors[0u] % neighbors[1u] %
                         neighbors[2u] % neighbors[3u]);
    }
  }

  molecule.put("<xmlattr>.formalCharge", mol_formal_charge);
  if (mol_num_radical_electrons < 2u) {
    molecule.put("<xmlattr>.spinMultiplicity", mol_num_radical_electrons + 1u);
  } else {
    BOOST_LOG(rdInfoLog)
        << "CMLWriter: Unable to determine molecule/@spinMultiplicity "
        << boost::format{"(%1% radical electrons)\n"} %
               mol_num_radical_electrons;
  }

  // bond/@id so that it can be referenced
  // http://www.xml-cml.org/convention/molecular#bond-id
  const auto bond_id_prefix = "b";
  unsigned bond_id = 0u;

  auto& bondArray = molecule.add("bondArray", "");
  for (auto atom_itr = rwmol.beginAtoms(), atom_itr_end = rwmol.endAtoms();
       atom_itr != atom_itr_end; ++atom_itr) {
    const auto& atom = *atom_itr;
    PRECONDITION(atom, "bad atom");
    const auto src = atom->getIdx();
    for (auto bond_itrs = rwmol.getAtomBonds(atom);
         bond_itrs.first != bond_itrs.second; ++bond_itrs.first) {
      auto* bptr = rwmol[*bond_itrs.first];
      auto* nptr = bptr->getOtherAtom(atom);
      const auto dst = nptr->getIdx();
      if (dst < src) {
        continue;
      }

      auto& bond = bondArray.add("bond", "");
      bond.put("<xmlattr>.atomRefs2",
               boost::format{"%1%%2% %1%%3%"} % atom_id_prefix % src % dst);

      bond.put("<xmlattr>.id",
               boost::format{"%1%%2%"} % bond_id_prefix % bond_id++);

      const auto btype = bptr->getBondType();
      switch (btype) {
        case Bond::SINGLE:
          bond.put("<xmlattr>.order", 'S');
          break;
        case Bond::DOUBLE:
          bond.put("<xmlattr>.order", 'D');
          break;
        case Bond::TRIPLE:
          bond.put("<xmlattr>.order", 'T');
          break;
        case Bond::AROMATIC:
          bond.put("<xmlattr>.order", 'A');
          break;
        case Bond::DATIVEONE:
          FALLTHROUGH;
        case Bond::DATIVE:
          FALLTHROUGH;
        case Bond::DATIVEL:
          FALLTHROUGH;
        case Bond::DATIVER:
          bond.put("<xmlattr>.order", 'S');
          break;
        // XXX RDKit extension: bond orders greater than 3
        case Bond::QUADRUPLE:
          bond.put("<xmlattr>.order", 4);
          break;
        case Bond::QUINTUPLE:
          bond.put("<xmlattr>.order", 5);
          break;
        case Bond::HEXTUPLE:
          bond.put("<xmlattr>.order", 6);
          break;
        // XXX RDKit extension: half-integer orders
        case Bond::THREECENTER:
          bond.put("<xmlattr>.order", "0.5");
          break;
        case Bond::ONEANDAHALF:
          bond.put("<xmlattr>.order", "1.5");
          break;
        case Bond::TWOANDAHALF:
          bond.put("<xmlattr>.order", "2.5");
          break;
        case Bond::THREEANDAHALF:
          bond.put("<xmlattr>.order", "3.5");
          break;
        case Bond::FOURANDAHALF:
          bond.put("<xmlattr>.order", "4.5");
          break;
        case Bond::FIVEANDAHALF:
          bond.put("<xmlattr>.order", "5.5");
          break;
        default:
          BOOST_LOG(rdInfoLog)
              << boost::format{"CMLWriter: Unsupported BondType %1%\n"} % btype;
          bond.put("<xmlattr>.order", "");
      }
      // bond/@BondStereo if appropriate
      // http://www.xml-cml.org/convention/molecular#bondStereo-element
      auto bdir = bptr->getBondDir();
      switch (bdir) {
        case Bond::BondDir::BEGINDASH:
          bond.put("<xmlattr>.bondStereo", "H");
          break;
        case Bond::BondDir::BEGINWEDGE:
          bond.put("<xmlattr>.bondStereo", "W");
          break;
        default:
          break;
      }
    }
  }

  return pt;
}
}  // namespace

void MolToCMLBlock(std::ostream& os, const ROMol& mol, int confId,
                   bool kekulize) {
  auto pt = molToPTree(mol, confId, kekulize);
  if (pt.empty()) {
    return;
  }
  boost::property_tree::write_xml(
      os, pt,
      boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
}

std::string MolToCMLBlock(const ROMol& mol, int confId, bool kekulize) {
  std::ostringstream ss;
  MolToCMLBlock(ss, mol, confId, kekulize);
  return ss.str();
}

void MolToCMLFile(const ROMol& mol, const std::string& fName, int confId,
                  bool kekulize) {
  std::ofstream ofs{fName};
  MolToCMLBlock(ofs, mol, confId, kekulize);
}

}  // namespace RDKit
