//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/types.h>
#include "RDFreeSASA.h"
#include "boost/format.hpp"

extern "C" {
#include "freesasa.h"
}

namespace RDKit {
  namespace common_properties {
    namespace Atom {
      const std::string SASA = "SASA";; // Solvent Accessible Surface Area for atom- double
      const std::string SASAClass = "SASAClass";     // Class type, 0,1,2... etc
      const std::string SASAClassName = "SASAClassName"; // Class name, Polar, APolar etc...
    }
    namespace Molecule {
      const std::string SASA = "SASA"; // Total Solvent Accessible Surface area for molecule;
    }
  }
}

namespace FreeSASA {
using namespace RDKit;

bool classifyAtoms(ROMol &mol, std::vector<double> &radii,
                   const SASAOpts &opts) {
  radii.clear();
  const freesasa_classifier *classifier = 0;
  switch (opts.d_classifier) {
    case SASAOpts::Protor:
      classifier = &freesasa_protor_classifier;
      break;
    case SASAOpts::NACCESS:
      classifier = &freesasa_naccess_classifier;
      break;
    case SASAOpts::OONS:
      classifier = &freesasa_oons_classifier;
      break;
    default:
      throw ValueErrorException("unknown FreeSASA classifier specified");
      return false;
  }

  bool success = true;
  for (ROMol::AtomIterator at = mol.beginAtoms(); at != mol.endAtoms(); ++at) {
    Atom *atom = *at;
    int cls = -1;
    std::string classification = "Unclassified";
    double radius = 0.0;

    const AtomMonomerInfo *info = atom->getMonomerInfo();
    if (info) {
      const char *atom_name = info->getName().c_str();
      const char *res_name = 0;

      if (info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
        res_name = ((AtomPDBResidueInfo *)info)->getResidueName().c_str();
        radius = classifier->radius(res_name, atom_name, classifier);

        if (radius == 0.0) {
          BOOST_LOG(rdWarningLog) << "Atom " << atom->getIdx()
                                  << " has zero radius" << std::endl;
        }

        cls = classifier->sasa_class(res_name, atom_name, classifier);
        if (cls == FREESASA_FAIL) {
          std::string err = boost::str(
              boost::format("Atom %d% failed classification") % atom->getIdx());
          throw ValueErrorException(err);
        } else if (cls == FREESASA_WARN) {
          BOOST_LOG(rdWarningLog) << "Atom " << atom->getIdx()
                                  << " could not be classified" << std::endl;
          success = false;
        } else {
          classification = classifier->class2str(cls, classifier);
        }
      }
    }

    radii.push_back(radius);
    atom->setProp(common_properties::Atom::SASAClass, cls);
    atom->setProp(common_properties::Atom::SASAClassName, classification);
  }

  return success;
}

double calcSASA(const ROMol &mol, const std::vector<double> &radii,
                const SASAOpts &opts) {
  return calcSASA(mol, radii, -1, opts);
}

double calcSASA(const ROMol &mol, const std::vector<double> &radii, int confIdx,
                const SASAOpts &opts) {
  PRECONDITION(mol.getNumConformers(), "No conformers in molecule");
  PRECONDITION(confIdx < rdcast<int>(mol.getNumConformers()),
               "Conformer index out of range");
  PRECONDITION(mol.getNumAtoms(), "Empty molecule");

  freesasa_parameters params = freesasa_default_parameters;
  params.n_threads = 1;
  switch (opts.d_alg) {
    case SASAOpts::LeeRichards:
      params.alg = FREESASA_LEE_RICHARDS;
      break;
    case SASAOpts::ShrakeRupley:
      params.alg = FREESASA_SHRAKE_RUPLEY;
      break;
    default:
      throw ValueErrorException("Unknown freesasa algorithm");
  }

  // sneaky, but legal :)
  std::vector<double> coords(mol.getNumAtoms() * 3);
  const RDGeom::POINT3D_VECT &vec = mol.getConformer(confIdx).getPositions();
  for (size_t i = 0; i < mol.getNumAtoms(); ++i) {
    coords[i * 3] = vec[i].x;
    coords[i * 3 + 1] = vec[i].y;
    coords[i * 3 + 2] = vec[i].z;
  }

  freesasa_result *res =
      freesasa_calc_coord(&coords[0], &radii[0], mol.getNumAtoms(), &params);
  if (!res) return 0.0;
  CHECK_INVARIANT(res->n_atoms == rdcast<int>(mol.getNumAtoms()),
                  "freesasa didn't return the correct number of atoms");

  double sasa = res->total;
  mol.setProp(common_properties::Molecule::SASA, sasa);
  size_t i = 0;
  for (ROMol::ConstAtomIterator at = mol.beginAtoms(); at != mol.endAtoms();
       ++at, ++i) {
    (*at)->setProp(common_properties::Atom::SASA, res->sasa[i]);
  }

  freesasa_result_free(res);
  return sasa;
}

double calcSASA(const RDKit::ROMol &mol, const std::vector<double> &radii,
                const RDKit::QueryAtom *query,
                const SASAOpts &opts) {
  return calcSASA(mol, radii, query, -1, opts);
}

double calcSASA(const RDKit::ROMol &mol, const std::vector<double> &radii,
                const RDKit::QueryAtom *query,
                int confIdx,
                const SASAOpts &opts) {
  calcSASA(mol, radii, confIdx, opts);
  double total = 0.0f;
  for (ROMol::ConstQueryAtomIterator at = mol.beginQueryAtoms(query);
       at != mol.endQueryAtoms(); ++at) {
    const Atom *atom = *at;
    total += atom->getProp<double>("SASA");
  }
  return total;
}

const RDKit::QueryAtom * makeFreeSasaAPolarAtomQuery() {
  QueryAtom *qa = new QueryAtom;
  qa->setQuery(makePropQuery<Atom, std::string>("SASAClassName", "Apolar"));
  return qa;
}

const RDKit::QueryAtom * makeFreeSasaPolarAtomQuery() {
  QueryAtom *qa = new QueryAtom;
  qa->setQuery(makePropQuery<Atom, std::string>("SASAClassName", "Polar"));
  return qa;
}
}



