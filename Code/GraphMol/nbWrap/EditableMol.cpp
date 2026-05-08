//
//  Copyright (C) 2007-2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>
#include <nanobind/nanobind.h>

// ours
#include <GraphMol/RDKitBase.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace {
class EditableMol {
 public:
  EditableMol(const ROMol &m) { dp_mol = new RWMol(m); };
  ~EditableMol() noexcept { delete dp_mol; };

  void RemoveAtom(unsigned int idx) {
    PRECONDITION(dp_mol, "no molecule");
    dp_mol->removeAtom(idx);
  };
  void RemoveBond(unsigned int idx1, unsigned int idx2) {
    PRECONDITION(dp_mol, "no molecule");
    dp_mol->removeBond(idx1, idx2);
  };
  int AddBond(unsigned int begAtomIdx, unsigned int endAtomIdx,
              Bond::BondType order = Bond::UNSPECIFIED) {
    PRECONDITION(dp_mol, "no molecule");
    return dp_mol->addBond(begAtomIdx, endAtomIdx, order);
  };
  int AddAtom(Atom *atom) {
    PRECONDITION(dp_mol, "no molecule");
    PRECONDITION(atom, "bad atom");
    return dp_mol->addAtom(atom, true, false);
  };
  void ReplaceAtom(unsigned int idx, Atom *atom, bool updateLabels,
                   bool preserveProps) {
    PRECONDITION(dp_mol, "no molecule");
    PRECONDITION(atom, "bad atom");
    dp_mol->replaceAtom(idx, atom, updateLabels, preserveProps);
  };
  void ReplaceBond(unsigned int idx, Bond *bond, bool preserveProps) {
    PRECONDITION(dp_mol, "no molecule");
    PRECONDITION(bond, "bad bond");
    dp_mol->replaceBond(idx, bond, preserveProps);
  };
  void BeginBatchEdit() {
    PRECONDITION(dp_mol, "no molecule");
    dp_mol->beginBatchEdit();
  };
  void RollbackBatchEdit() {
    PRECONDITION(dp_mol, "no molecule");
    dp_mol->rollbackBatchEdit();
  };
  void CommitBatchEdit() {
    PRECONDITION(dp_mol, "no molecule");
    dp_mol->commitBatchEdit();
  };

  ROMol *GetMol() const {
    PRECONDITION(dp_mol, "no molecule");
    auto *res = new ROMol(*dp_mol);
    return res;
  };

 private:
  RWMol *dp_mol;
};
}  // namespace

struct EditableMol_wrapper {
  static void wrap(nb::module_ &m) {
    std::string molClassDoc =
        R"DOC(The EditableMol class.

   This class can be used to add/remove bonds and atoms to
   a molecule.
   In order to use it, you need to first construct an EditableMol
   from a standard Mol:

   >>> m = Chem.MolFromSmiles('CCC')
   >>> em = Chem.EditableMol(m)
   >>> em.AddAtom(Chem.Atom(8))
   >>> em.AddBond(0,3,Chem.BondType.SINGLE)
   >>> m2 = em.GetMol()
   >>> Chem.SanitizeMol(m2)
   >>> Chem.MolToSmiles(m2)
   'CCCO'

   *Note*: It is very, very easy to shoot yourself in the foot with
     this class by constructing an unreasonable molecule.
)DOC";
    nb::class_<EditableMol>(m, "EditableMol", molClassDoc.c_str())
        .def(nb::init<const ROMol &>(), "m"_a, R"DOC(Construct from a Mol)DOC")
        .def("RemoveAtom", &EditableMol::RemoveAtom, "idx"_a,
             R"DOC(Remove the specified atom from the molecule)DOC")
        .def("RemoveBond", &EditableMol::RemoveBond, "idx1"_a, "idx2"_a,
             R"DOC(Remove the specified bond from the molecule)DOC")

        .def("AddBond", &EditableMol::AddBond, "beginAtomIdx"_a, "endAtomIdx"_a,
             "order"_a = Bond::UNSPECIFIED,
             R"DOC(add a bond, returns the total number of bonds)DOC")

        .def("AddAtom", &EditableMol::AddAtom, "atom"_a,
             R"DOC(add an atom, returns the index of the newly added atom)DOC")

        .def("ReplaceAtom", &EditableMol::ReplaceAtom, "index"_a, "newAtom"_a,
             "updateLabel"_a = false, "preserveProps"_a = false,
             R"DOC(replaces the specified atom with the provided one
If updateLabel is True, the new atom becomes the active atom
If preserveProps is True preserve keep the existing props unless explicit set on the new atom)DOC")
        .def("ReplaceBond", &EditableMol::ReplaceBond, "index"_a, "newBond"_a,
             "preserveProps"_a = false,
             R"DOC(replaces the specified bond with the provided one.
If preserveProps is True preserve keep the existing props unless explicit set on the new bond)DOC")

        .def("BeginBatchEdit", &EditableMol::BeginBatchEdit,
             R"DOC(starts batch editing)DOC")
        .def("RollbackBatchEdit", &EditableMol::RollbackBatchEdit,
             R"DOC(cancels batch editing)DOC")
        .def("CommitBatchEdit", &EditableMol::CommitBatchEdit,
             R"DOC(finishes batch editing and makes the actual edits)DOC")

        .def("GetMol", &EditableMol::GetMol, nb::rv_policy::take_ownership,
             R"DOC(Returns a Mol (a normal molecule))DOC");
  };
};

}  // namespace RDKit
void wrap_EditableMol(nb::module_ &m) { RDKit::EditableMol_wrapper::wrap(m); }
