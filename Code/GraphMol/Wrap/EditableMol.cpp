// $Id$
//
//  Copyright (C) 2007 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// there's a compiler bug in some versions of g++ that causes this file to not compile unless
// we skip the docstrings:
//#define BOOST_PYTHON_NO_PY_SIGNATURES

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

#include "rdchem.h"
// ours
#include <GraphMol/RDKitBase.h>

namespace python = boost::python;

namespace RDKit {

namespace {
  class EditableMol : boost::noncopyable {
  public:
    EditableMol(const ROMol &m) {
      dp_mol= new RWMol(m);
    };
    ~EditableMol() {
      PRECONDITION(dp_mol,"no molecule");
      delete dp_mol;
    };

    void RemoveAtom(unsigned int idx){
      PRECONDITION(dp_mol,"no molecule");
      dp_mol->removeAtom(idx);
    };
    void RemoveBond(unsigned int idx1,unsigned int idx2){
      PRECONDITION(dp_mol,"no molecule");
      dp_mol->removeBond(idx1,idx2);
    };
    int AddBond(unsigned int begAtomIdx,
                 unsigned int endAtomIdx,
                 Bond::BondType order=Bond::UNSPECIFIED)
    {
      PRECONDITION(dp_mol,"no molecule");
      return dp_mol->addBond(begAtomIdx,endAtomIdx,order);
    };
    int AddAtom(Atom *atom){
      PRECONDITION(dp_mol,"no molecule");
      PRECONDITION(atom,"bad atom");
      return dp_mol->addAtom(atom,true,false);
    };
    void ReplaceAtom(unsigned int idx,Atom *atom){
      PRECONDITION(dp_mol,"no molecule");
      PRECONDITION(atom,"bad atom");
      dp_mol->replaceAtom(idx,atom);
    };
    ROMol *GetMol() const{
      PRECONDITION(dp_mol,"no molecule");
      ROMol *res=new ROMol(*dp_mol);
      return res;
    };
  private:
    RWMol *dp_mol;
  };
}
  
struct EditableMol_wrapper {
  static void wrap(){
    std::string molClassDoc = "The EditableMol class.\n\n\
   This class can be used to add/remove bonds and atoms to\n\
   a molecule.\n\
   In order to use it, you need to first construct an EditableMol\n\
   from a standard Mol:\n\
   >>> m = Chem.MolFromSmiles('CCC')\n\
   >>> em = Chem.EditableMol(m)\n\
   >>> em.AddAtom(Chem.Atom(8))\n\
   >>> em.AddBond(0,3,Chem.BondType.SINGLE)\n\
   >>> m2 = em.GetMol()\n\
   >>> Chem.SanitizeMol(m2)\n\
   >>> Chem.MolToSmiles(m2)\n\
   'CCCO'\n\
\n\
   *Note*: It is very, very easy to shoot yourself in the foot with\n\
           this class by constructing an unreasonable molecule.\n\
";
    python::class_<EditableMol,boost::noncopyable >("EditableMol",
        "an editable molecule class",
        python::init<const ROMol &>("Construct from a Mol"))
      .def("RemoveAtom",&EditableMol::RemoveAtom,
      "Remove the specified atom from the molecule")
      .def("RemoveBond",&EditableMol::RemoveBond,
      "Remove the specified bond from the molecule")
      
      .def("AddBond",&EditableMol::AddBond,
                     (python::arg("mol"),python::arg("beginAtomIdx"),python::arg("endAtomIdx"),
                      python::arg("order")=Bond::UNSPECIFIED),
      "add a bond, returns the index of the newly added bond")
      
      .def("AddAtom",&EditableMol::AddAtom,
                     (python::arg("mol"),python::arg("atom")),
      "add an atom, returns the index of the newly added atom")
      .def("ReplaceAtom",&EditableMol::ReplaceAtom,
                     (python::arg("mol"),python::arg("index"),python::arg("newAtom")),
      "replaces the specified atom with the provided one")
      .def("GetMol",&EditableMol::GetMol,
           "Returns a Mol (a normal molecule)",
           python::return_value_policy<python::manage_new_object>())
      ;

  };
};

}// end of namespace
void wrap_EditableMol() {
  RDKit::EditableMol_wrapper::wrap();
}
