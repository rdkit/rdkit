// $Id$
//
//  Copyright (C) 2009 Greg Landrum
//
//   @@ All Rights Reserved  @@
//

#include <boost/python.hpp>
#include <RDGeneral/types.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

namespace python = boost::python;
using namespace RDKit;

namespace {
  PyObject* replaceSubstructures(const ROMol &orig,
                                 const ROMol &query,
                                 const ROMol &replacement,
                                 bool replaceAll=false) {
    std::vector<ROMOL_SPTR> v=replaceSubstructs(orig, query,
                                                replacement, replaceAll);
    PyObject *res=PyTuple_New(v.size());
    for(unsigned int i=0;i<v.size();++i){
      PyTuple_SetItem(res,i,
                      python::converter::shared_ptr_to_python(v[i]));
    }
    return res;
  }


}

BOOST_PYTHON_MODULE(rdchemtransforms)
{
  python::scope().attr("__doc__") =
    "Module containing RDKit functionality for doing chemical transformations on molecules."
    ;

  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  std::string docString;

      // ------------------------------------------------------------------------
      docString="Replaces sidechains in a molecule with dummy atoms for their attachment points.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
  RETURNS: a new molecule with the sidechains removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSidechains('CCC1CCC1','C1CCC1') -> '[Xa]C1CCC1'\n\
\n\
    - ReplaceSidechains('CCC1CC1','C1CCC1') -> ''\n\
\n\
    - ReplaceSidechains('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'\n\
\n";
      python::def("ReplaceSidechains", replaceSidechains,
                  (python::arg("mol"),python::arg("coreQuery")),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Removes the core of a molecule and labels the sidechains with dummy atoms.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
    - replaceDummies: toggles replacement of atoms that match dummies in the query\n\
\n\
  RETURNS: a new molecule with the core removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceCore('CCC1CCC1','C1CCC1') -> 'CC[Xa]'\n\
\n\
    - ReplaceCore('CCC1CC1','C1CCC1') -> ''\n\
\n\
    - ReplaceCore('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'\n\
\n\
    - ReplaceCore('C1CNCC1','N') -> '[Xa]CCCC[Xb]'\n\
\n\
    - ReplaceCore('C1CCC1CN','C1CCC1[*]',False) -> '[Xa]CN'\n\
\n";
      python::def("ReplaceCore", replaceCore,
                  (python::arg("mol"),python::arg("coreQuery"),
                   python::arg("replaceDummies")=true,
                   python::arg("labelByIndex")=false),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Removes atoms matching a substructure query from a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as a substructure query\n\
\n\
    - onlyFrags: (optional) if this toggle is set, atoms will only be removed if\n\
      the entire fragment in which they are found is matched by the query.\n\
      See below for examples.\n\
      Default value is 0 (remove the atoms whether or not the entire fragment matches)\n\
\n\
  RETURNS: a new molecule with the substructure removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - DeleteSubstructs('CCOC','OC') -> 'CC'\n\
\n\
    - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'\n\
\n\
    - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'\n\
\n\
    - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'\n\
\n";
      python::def("DeleteSubstructs", deleteSubstructs,
                  (python::arg("mol"),python::arg("query"),
                   python::arg("onlyFrags")=false),
                  docString.c_str(),
                  python::return_value_policy<python::manage_new_object>());

      // ------------------------------------------------------------------------
      docString="Replaces atoms matching a substructure query in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as a substructure query\n\
\n\
    - replacement: the molecule to be used as the replacement\n\
\n\
    - replaceAll: (optional) if this toggle is set, all substructures matching\n\
      the query will be replaced in a single result, otherwise each result will\n\
      contain a separate replacement.\n\
      Default value is False (return multiple replacements)\n\
\n\
  RETURNS: a tuple of new molecules with the substructures replaced removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSubstructs('CCOC','OC','NC') -> ('CCNC',)\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','NC') -> ('COCCNC','CNCCOC')\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','NC',True) -> ('CNCCNC',)\n\
\n";
      python::def("ReplaceSubstructs", replaceSubstructures,
                  (python::arg("mol"),python::arg("query"),
                   python::arg("replacement"),
                   python::arg("replaceAll")=false),
                  docString.c_str());
}


