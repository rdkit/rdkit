//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmoldraw2d_array_API
#include <boost/python.hpp>
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <Geometry/point.h>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

namespace python = boost::python;

namespace RDKit {
  void drawMoleculeHelper1(MolDraw2D &self,
                          const ROMol &mol ,
                          python::object highlight_atoms,
                          python::object highlight_atom_map,
                          int confId=-1){
    std::vector<int> *highlightAtoms=pythonObjectToVect(highlight_atoms,static_cast<int>(mol.getNumAtoms()));
    // FIX: support these
    std::map<int,DrawColour> *ham=NULL;

    self.drawMolecule(mol,highlightAtoms,ham,confId);
    
    delete highlightAtoms;
    delete ham;
  }
  void drawMoleculeHelper2(MolDraw2D &self,
                          const ROMol &mol ,
                          python::object highlight_atoms,
                          python::object highlight_bonds,
                          python::object highlight_atom_map,
                          python::object highlight_bond_map,
                          int confId=-1){
    std::vector<int> *highlightAtoms=pythonObjectToVect(highlight_atoms,static_cast<int>(mol.getNumAtoms()));
    std::vector<int> *highlightBonds=pythonObjectToVect(highlight_bonds,static_cast<int>(mol.getNumBonds()));
    // FIX: support these
    std::map<int,DrawColour> *ham=NULL;
    std::map<int,DrawColour> *hbm=NULL;

    self.drawMolecule(mol,highlightAtoms,highlightBonds,ham,hbm,confId);
    
    delete highlightAtoms;
    delete highlightBonds;
    delete ham;
    delete hbm;
  }
}

BOOST_PYTHON_MODULE(rdMolDraw2D) {
  
  python::scope().attr("__doc__") =
    "Module containing a C++ implementation of 2D molecule drawing";
   
  std::string docString;

  python::class_<std::map<int,std::string> >("IntStringMap")
    .def(python::map_indexing_suite<std::map<int,std::string>, true >())
    ;


  docString="Drawing options";
  python::class_<RDKit::MolDrawOptions,boost::noncopyable>("MolDrawOptions",docString.c_str())
    .def_readwrite("dummiesAreAttachments",&RDKit::MolDrawOptions::dummiesAreAttachments)
    .def_readwrite("circleAtoms",&RDKit::MolDrawOptions::circleAtoms)
    //.def_readwrite("highlightColour",&RDKit::MolDrawOptions::highlightColour)
    .def_readwrite("atomLabels",&RDKit::MolDrawOptions::atomLabels,"maps indices to atom labels")
    .def_readwrite("continuousHighlight",&RDKit::MolDrawOptions::continuousHighlight)
    .def_readwrite("flagCloseContactsDist",&RDKit::MolDrawOptions::flagCloseContactsDist)

    ;
  docString="Drawer abstract base class";
  python::class_<RDKit::MolDraw2D,boost::noncopyable>("MolDraw2D",docString.c_str(),python::no_init)
    .def("DrawMolecule",RDKit::drawMoleculeHelper1,
         (python::arg("self"),python::arg("mol"),
          python::arg("highlightAtoms")=python::object(),
          python::arg("highlightAtomMap")=python::object(),
          python::arg("confId")=-1
          ),
         "renders a molecule\n")
    .def("DrawMolecule",RDKit::drawMoleculeHelper2,
         (python::arg("self"),python::arg("mol"),
          python::arg("highlightAtoms"),
          python::arg("highlightBonds"),
          python::arg("highlightAtomMap")=python::object(),
          python::arg("highlightBondMap")=python::object(),
          python::arg("confId")=-1
          ),
         "renders a molecule\n")
    .def("drawOptions",(RDKit::MolDrawOptions &(RDKit::MolDraw2D::*)())&RDKit::MolDraw2D::drawOptions,
         python::return_internal_reference<1,
         python::with_custodian_and_ward_postcall<0,1> >(),
         "Returns a modifiable version of the current drawing options"
         )
    ;
  docString="SVG molecule drawer";
  python::class_<RDKit::MolDraw2DSVG, python::bases<RDKit::MolDraw2D>, boost::noncopyable >("MolDraw2DSVG",docString.c_str(),
                                                                                            python::init<int,int>())
    .def("FinishDrawing",&RDKit::MolDraw2DSVG::finishDrawing,
         "add the last bits of SVG to finish the drawing")
    .def("GetDrawingText",&RDKit::MolDraw2DSVG::getDrawingText,
         "return the SVG")
    ;


}
