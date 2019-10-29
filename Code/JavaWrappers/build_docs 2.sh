#! /bin/bash -x
python parse_doxy_html.py atom.doxytext RDKit::Atom > Atom_doc.i
python parse_doxy_html.py bond.doxytext RDKit::Bond > Bond_doc.i
python parse_doxy_html.py bonditerators.doxytext RDKit::BondIterator_ >  BondIterators_doc.i
python parse_doxy_html.py chemreactions.doxytext RDKit::ChemicalReaction > ChemReactions_doc.i
python parse_doxy_html.py conformer.doxytext RDKit::Conformer > Conformer_doc.i
python parse_doxy_html.py explicitbitvect.doxytext RDKit::ExplicitBitVect > ExplicitBitVect_doc.i
python parse_doxy_html.py molops.doxytext RDKit::MolOps > RDKFuncs_doc.i
python parse_doxy_html.py periodictable.doxytext RDKit::PeriodicTable > PeriodicTable_doc.i
python parse_doxy_html.py queryatom.doxytext RDKit::QueryAtom > QueryAtom_doc.i
python parse_doxy_html.py point.doxytext RDGeom::Point3D > Point3D_doc.i
python parse_doxy_html.py querybond.doxytext RDKit::QueryBond > QueryBond_doc.i
python parse_doxy_html.py ringinfo.doxytext RDKit::RingInfo > RingInfo_doc.i
python parse_doxy_html.py romol.doxytext RDKit::ROMol > ROMol_doc.i
python parse_doxy_html.py rwmol.doxytext RDKit::RWMol > RWMol_doc.i
python parse_doxy_html.py sdmolsupplier.doxytext RDKit::SDMolSupplier > SDMolSupplier_doc.i
python parse_doxy_html.py smilesmolsupplier.doxytext RDKit::SmilesMolSupplier > SmilesMolSupplier_doc.i
python parse_doxy_html.py smileswriter.doxytext RDKit::SmilesWriter > SmilesWriter_doc.i
python parse_doxy_html.py tdtmolsupplier.doxytext RDKit::TDTMolSupplier > TDTMolSupplier_doc.i
python parse_doxy_html.py tdtwriter.doxytext RDKit::TDTWriter > TDTWriter_doc.i
python parse_doxy_html.py transform2d.doxytext RDGeom::Transform2D > Transform2D_doc.i
python parse_doxy_html.py transform3d.doxytext RDGeom::Transform3D > Transform3D_doc.i

