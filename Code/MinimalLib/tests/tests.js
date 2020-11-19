//
//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

const assert = require('assert');
var initRDKitModule = require("../demo/RDKit_minimal.js");
var RDKitModule;

// the goal here isn't to be comprehensive (the RDKit has tests for that),
// just to make sure that the wrappers are working as expected
function test_basics(){
    var bmol = RDKitModule.get_mol("c1ccccc");
    assert.equal(bmol.is_valid(),0);
    
    var mol = RDKitModule.get_mol("c1ccccc1O");
    assert.equal(mol.is_valid(),1);
    assert.equal(mol.get_smiles(),"Oc1ccccc1");
    assert.equal(mol.get_inchi(),"InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H");
    assert.equal(RDKitModule.get_inchikey_for_inchi(mol.get_inchi()),"ISWSIDIOOBJBQZ-UHFFFAOYSA-N");

    var mb = mol.get_molblock();
    assert(mb.search("M  END")>0);
    var mol2 = RDKitModule.get_mol(mb);
    assert.equal(mol2.is_valid(),1);
    assert.equal(mol2.get_smiles(),"Oc1ccccc1");
    
    
    var descrs = JSON.parse(mol.get_descriptors());
    assert.equal(descrs.NumAromaticRings,1);
    assert.equal(descrs.NumRings,1);
    assert.equal(descrs.amw,94.11299);

    var fp1 = mol.get_morgan_fp();
    assert.equal(fp1.length,2048);
    assert.equal((fp1.match(/1/g)||[]).length,11);
    var fp2 = mol.get_morgan_fp(0,512);
    assert.equal(fp2.length,512);
    assert.equal((fp2.match(/1/g)||[]).length,3);
    
    var svg = mol.get_svg();
    assert(svg.search("svg")>0);

    var qmol = RDKitModule.get_qmol("Oc(c)c");
    assert.equal(qmol.is_valid(),1);
    var match = mol.get_substruct_match(qmol);
    var pmatch = JSON.parse(match);
    assert.equal(pmatch.atoms.length,4);
    assert.equal(pmatch.atoms[0],6);
    var svg2 = mol.get_svg_with_highlights(match);
    assert(svg2.search("svg")>0);
    assert(svg.search("#FF7F7F")<0);
    assert(svg2.search("#FF7F7F")>0);
}

function test_sketcher_services(){
    var mol = RDKitModule.get_mol("C[C@](F)(Cl)/C=C/C(F)Br");
    assert.equal(mol.is_valid(),1);
    var tags = mol.get_stereo_tags();
    assert.equal(tags,'{"CIP_atoms":[[1,"(S)"],[6,"(?)"]],"CIP_bonds":[[4,5,"(E)"]]}');
}

function test_sketcher_services2(){
    var mol = RDKitModule.get_mol("c1ccccc1");
    assert.equal(mol.is_valid(),1);
    var molb = mol.add_hs();
    assert(molb.search(" H ")>0);
    assert.equal((molb.match(/ H /g) || []).length,6);

    var mol2 = RDKitModule.get_mol(molb);
    assert.equal(mol2.is_valid(),1);
    var molb2 = mol2.get_molblock();
    assert(molb2.search(" H ")>0); 
    assert.equal((molb2.match(/ H /g) || []).length,6);

    molb2 = mol2.remove_hs();
    assert(molb2.search(" H ")<0); 
}


function test_abbreviations(){
    var bmol = RDKitModule.get_mol("C1CCC1C(F)(F)F");
    assert.equal(bmol.is_valid(),1);
    bmol.condense_abbreviations();
    assert.equal(bmol.get_cxsmiles(),"FC(F)(F)C1CCC1");
    bmol.condense_abbreviations(1.0,false);
    assert.equal(bmol.get_cxsmiles(),"*C1CCC1 |$CF3;;;;$|");
}


function test_generate_aligned_coords(){
    var smiles = "CCC";
    var mol = RDKitModule.get_mol(smiles);
    var template = "CC";
    var qmol = RDKitModule.get_mol(template);
    assert.equal(mol.generate_aligned_coords(qmol, true), "");
}


initRDKitModule().then(function(instance) {
    RDKitModule = instance;
    console.log(RDKitModule.version());
    test_basics();
    test_sketcher_services();
    test_sketcher_services2();
    test_abbreviations();
    test_generate_aligned_coords();
    console.log("Tests finished successfully");
});
