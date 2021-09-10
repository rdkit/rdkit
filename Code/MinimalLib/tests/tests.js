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
const {
    performance
  } = require('perf_hooks');
var initRDKitModule = require("../demo/RDKit_minimal.js");
var RDKitModule;
const fs       = require('fs');
const readline = require('readline');

// the goal here isn't to be comprehensive (the RDKit has tests for that),
// just to make sure that the wrappers are working as expected
function test_basics() {
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
    
    var mjson = mol.get_json();
    assert(mjson.search("commonchem")>0);
    var mol3 = RDKitModule.get_mol(mjson);
    assert.equal(mol3.is_valid(),1);
    assert.equal(mol3.get_smiles(),"Oc1ccccc1");
    
    var descrs = JSON.parse(mol.get_descriptors());
    assert.equal(descrs.NumAromaticRings,1);
    assert.equal(descrs.NumRings,1);
    assert.equal(descrs.amw,94.11299);

    var checkStringBinaryFpIdentity = (stringFp, binaryFp) => {
        assert.equal(binaryFp.length, stringFp.length / 8);
        for (var i = 0, c = 0; i < binaryFp.length; ++i) {
            var byte = 0;
            for (var j = 0; j < 8; ++j, ++c) {
                if (stringFp[c] === "1") {
                    byte |= (1 << j);
                }
            }
            assert.equal(byte, binaryFp[i]);
        }
    };

    var fp1 = mol.get_morgan_fp();
    assert.equal(fp1.length,2048);
    assert.equal((fp1.match(/1/g)||[]).length,11);
    var fp1Uint8Array = mol.get_morgan_fp_as_uint8array();
    checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
    var fp2 = mol.get_morgan_fp(0,512);
    assert.equal(fp2.length,512);
    assert.equal((fp2.match(/1/g)||[]).length,3);
    var fp2Uint8Array = mol.get_morgan_fp_as_uint8array(0, 512);
    checkStringBinaryFpIdentity(fp2, fp2Uint8Array);
    
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

function test_molblock_nostrict() {
    var molblock = `
  MJ201100                      

 10 10  0  0  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  8  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
M  STY  1   1 SUP
M  SAL   1  4   7   8   9  10
M  SMT   1 CF3
M  SBL   1  1   7
M  SAP   1  1   8
M  END`;
    var mol = RDKitModule.get_mol(molblock);
    assert.equal(mol.is_valid(),1);
    var mb = mol.get_molblock();
    assert.equal(mb.includes("M  SAP   1  1   8   6"), true);
    var qmol = RDKitModule.get_qmol(molblock);
    assert.equal(qmol.is_valid(),1);
    var qmb = qmol.get_molblock();
    assert.equal(qmb.includes("M  SAP   1  1   8   6"), true);
}

function test_molblock_rgp() {
    var molblock = `
  MJ190400                      

  9  9  0  0  0  0  0  0  0  0999 V2000
   -6.5623    0.3105    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -5.8478   -0.1019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1333    0.3105    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4188   -0.1019    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4188   -0.9269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1333   -1.3394    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8478   -0.9269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7043   -1.3394    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9898   -0.9268    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  7  4  0  0  0  0
  2  3  4  0  0  0  0
  2  7  4  0  0  0  0
  5  8  1  0  0  0  0
  8  9  1  0  0  0  0
  1  2  1  0  0  0  0
M  RGP  2   1   2   9   1
M  END`;
    var mol = RDKitModule.get_mol(molblock);
    assert.equal(mol.is_valid(),1);
}

function test_sketcher_services() {
    var mol = RDKitModule.get_mol("C[C@](F)(Cl)/C=C/C(F)Br");
    assert.equal(mol.is_valid(),1);
    var tags = mol.get_stereo_tags();
    assert.equal(tags,'{"CIP_atoms":[[1,"(S)"],[6,"(?)"]],"CIP_bonds":[[4,5,"(E)"]]}');
}

function test_sketcher_services2() {
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


function test_abbreviations() {
    var bmol = RDKitModule.get_mol("C1CCC1C(F)(F)F");
    assert.equal(bmol.is_valid(),1);
    bmol.condense_abbreviations();
    assert.equal(bmol.get_cxsmiles(),"FC(F)(F)C1CCC1");
    bmol.condense_abbreviations(1.0,false);
    assert.equal(bmol.get_cxsmiles(),"*C1CCC1 |$CF3;;;;$|");
}

function test_substruct_library(done) {
    done.test_substruct_library = false;
    var smiReader = readline.createInterface({
        input: fs.createReadStream(__dirname + '/../../GraphMol/test_data/compounds.smi')
    });
    var sslib = new RDKitModule.SubstructLibrary();
    // var t0 = performance.now()
    // console.log('Started adding trusted SMILES');
    smiReader.on('line', (smi) => {
        sslib.add_trusted_smiles(smi);
    });
    smiReader.on('close', () => {
        var query = RDKitModule.get_qmol("N");
        // var t1 = performance.now();
        // console.log('Finished adding trusted SMILES took ' + (t1 - t0) / 1000 + ' seconds');
        assert.equal(sslib.count_matches(query), 52);
        assert.equal(sslib.get_matches(query), JSON.stringify([
            12,13,19,22,24,30,31,32,35,36,39,41,43,44,55,56,58,64,72,80,
            85,95,96,101,105,113,124,127,128,131,143,150,151,185,201,202,
            203,214,215,223,232,234,238,240,241,246,258,261,263,265,266,284
        ]));
        done.test_substruct_library = true;
    });
}

function test_substruct_library_merge_hs() {
    var sslib = new RDKitModule.SubstructLibrary();
    var mol1 = RDKitModule.get_mol("c1ccccc1");
    var mol2 = RDKitModule.get_mol("Cc1ccccc1");
    sslib.add_trusted_smiles(mol1.get_smiles());
    sslib.add_trusted_smiles(mol2.get_smiles());
    var query = RDKitModule.get_mol(`
  MJ201100          2D          

  6  6  0  0  0  0  0  0  0  0999 V2000
   -1.0491    0.7134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7635    0.3009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7635   -0.5241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0491   -0.9366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3346   -0.5241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3346    0.3009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  5  6  2  0  0  0  0
  1  2  2  0  0  0  0
  6  1  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
M  END`);
    assert.equal(sslib.get_matches(query), JSON.stringify([0, 1]));
    var query_hs = RDKitModule.get_mol(query.add_hs());
    assert.equal(sslib.get_matches(query_hs), JSON.stringify([]));
    query_hs = RDKitModule.get_mol(query_hs.get_molblock(), JSON.stringify({ mergeQueryHs: true }));
    assert.equal(sslib.get_matches(query_hs), JSON.stringify([0]));
}

function test_substruct_library_empty_mols() {
    var sslib = new RDKitModule.SubstructLibrary();
    var mol1 = RDKitModule.get_mol("");
    var mol2 = RDKitModule.get_mol("");
    sslib.add_trusted_smiles(mol1.get_smiles());
    sslib.add_trusted_smiles(mol2.get_smiles());
    var query = RDKitModule.get_mol("C");
    assert.equal(sslib.get_matches(query), JSON.stringify([]));
    var empty_query = RDKitModule.get_mol("");
    assert.equal(sslib.get_matches(empty_query), JSON.stringify([]));
}

function test_substruct_library_empty_query() {
    var sslib = new RDKitModule.SubstructLibrary();
    var mol1 = RDKitModule.get_mol("C");
    var mol2 = RDKitModule.get_mol("CC");
    sslib.add_trusted_smiles(mol1.get_smiles());
    sslib.add_trusted_smiles(mol2.get_smiles());
    var query = RDKitModule.get_mol("");
    assert.equal(sslib.get_matches(query), JSON.stringify([]));
}

function test_substruct_library_empty_lib() {
    var sslib = new RDKitModule.SubstructLibrary();
    var query = RDKitModule.get_mol("C");
    assert.equal(sslib.get_matches(query), JSON.stringify([]));
    var empty_query = RDKitModule.get_mol("");
    assert.equal(sslib.get_matches(empty_query), JSON.stringify([]));
}

function test_generate_aligned_coords() {
    var smiles = "CCC";
    var mol = RDKitModule.get_mol(smiles);
    var template = "CC";
    var qmol = RDKitModule.get_mol(template);
    assert.equal(mol.generate_aligned_coords(qmol, true), "");
}


function test_isotope_labels() {
    var mol = RDKitModule.get_mol("[1*]c1cc([2*])c([3*])c[14c]1");
    assert.equal(mol.is_valid(), 1);

    var textIsoDummyIso = mol.get_svg_with_highlights(JSON.stringify({}));
    var nLinesIsoDummyIso = textIsoDummyIso.split("\n").length;

    var textNoIsoDummyIso = mol.get_svg_with_highlights(JSON.stringify({ isotopeLabels: false }));
    var nLinesNoIsoDummyIso = textNoIsoDummyIso.split("\n").length;

    var textIsoNoDummyIso = mol.get_svg_with_highlights(JSON.stringify({ dummyIsotopeLabels: false }));
    var nLinesIsoNoDummyIso = textIsoNoDummyIso.split("\n").length;

    var textNoIsoNoDummyIso = mol.get_svg_with_highlights(JSON.stringify({ isotopeLabels: false, dummyIsotopeLabels: false }));
    var nLinesNoIsoNoDummyIso = textNoIsoNoDummyIso.split("\n").length;

    var res = [nLinesNoIsoNoDummyIso, nLinesIsoNoDummyIso, nLinesNoIsoDummyIso, nLinesIsoDummyIso];
    var resSorted = [...res];
    resSorted.sort((a, b) => (a - b));
    assert.ok(res.every((resItem, i) => (resItem === resSorted[i])));
}


function test_generate_aligned_coords_allow_rgroups() {
    var template_molblock = `
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  RGP  3   7   1   8   2   9   3
M  END`;
    var ortho_meta_smiles = 'c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1';
    var ortho_smiles = 'c1ccc(-c2ccccc2-c2ccccc2)cc1';
    var meta_smiles = 'c1ccc(-c2cccc(-c3ccccc3)c2)cc1';
    var biphenyl_smiles = 'c1ccccc1-c1ccccc1';
    var phenyl_smiles = 'c1ccccc1';
    var template_ref = RDKitModule.get_mol(template_molblock);
    var ortho_meta = RDKitModule.get_mol(ortho_meta_smiles);
    var ortho = RDKitModule.get_mol(ortho_smiles);
    var meta = RDKitModule.get_mol(meta_smiles);
    var biphenyl = RDKitModule.get_mol(biphenyl_smiles);
    var phenyl = RDKitModule.get_mol(phenyl_smiles);
    assert.equal(JSON.parse(ortho_meta.generate_aligned_coords(
        template_ref, false, true)).atoms.length, 9);
    assert.equal(JSON.parse(ortho.generate_aligned_coords(
        template_ref, false, true)).atoms.length, 8);
    assert.equal(JSON.parse(meta.generate_aligned_coords(
        template_ref, false, true)).atoms.length, 8);
    assert.equal(JSON.parse(biphenyl.generate_aligned_coords(
        template_ref, false, true)).atoms.length, 7);
    assert.equal(JSON.parse(phenyl.generate_aligned_coords(
        template_ref, false, true)).atoms.length, 6);
}


function test_accept_failure() {
    var template_molblock = `
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  RGP  3   7   1   8   2   9   3
M  END
`;
    var mol_molblock = `
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  END
`;
    var template_ref = RDKitModule.get_mol(template_molblock);
    var mol = RDKitModule.get_mol(mol_molblock);
    var hasThrown = false;
    try {
        mol.generate_aligned_coords(template_ref, false, true, false);
    } catch (e) {
        hasThrown = true;
    }
    assert(hasThrown);
    assert.equal(mol.get_molblock(), mol_molblock);
}

function test_get_mol_no_kekulize() {
    var mol = RDKitModule.get_mol("c");
    assert(!mol.is_valid());
    mol = RDKitModule.get_mol("c", JSON.stringify({kekulize: false}));
    assert(mol.is_valid());
}


initRDKitModule().then(function(instance) {
    var done = {};
    const waitAllTestsFinished = () => {
        const poll = resolve => {
            if (Object.values(done).every(v => v)) {
                resolve();
            } else {
                setTimeout(() => poll(resolve), 100);
            }
        }
        return new Promise(poll);
    }
    RDKitModule = instance;
    console.log(RDKitModule.version());
    test_basics();
    test_molblock_nostrict();
    test_molblock_rgp();
    test_sketcher_services();
    test_sketcher_services2();
    test_abbreviations();
    test_substruct_library(done);
    test_substruct_library_merge_hs();
    test_substruct_library_empty_mols();
    test_substruct_library_empty_lib();
    test_substruct_library_empty_query();
    test_generate_aligned_coords();
    test_isotope_labels();
    test_generate_aligned_coords_allow_rgroups();
    test_accept_failure();
    test_get_mol_no_kekulize();
    waitAllTestsFinished().then(() =>
        console.log("Tests finished successfully")
    );
});
