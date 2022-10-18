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
    var bmol = null;
    try {
        bmol = RDKitModule.get_mol("c1ccccc");
        assert.equal(bmol.is_valid(),0);
    } catch {
        assert(bmol === null);
    }
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

    {
        var fp1 = mol.get_morgan_fp();
        assert.equal(fp1.length, 2048);
        assert.equal((fp1.match(/1/g)||[]).length, 11);
        var fp1Uint8Array = mol.get_morgan_fp_as_uint8array();
        checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
        var fp2 = mol.get_morgan_fp(JSON.stringify({ radius: 0, nBits: 512 }));
        assert.equal(fp2.length, 512);
        assert.equal((fp2.match(/1/g)||[]).length, 3);
        var fp2Uint8Array = mol.get_morgan_fp_as_uint8array(JSON.stringify({ radius: 0, nBits: 512 }));
        checkStringBinaryFpIdentity(fp2, fp2Uint8Array);
        var fp3 = mol.get_morgan_fp(0, 512);
        assert.equal(fp3.length, 512);
        assert.equal((fp3.match(/1/g)||[]).length, 3);
        var fp3Uint8Array = mol.get_morgan_fp_as_uint8array(0, 512);
        checkStringBinaryFpIdentity(fp3, fp3Uint8Array);
    }

    {
        var fp1 = mol.get_pattern_fp();
        assert.equal(fp1.length, 2048);
        assert.equal((fp1.match(/1/g)||[]).length, 73);
        var fp1Uint8Array = mol.get_pattern_fp_as_uint8array();
        checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
        var fp2 = mol.get_pattern_fp(JSON.stringify({ nBits: 256 }));
        assert.equal(fp2.length, 256);
        assert.equal((fp2.match(/1/g)||[]).length, 65);
        var fp2Uint8Array = mol.get_pattern_fp_as_uint8array(JSON.stringify({ nBits: 256 }));
        checkStringBinaryFpIdentity(fp2, fp2Uint8Array);
        var fp3 = mol.get_pattern_fp(256);
        assert.equal(fp3.length, 256);
        assert.equal((fp3.match(/1/g)||[]).length, 65);
        var fp3Uint8Array = mol.get_pattern_fp_as_uint8array(256);
        checkStringBinaryFpIdentity(fp3, fp3Uint8Array);
    }

    {
        var fp1 = mol.get_topological_torsion_fp();
        assert.equal(fp1.length, 2048);
        assert.equal((fp1.match(/1/g)||[]).length, 8);
        var fp1Uint8Array = mol.get_topological_torsion_fp_as_uint8array();
        checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
        var fp2 = mol.get_topological_torsion_fp(JSON.stringify({ nBits: 512 }));
        assert.equal(fp2.length, 512);
        assert.equal((fp2.match(/1/g)||[]).length, 8);
        var fp2Uint8Array = mol.get_topological_torsion_fp_as_uint8array(JSON.stringify({ nBits: 512 }));
        checkStringBinaryFpIdentity(fp2, fp2Uint8Array);
    }

    {
        var fp1 = mol.get_rdkit_fp();
        assert.equal(fp1.length, 2048);
        assert.equal((fp1.match(/1/g)||[]).length, 38);
        var fp1Uint8Array = mol.get_rdkit_fp_as_uint8array();
        checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
        var fp2 = mol.get_rdkit_fp(JSON.stringify({ nBits: 512 }));
        assert.equal(fp2.length, 512);
        assert.equal((fp2.match(/1/g)||[]).length, 37);
        var fp2Uint8Array = mol.get_rdkit_fp_as_uint8array(JSON.stringify({ nBits: 512 }));
        checkStringBinaryFpIdentity(fp2, fp2Uint8Array);
        var fp3 = mol.get_rdkit_fp(JSON.stringify({ nBits: 512, minPath: 2 }));
        assert.equal(fp3.length, 512);
        assert.equal((fp3.match(/1/g)||[]).length, 34);
        var fp3Uint8Array = mol.get_rdkit_fp_as_uint8array(JSON.stringify({ nBits: 512, minPath: 2 }));
        checkStringBinaryFpIdentity(fp3, fp3Uint8Array);
        var fp4 = mol.get_rdkit_fp(JSON.stringify({ nBits: 512, minPath: 2, maxPath: 4 }));
        assert.equal(fp4.length, 512);
        assert.equal((fp4.match(/1/g)||[]).length, 16);
        var fp4Uint8Array = mol.get_rdkit_fp_as_uint8array(JSON.stringify({ nBits: 512, minPath: 2, maxPath: 4 }));
        checkStringBinaryFpIdentity(fp4, fp4Uint8Array);
    }

    {
        var fp1 = mol.get_atom_pair_fp();
        assert.equal(fp1.length, 2048);
        assert.equal((fp1.match(/1/g)||[]).length, 19);
        var fp1Uint8Array = mol.get_atom_pair_fp_as_uint8array();
        checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
        var fp2 = mol.get_atom_pair_fp(JSON.stringify({ nBits: 512 }));
        assert.equal(fp2.length, 512);
        assert.equal((fp2.match(/1/g)||[]).length, 19);
        var fp2Uint8Array = mol.get_atom_pair_fp_as_uint8array(JSON.stringify({ nBits: 512 }));
        checkStringBinaryFpIdentity(fp2, fp2Uint8Array);
        var fp3 = mol.get_atom_pair_fp(JSON.stringify({ nBits: 512, minLength: 2 }));
        assert.equal(fp3.length, 512);
        assert.equal((fp3.match(/1/g)||[]).length, 13);
        var fp3Uint8Array = mol.get_atom_pair_fp_as_uint8array(JSON.stringify({ nBits: 512, minLength: 2 }));
        checkStringBinaryFpIdentity(fp3, fp3Uint8Array);
        var fp4 = mol.get_atom_pair_fp(JSON.stringify({ nBits: 512, minLength: 2, maxLength: 3 }));
        assert.equal(fp4.length, 512);
        assert.equal((fp4.match(/1/g)||[]).length, 12);
        var fp4Uint8Array = mol.get_atom_pair_fp_as_uint8array(JSON.stringify({ nBits: 512, minLength: 2, maxLength: 3 }));
        checkStringBinaryFpIdentity(fp4, fp4Uint8Array);
    }

    if (typeof Object.getPrototypeOf(mol).get_avalon_fp === 'function') {
        var fp1 = mol.get_avalon_fp();
        assert.equal(fp1.length, 512);
        assert.equal((fp1.match(/1/g)||[]).length, 19);
        var fp1Uint8Array = mol.get_avalon_fp_as_uint8array();
        checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
        var fp2 = mol.get_avalon_fp(JSON.stringify({ nBits: 2048 }));
        assert.equal(fp2.length, 2048);
        assert.equal((fp2.match(/1/g)||[]).length, 20);
        var fp2Uint8Array = mol.get_avalon_fp_as_uint8array(JSON.stringify({ nBits: 2048 }));
        checkStringBinaryFpIdentity(fp2, fp2Uint8Array);
    }

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

function test_get_aromatic_kekule_form() {
    const aromRegExp = /  \d  \d  4  \d\n/g;
    const kekRegExp = /  \d  \d  [12]  \d\n/g;
    var mol = RDKitModule.get_mol("c1ccccc1");
    var molblock = mol.get_molblock();
    assert (molblock.match(aromRegExp) === null);
    assert (molblock.match(kekRegExp).length === 6);
    var molblock = mol.get_kekule_form();
    assert (molblock.match(aromRegExp) === null);
    assert (molblock.match(kekRegExp).length === 6);
    molblock = mol.get_aromatic_form();
    assert (molblock.match(aromRegExp).length === 6);
    assert (molblock.match(kekRegExp) === null);
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
    var bmol = RDKitModule.get_mol("C1CC(C(F)(F)F)C1");
    assert.equal(bmol.is_valid(),1);
    var mapping = bmol.condense_abbreviations();
    assert.equal(mapping, JSON.stringify({
        atoms: [0, 1, 2, 3, 4, 5, 6, 7],
        bonds: [0, 1, 2, 3, 4, 5, 6, 7],
    }));
    assert.equal(bmol.get_cxsmiles(),"FC(F)(F)C1CCC1");
    mapping = bmol.condense_abbreviations(1.0,false);
    assert.equal(mapping, JSON.stringify({
        atoms: [0, 1, 2, 3, 7],
        bonds: [0, 1, 2, 6, 7],
    }));
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
    [ true, false ].forEach(alignOnly => {
        var opts = JSON.stringify({ useCoordGen: false, allowRGroups: true, alignOnly })
        assert.equal(JSON.parse(ortho_meta.generate_aligned_coords(
            template_ref, opts)).atoms.length, 9);
        assert.equal(JSON.parse(ortho.generate_aligned_coords(
            template_ref, opts)).atoms.length, 8);
        assert.equal(JSON.parse(meta.generate_aligned_coords(
            template_ref, opts)).atoms.length, 8);
        assert.equal(JSON.parse(biphenyl.generate_aligned_coords(
            template_ref, opts)).atoms.length, 7);
        assert.equal(JSON.parse(phenyl.generate_aligned_coords(
            template_ref, opts)).atoms.length, 6);
    });
}

function test_generate_aligned_coords_accept_failure() {
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
    var res = mol.generate_aligned_coords(template_ref, false, true, false);
    assert(res === "");
    assert.equal(mol.get_molblock(), mol_molblock);
    res = mol.generate_aligned_coords(template_ref, JSON.stringify({
        useCoordGen: false,
        allowRGroups: true,
        acceptFailure: false,
    }));
    assert(res === "");
    assert.equal(mol.get_molblock(), mol_molblock);
    res = mol.generate_aligned_coords(template_ref, JSON.stringify({
        useCoordGen: false,
        allowRGroups: true,
        acceptFailure: true,
    }));
    assert(res === "{}");
    assert.notEqual(mol.get_molblock(), mol_molblock);
    var molNoCoords = RDKitModule.get_mol(mol.get_smiles());
    assert(!molNoCoords.has_coords());
    res = molNoCoords.generate_aligned_coords(template_ref, JSON.stringify({
        useCoordGen: false,
        allowRGroups: true,
        acceptFailure: false,
    }));
    assert(res === "");
    assert(!molNoCoords.has_coords());
    res = molNoCoords.generate_aligned_coords(template_ref, JSON.stringify({
        useCoordGen: false,
        allowRGroups: true,
        acceptFailure: true,
    }));
    assert(res === "{}");
    assert(molNoCoords.has_coords());
}

function test_generate_aligned_coords_align_only() {
    const template_molblock = `
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
  -13.7477    6.3036    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  -13.7477    4.7567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.6540    3.6628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.7477    2.5691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.8414    3.6628    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  -11.1071    3.6628    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  2  5  1  0
  3  6  1  0
M  RGP  2   1   1   6   2
M  END
`;
    const mol_molblock = `
     RDKit          2D

 18 22  0  0  0  0  0  0  0  0999 V2000
    4.3922   -1.5699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9211   -2.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5995   -0.5349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3731    0.8046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8441    1.2825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0704   -0.0568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8666    0.7748    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7736   -0.3197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7749   -1.8666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7718   -1.8679    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7731   -0.3208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8679    0.7718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0718   -0.0598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3933   -1.5729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9222   -2.0509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6008   -0.5379    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3744    0.8016    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8454    1.2795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  9 10  1  0
 11 10  1  0
 11  8  1  0
  8  9  1  0
  4  5  1  0
  6  5  1  0
  7  6  1  0
  3  4  1  0
  3  7  1  0
  1  6  1  0
  2  3  1  0
  2  1  1  0
 17 18  1  0
 13 18  1  0
 12 13  1  0
 16 17  1  0
 16 12  1  0
 14 13  1  0
 15 16  1  0
 15 14  1  0
 12 11  1  0
  8  7  1  0
M  END
`;
    const getCoordArray = (mol) => {
        const molblock = mol.get_molblock().split("\n");
        const numHeavyAtoms = parseInt(molblock[3].substr(0, 3));
        return molblock.filter((_, idx) => idx > 3 && idx < 4 + numHeavyAtoms).map(
            line => [0, 1, 2].map(i => parseFloat(line.substr(i * 10, (i + 1) * 10)))
        );
    };
    const arraySum = (a) => {
        const s = 0.;
        return a.reduce((prev, curr) => prev + curr, s);
    };
    const sqDist = (xyz1, xyz2) => arraySum(xyz1.map((c, i) => (c - xyz2[i]) * (c - xyz2[i])))

    const template_ref = RDKitModule.get_mol(template_molblock);
    const mol = RDKitModule.get_mol(mol_molblock);
    const molCoords = getCoordArray(mol);
    const bondLength11_12 = Math.sqrt(sqDist(molCoords[11], molCoords[12]));
    const bondLength5_6 = Math.sqrt(sqDist(molCoords[5], molCoords[6]));
    assert(Math.abs(bondLength11_12 - bondLength5_6) < 1.e-4);
    assert(bondLength11_12 > 2.3);
    [ false, true ].forEach(alignOnly => {
        const mol = RDKitModule.get_mol(mol_molblock);
        const res = JSON.parse(mol.generate_aligned_coords(template_ref, JSON.stringify({ allowRGroups: true, alignOnly })));
        assert([ 11, 10, 7, 8, 9, 6 ].every((idx, i) => res.atoms[i] === idx));
        assert(mol.get_smiles() === "C1CC2CCC1N2C1CNC1N1C2CCC1CC2");
        const molAliCoords = getCoordArray(mol);
        const templateCoords = getCoordArray(template_ref);
        res.atoms.forEach((molIdx, templateIdx) => {
            assert(sqDist(molAliCoords[molIdx], templateCoords[templateIdx]) < 1.e-4);
        });
        const bondLengthAli11_12 = Math.sqrt(sqDist(molAliCoords[11], molAliCoords[12]));
        const bondLengthAli5_6 = Math.sqrt(sqDist(molAliCoords[5], molAliCoords[6]));
        assert(Math.abs(bondLengthAli11_12 - bondLengthAli5_6) < 1.e-4);
        if (alignOnly) {
            assert(bondLengthAli11_12 > 2.3);
        } else {
            assert(bondLengthAli11_12 < 1.6);
        }
    });
}

function test_get_mol_no_kekulize() {
    var molIsValid = true;
    try {
        mol = RDKitModule.get_mol("c");
        molIsValid = mol.is_valid();
    } catch (e) {
        molIsValid = false;
    }
    assert(!molIsValid);
    mol = RDKitModule.get_mol("c", JSON.stringify({kekulize: false}));
    assert(mol.is_valid());
}

function test_get_smarts() {
    var mol = RDKitModule.get_mol(`
  MJ201100

  8  8  0  0  0  0  0  0  0  0999 V2000
   -1.0491    1.5839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7635    1.1714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7635    0.3463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0491   -0.0661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3346    0.3463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3346    1.1714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3798    1.5839    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    0.3798   -0.0661    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  2  0
  5  8  1  0  0  2  0
M  RGP  2   7   1   8   2
M  END
`);
    assert(mol.is_valid());
    smarts = mol.get_smarts();
    assert(smarts == "[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1-&!@*)-&!@*");
}

function test_get_cxsmarts() {
    var mol = RDKitModule.get_mol(`
  MJ201100

  8  8  0  0  0  0  0  0  0  0999 V2000
   -1.0491    1.5839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7635    1.1714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7635    0.3463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0491   -0.0661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3346    0.3463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3346    1.1714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3798    1.5839    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    0.3798   -0.0661    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  2  0
  5  8  1  0  0  2  0
M  RGP  2   7   1   8   2
M  END
`);
    assert(mol.is_valid());
    cxsmarts = mol.get_cxsmarts();
    assert(cxsmarts == "[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1-&!@*)-&!@* |" +
        "(-1.0491,1.5839,;-1.7635,1.1714,;-1.7635,0.3463,;-1.0491,-0.0661,;" +
        "-0.3346,0.3463,;-0.3346,1.1714,;0.3798,1.5839,;0.3798,-0.0661,)," +
        "atomProp:6.dummyLabel.R1:7.dummyLabel.R2|");
}

function test_normalize_depiction() {
    var mol = RDKitModule.get_mol(`
  MJ201100

  9 10  0  0  0  0  0  0  0  0999 V2000
   -1.5402    1.3161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2546    0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2546    0.0785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5402   -0.3339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8257    0.0785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8257    0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0411   -0.1764    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4438    0.4909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0410    1.1583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  8  9  2  0  0  0  0
  6  9  1  0  0  0  0
  7  8  1  0  0  0  0
  5  7  1  0  0  0  0
M  END
`);
    const scalingFactor = mol.normalize_depiction();
    assert(scalingFactor > 1.8 && scalingFactor < 1.9);
}

function test_straighten_depiction() {
    var mol1 = RDKitModule.get_mol(`
  MJ201900

  2  1  0  0  0  0  0  0  0  0999 V2000
   -0.3904    2.1535    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1049    1.7410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
M  END
`);
    var mol2 = RDKitModule.get_mol(`
  MJ201900

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.1899    1.9526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5245    1.5401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
M  END
`);
    mol1.normalize_depiction();
    mol1Copy1 = RDKitModule.get_mol_copy(mol1)
    mol1Copy2 = RDKitModule.get_mol_copy(mol1)
    mol1.straighten_depiction();
    mol2.normalize_depiction();
    mol2.straighten_depiction();
    assert(mol1.get_molblock() === mol2.get_molblock());
    mol1Copy1.straighten_depiction(true);
    assert(mol1Copy1.get_molblock() !== mol2.get_molblock());
    assert(mol1Copy1.get_molblock() === mol1Copy2.get_molblock());
}

function test_has_coords() {
    var mol = RDKitModule.get_mol('CC');
    assert(!mol.has_coords());
    var mol2 = RDKitModule.get_mol(mol.get_new_coords());
    assert(mol2.has_coords());
    assert(!mol.has_coords());
    mol.set_new_coords();
    assert(mol.has_coords());
}

function test_kekulize() {
    const badAromaticSmiles = 'c1cccc1';
    var mol = null;
    try {
        mol = RDKitModule.get_mol(badAromaticSmiles);
        assert(!mol.is_valid());
    } catch {
        assert(mol === null);
    }
    mol = RDKitModule.get_mol(badAromaticSmiles, JSON.stringify({ kekulize: false }));
    assert(mol.is_valid());
}

function test_sanitize() {
    const badValenceSmiles = 'N(C)(C)(C)C';
    var mol = null;
    try {
        mol = RDKitModule.get_mol(badValenceSmiles);
        assert(!mol.is_valid());
    } catch {
        assert(mol === null);
    }
    try {
        mol = RDKitModule.get_mol(badValenceSmiles, JSON.stringify({ kekulize: false }));
        assert(!mol.is_valid());
    } catch {
        assert(mol === null);
    }
    mol = RDKitModule.get_mol(badValenceSmiles, JSON.stringify({ sanitize: false }));
    assert(mol.is_valid());
}

function test_flexicanvas() {
    var mol = RDKitModule.get_mol("CCCC");
    assert.equal(mol.is_valid(),1);

    var svg = mol.get_svg(-1,-1);
    assert(svg.search("svg")>0);
    assert(svg.search("width='95px'")>0);
    assert(svg.search("height='21px'")>0);
}

function test_rxn_drawing() {
    {
        var rxn = RDKitModule.get_rxn("[CH3:1][OH:2]>>[CH2:1]=[OH0:2]");
        var svg = rxn.get_svg();
        assert(svg.search("svg") > 0);
        var rxn_from_block = RDKitModule.get_rxn(`$RXN

      RDKit

  1  1
$MOL

     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0
  1  2  6  0
V    1 [C&H3:1]
V    2 [O&H1:2]
M  END
$MOL

     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0
  1  2  2  0
V    1 [C&H2:1]
V    2 [O&H0:2]
M  END`);
        var svg_from_block = rxn_from_block.get_svg();
        assert(svg_from_block.search("svg") > 0);
        assert(svg.match(/<path/g).length > 0);
        assert(svg.match(/<path/g).length === svg_from_block.match(/<path/g).length);
    }
    {
        var rxn = RDKitModule.get_rxn("[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])O.[N-:13]=[N+:14]=[N-:15]>C(Cl)Cl.C(=O)(C(=O)Cl)Cl>[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])[N:13]=[N+:14]=[N-:15]",
                                        JSON.stringify({ useSmiles: true }));
        {
            var svg = rxn.get_svg();
            assert(svg.search("svg") > 0);
            assert(svg.search("stroke-dasharray" < 0));
            assert(svg.search("ellipse") < 0);
        }
        {
            var svg = rxn.get_svg_with_highlights(JSON.stringify({ kekulize: false }));
            assert(svg.search("svg") > 0);
            assert(svg.search("stroke-dasharray" > 0));
            assert(svg.search("ellipse") < 0);
        }
        {
            var svg = rxn.get_svg_with_highlights(JSON.stringify({ highlightByReactant: true }));
            assert(svg.search("svg") > 0);
            assert(svg.search("stroke-dasharray" < 0));
            assert(svg.search("ellipse") > 0);
            assert(svg.search("#BCD6ED") < 0);
            assert(svg.search("#A8D08D") < 0);
        }
        {
            var svg = rxn.get_svg_with_highlights(JSON.stringify({
                highlightColorsReactants: [[.741, .843, .933], [.659, .816, .553]]
            }));
            assert(svg.search("svg") > 0);
            assert(svg.search("stroke-dasharray" < 0));
            assert(svg.search("ellipse") < 0);
        }
        {
            var svg = rxn.get_svg_with_highlights(JSON.stringify({
                highlightByReactant: true,
                highlightColorsReactants: [[.741, .843, .933], [.659, .816, .553]]
            }));
            assert(svg.search("svg") > 0);
            assert(svg.search("stroke-dasharray" < 0));
            assert(svg.search("ellipse") > 0);
            assert(svg.search("#BCD6ED") > 0);
            assert(svg.search("#A8D08D") > 0);
        }
        rxn.delete();
    }
    {
        var rxn = RDKitModule.get_rxn("[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])O.[N-:13]=[N+:14]=[N-:15]>C(Cl)Cl.C(=O)(C(=O)Cl)Cl>[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])[N:13]=[N+:14]=[N-:15]");
        var hasThrown = false;
        var isAromatic = false;
        // if the code is built with exception support it will not throw
        // but rather display molecules in their aromatic form rather
        // than kekulized
        try {
            var svg = rxn.get_svg();
            isAromatic = svg.search("stroke-dasharray" > 0);
        } catch {
            hasThrown = true;
        }
        assert(isAromatic || hasThrown);
    }
}

function test_legacy_stereochem() {
    RDKitModule.use_legacy_stereo_perception(true);
    var mol = RDKitModule.get_mol("C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2");
    assert.equal(mol.is_valid(),1);
    assert.equal(mol.get_smiles(),"CC1CCC2(CC1)CC[C@H](C)C(C)C2");

    RDKitModule.use_legacy_stereo_perception(false);
    mol = RDKitModule.get_mol("C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2");
    assert.equal(mol.is_valid(),1);
    //assert.equal(mol.get_smiles(),"CC1CC2(CC[C@H](C)CC2)CC[C@@H]1C");
}

function test_prop() {
    var mol = RDKitModule.get_mol(`
  MJ201900

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.1899    1.9526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5245    1.5401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
M  END
`);
    assert.equal(mol.has_prop("test1"), false);
    assert.equal(mol.set_prop("test1","val"), true);
    assert.equal(mol.has_prop("test1"), true);
    assert.equal(mol.get_prop("test1"),"val");
    assert.equal(mol.set_prop("test2","val"), true);
    props = mol.get_prop_list(false, false);
    assert.equal(props.get(0), "test1");
    assert.equal(props.get(1), "test2");
}

function test_highlights() {
    var mol = RDKitModule.get_mol('c1cc(O)ccc1');
    var svg = mol.get_svg_with_highlights(JSON.stringify({
        highlightAtomColors: {
            0: [1.0, 0.0, 0.0],
            1: [1.0, 0.0, 0.0],
            2: [1.0, 0.0, 0.0],
            3: [0.0, 1.0, 0.0],
            4: [1.0, 0.0, 0.0],
            5: [1.0, 0.0, 0.0],
            6: [1.0, 0.0, 0.0],
        }, highlightBondColors: {
            2: [0.0, 0.7, 0.9],
        }, highlightAtomRadii: {
            0: 0.1,
            1: 0.1,
            2: 0.1,
            3: 0.8,
            4: 0.1,
            5: 0.1,
            6: 0.1,
        }, atoms: [0, 1, 2, 3, 4, 5, 6],  bonds: [2], width: 127
    }));
    assert(!svg.includes('fill:#FF7F7F'));
    assert(svg.includes('fill:#FF0000'));
    assert(svg.includes('fill:#00FF00'));
    assert(svg.includes('fill:#00B2E5'));
    assert(svg.includes("width='127px'"));
    assert(svg.includes('</svg>'));
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
    test_get_aromatic_kekule_form();
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
    test_generate_aligned_coords_accept_failure();
    test_generate_aligned_coords_align_only();
    test_get_mol_no_kekulize();
    test_get_smarts();
    test_get_cxsmarts();
    test_has_coords();
    test_kekulize();
    test_sanitize();
    test_normalize_depiction();
    test_straighten_depiction();
    test_flexicanvas();
    test_rxn_drawing();
    test_legacy_stereochem();
    test_prop();
    test_highlights();
    waitAllTestsFinished().then(() =>
        console.log("Tests finished successfully")
    );
});
