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
// the default path to RDKit_minimal.js can be overridden through
// the RDKIT_MINIMAL_JS variable if needed
const minimalLib = process.env.RDKIT_MINIMAL_JS || '../demo/RDKit_minimal.js';
console.log('Loading ' + minimalLib);
var initRDKitModule = require(minimalLib);
var RDKitModule;
const fs       = require('fs');
const readline = require('readline');

const extractBondCoords = (svg, bondDetail) => {
    const getStartEndCoords = (bond) => {
        const m = bond.match(/^.*\s+d='M\s+([^,]+),([^ ]+)\s+L\s+([^,]+),([^ ]+)'.*$/);
        return [[m[1], m[2]], [m[3], m[4]]];
    };
    const bond = svg.split('\n').filter(line => line.includes(bondDetail));
    assert.equal(bond.length, 1);
    return getStartEndCoords(bond[0]);
}
const angleDegBetweenVectors = (v1, v2) => 180 / Math.PI * Math.acos((v1[0] * v2[0] + v1[1] * v2[1])
    / Math.sqrt((v1[0] * v1[0] + v1[1] * v1[1]) * (v2[0] * v2[0] + v2[1] * v2[1])));

// the goal here isn't to be comprehensive (the RDKit has tests for that),
// just to make sure that the wrappers are working as expected
function test_basics() {
    var bmol = null;
    try {
        bmol = RDKitModule.get_mol("c1ccccc");
        assert(bmol === null);
    } catch {
        assert(bmol === null);
    }
    var mol = RDKitModule.get_mol("c1ccccc1O");
    assert(mol !== null);
    assert.equal(mol.get_smiles(),"Oc1ccccc1");
    if (typeof Object.getPrototypeOf(mol).get_inchi === 'function') {
        assert.equal(mol.get_inchi(),"InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H");
        assert.equal(mol.get_inchi("-FixedH"),"InChI=1/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H");
    }
    if (typeof RDKitModule.get_inchikey_for_inchi === 'function') {
        assert.equal(RDKitModule.get_inchikey_for_inchi(mol.get_inchi("-FixedH")),"ISWSIDIOOBJBQZ-UHFFFAOYNA-N");
        assert.equal(RDKitModule.get_inchikey_for_inchi(mol.get_inchi()),"ISWSIDIOOBJBQZ-UHFFFAOYSA-N");
    }

    var mb = mol.get_molblock();
    assert(mb.search("M  END")>0);
    var mol2 = RDKitModule.get_mol(mb);
    assert(mol2 !== null);
    assert.equal(mol2.get_smiles(),"Oc1ccccc1");

    var mjson = mol.get_json();
    assert(mjson.search("rdkitjson")>0);
    var mol3 = RDKitModule.get_mol(mjson);
    assert(mol3 !== null);
    assert.equal(mol3.get_smiles(),"Oc1ccccc1");

    var descrs = JSON.parse(mol.get_descriptors());
    assert.equal(descrs.NumAromaticRings,1);
    assert.equal(descrs.NumRings,1);
    assert.equal(descrs.amw,94.11299999999999);

    var checkStringBinaryFpIdentity = (stringFp, binaryFp) => {
        assert.equal(binaryFp.length, Math.ceil(stringFp.length / 8));
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

    {
        var fp1 = mol.get_maccs_fp();
        assert.equal(fp1.length, 167);
        assert.equal((fp1.match(/1/g)||[]).length, 10);
        var fp1Uint8Array = mol.get_maccs_fp_as_uint8array();
        checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
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
    assert(qmol !== null);
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
    assert(mol !== null);
    var mb = mol.get_molblock();
    assert.equal(mb.includes("M  SAP   1  1   8   6"), true);
    var qmol = RDKitModule.get_qmol(molblock);
    assert(qmol !== null);
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
    assert(mol !== null);
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
    mol.delete();
    mol = RDKitModule.get_qmol(`
  MJ201100                      

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.6473    0.8696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3617    0.4571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3617   -0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6473   -0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0671   -0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0671    0.4571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  1  4  0  0  0  0
M  END
`);
    molblock = mol.get_aromatic_form();
    assert (molblock.match(aromRegExp).length === 6);
    assert (molblock.match(kekRegExp) === null);
    molblock = mol.get_kekule_form();
    assert (molblock.match(aromRegExp) === null);
    assert (molblock.match(kekRegExp).length === 6);
    mol.delete();
    mol = RDKitModule.get_qmol(`
  MJ201100                      

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.6473    0.8696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3617    0.4571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3617   -0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6473   -0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0671   -0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0671    0.4571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
`);
    molblock = mol.get_molblock(JSON.stringify({ kekulize: false }));
    assert (molblock.match(aromRegExp) === null);
    assert (molblock.match(kekRegExp).length === 6);
    mol.convert_to_aromatic_form();
    molblock = mol.get_molblock(JSON.stringify({ kekulize: false }));
    assert (molblock.match(aromRegExp).length === 6);
    assert (molblock.match(kekRegExp) === null);
    mol.convert_to_kekule_form();
    molblock = mol.get_molblock(JSON.stringify({ kekulize: false }));
    assert (molblock.match(aromRegExp) === null);
    assert (molblock.match(kekRegExp).length === 6);
    mol.delete();
}

function test_sketcher_services() {
    var mol = RDKitModule.get_mol("C[C@](F)(Cl)/C=C/C(F)Br");
    assert(mol !== null);
    var tags = mol.get_stereo_tags();
    assert.equal(tags,'{"CIP_atoms":[[1,"(S)"],[6,"(?)"]],"CIP_bonds":[[4,5,"(E)"]]}');
}

function test_sketcher_services2() {
    var mol = RDKitModule.get_mol("c1ccccc1");
    assert(mol !== null);
    var molb = mol.add_hs();
    assert(molb.search(" H ")>0);
    assert.equal((molb.match(/ H /g) || []).length,6);

    var mol2 = RDKitModule.get_mol(molb);
    assert(mol2 !== null);
    var molb2 = mol2.get_molblock();
    assert(molb2.search(" H ")>0);
    assert.equal((molb2.match(/ H /g) || []).length,6);

    molb2 = mol2.remove_hs();
    assert(molb2.search(" H ")<0);
}

function test_abbreviations() {
    var bmol = RDKitModule.get_mol("C1CC(C(F)(F)F)C1");
    assert(bmol !== null);
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
    var query = RDKitModule.get_qmol('C1CCCCN1');
    var nonExistingQuery = RDKitModule.get_mol('O=C(O)C(c1ccc(cc1)CCN4CCC(c2nc3ccccc3n2CCOCC)CC4)(C)C');
    const numBitOptions = [-1, 0];
    var patternFpArray = [];
    var sslibObjs = [];
    numBitOptions.forEach((numBits, optIdx) => {
        var smiReader = readline.createInterface({
            input: fs.createReadStream(__dirname + '/../../GraphMol/test_data/compounds.smi')
        });
        var sslib = numBits < 0 ? new RDKitModule.SubstructLibrary() : new RDKitModule.SubstructLibrary(numBits);
        assert(sslib);
        sslibObjs.push(sslib);
        assert.equal(sslib.count_matches(query), 0);
        assert.equal(sslib.get_matches(query), JSON.stringify([]));
        assert.equal(sslib.get_matches_as_uint32array(query).length, 0);
        // var t0 = performance.now()
        // console.log('Started adding trusted SMILES');
        var matches = [];
        var expectedMatches = [39, 64, 80, 127, 128, 234, 240];
        var i = 0;
        var trustedSmiArray = []
        smiReader.on('line', (smi) => {
            sslib.add_trusted_smiles(smi);
            var mol = RDKitModule.get_mol(smi);
            var res = JSON.parse(mol.get_substruct_match(query));
            if (res.atoms) {
                matches.push(i);
            }
            ++i;
            mol.delete();
        });
        smiReader.on('close', () => {
            // var t1 = performance.now();
            // console.log('Finished adding trusted SMILES took ' + (t1 - t0) / 1000 + ' seconds');
            for (var i = 0; i < sslib.size(); ++i) {
                trustedSmiArray.push(sslib.get_trusted_smiles(i));
                let excRaised = false;
                try {
                    var fp = sslib.get_pattern_fp_as_uint8array(i);
                    patternFpArray.push(fp);
                } catch (e) {
                    // this is expected to fail when numBits === 0
                    excRaised = true;
                }
                assert((excRaised && numBits === 0) || (!excRaised && numBits !== 0));
            }
            assert.equal(trustedSmiArray.length, sslib.size());
            assert.equal(patternFpArray.length, sslib.size());
            assert.equal(trustedSmiArray.length, patternFpArray.length);
            var sslib2 = new RDKitModule.SubstructLibrary();
            assert(sslib2);
            sslibObjs.push(sslib2);
            for (var i = 0; i < sslib.size(); ++i) {
                sslib2.add_trusted_smiles_and_pattern_fp(trustedSmiArray[i], patternFpArray[i]);
            }
            assert.equal(sslib.size(), sslib2.size());
            {
                assert.equal(sslib.count_matches(query, false), 7);
                var sslibMatches = sslib.get_matches(query);
                assert.equal(sslibMatches, JSON.stringify(expectedMatches));
                assert.equal(sslibMatches, JSON.stringify(matches));
                var sslibMatchesUInt32Array = sslib.get_matches_as_uint32array(query);
                assert.equal(sslibMatchesUInt32Array.length, expectedMatches.length);
                for (var i = 0; i < expectedMatches.length; ++i) {
                    assert.equal(sslibMatchesUInt32Array[i], expectedMatches[i]);
                }
            }
            {
                assert.equal(sslib.count_matches(nonExistingQuery, false), 0);
                var sslibMatches = sslib.get_matches(nonExistingQuery);
                assert.equal(sslibMatches, JSON.stringify([]));
                var sslibMatchesUInt32Array = sslib.get_matches_as_uint32array(nonExistingQuery);
                assert.equal(sslibMatchesUInt32Array.length, 0);
            }
            {
                assert.equal(sslib2.count_matches(query, false), 7);
                var sslib2Matches = sslib2.get_matches(query);
                assert.equal(sslib2Matches, JSON.stringify(expectedMatches));
                assert.equal(sslib2Matches, JSON.stringify(matches));
                var sslib2MatchesUInt32Array = sslib2.get_matches_as_uint32array(query);
                assert.equal(sslib2MatchesUInt32Array.length, expectedMatches.length);
                for (var i = 0; i < expectedMatches.length; ++i) {
                    assert.equal(sslib2MatchesUInt32Array[i], expectedMatches[i]);
                }
            }
            {
                assert.equal(sslib2.count_matches(nonExistingQuery, false), 0);
                var sslib2Matches = sslib2.get_matches(nonExistingQuery);
                assert.equal(sslib2Matches, JSON.stringify([]));
                var sslib2MatchesUInt32Array = sslib2.get_matches_as_uint32array(nonExistingQuery);
                assert.equal(sslib2MatchesUInt32Array.length, 0);
            }
            if (optIdx === numBitOptions.length - 1) {
                done.test_substruct_library = true;
                query.delete();
                nonExistingQuery.delete();
                sslibObjs.forEach((sslib) => {
                    sslib.delete();
                });
            }
        });
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
    assert.equal(mol.generate_aligned_coords(qmol, JSON.stringify({useCoordGen: true})), "");
}

function test_isotope_labels() {
    var mol = RDKitModule.get_mol("[1*]c1cc([2*])c([3*])c[14c]1");
    assert(mol !== null);

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
    var para_smiles = 'c1ccc(-c2ccc(-c3ccccc3)cc2)cc1';
    var biphenyl_smiles = 'c1ccccc1-c1ccccc1';
    var phenyl_smiles = 'c1ccccc1';
    var template_ref = RDKitModule.get_mol(template_molblock);
    var ortho_meta = RDKitModule.get_mol(ortho_meta_smiles);
    var ortho = RDKitModule.get_mol(ortho_smiles);
    var meta = RDKitModule.get_mol(meta_smiles);
    var para = RDKitModule.get_mol(para_smiles);
    var biphenyl = RDKitModule.get_mol(biphenyl_smiles);
    var phenyl = RDKitModule.get_mol(phenyl_smiles);
    assert.equal(JSON.parse(ortho_meta.generate_aligned_coords(
        template_ref, JSON.stringify({ useCoordGen: false, allowRGroups: true}))).atoms.length, 9);
    [ true, false ].forEach(alignOnly => {
        var opts = JSON.stringify({ useCoordGen: false, allowRGroups: true, alignOnly });
        assert.equal(JSON.parse(ortho_meta.generate_aligned_coords(
            template_ref, opts)).atoms.length, 9);
        assert.equal(JSON.parse(ortho.generate_aligned_coords(
            template_ref, opts)).atoms.length, 8);
        assert.equal(JSON.parse(meta.generate_aligned_coords(
            template_ref, opts)).atoms.length, 8);
        assert.equal(JSON.parse(para.generate_aligned_coords(
            template_ref, opts)).atoms.length, 8);
        assert.equal(JSON.parse(biphenyl.generate_aligned_coords(
            template_ref, opts)).atoms.length, 7);
        assert.equal(JSON.parse(phenyl.generate_aligned_coords(
            template_ref, opts)).atoms.length, 6);
    });
    const pyridineRef = RDKitModule.get_mol(`
  MJ201100                      

  8  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -0.4124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144   -0.4124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144    0.4125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  1  7  1  0  0  0  0
  4  8  1  0  0  0  0
M  END
`);
    const referenceSmarts = "[*:3]a1a([*:1])aa([*:2])aa1";
    [ true, false ].forEach(alignOnly => {
        var opts = JSON.stringify({ useCoordGen: false, allowRGroups: true, alignOnly, referenceSmarts });
        assert.equal(JSON.parse(ortho_meta.generate_aligned_coords(
            pyridineRef, opts)).atoms.length, 8);
        assert.equal(JSON.parse(ortho.generate_aligned_coords(
            pyridineRef, opts)).atoms.length, 7);
        assert.equal(JSON.parse(meta.generate_aligned_coords(
            pyridineRef, opts)).atoms.length, 7);
        assert.equal(JSON.parse(para.generate_aligned_coords(
            pyridineRef, opts)).atoms.length, 8);
        assert.equal(JSON.parse(biphenyl.generate_aligned_coords(
            pyridineRef, opts)).atoms.length, 7);
        assert.equal(JSON.parse(phenyl.generate_aligned_coords(
            pyridineRef, opts)).atoms.length, 6);
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
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  1  2  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  END
`;
    var template_ref = RDKitModule.get_mol(template_molblock);
    var mol = RDKitModule.get_mol(mol_molblock);
    var res = mol.generate_aligned_coords(template_ref, JSON.stringify({
        useCoordGen: false,
        allowRGroups: true,
        acceptFailure: false
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
        molIsValid = (mol !== null);
    } catch (e) {
        molIsValid = false;
    }
    assert(!molIsValid);
    mol = RDKitModule.get_mol("c", JSON.stringify({kekulize: false}));
    assert(mol !== null);
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
    assert(mol !== null);
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
    assert(mol !== null);
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
    var benzeneHoriz = RDKitModule.get_mol(`
  MJ201100                      

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.0785    1.6073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9035    1.6073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3160    0.8928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9036    0.1783    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0786    0.1783    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3339    0.8929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
`);
    var benzeneVert = RDKitModule.get_mol(`
  MJ201100                      

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.2234    1.3054    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4910    1.7178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2055    1.3053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2056    0.4803    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4911    0.0678    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2234    0.4804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
`);
    var benzeneHorizCopy = RDKitModule.get_mol_copy(benzeneHoriz);
    var benzeneVertCopy = RDKitModule.get_mol_copy(benzeneVert);
    benzeneHoriz.straighten_depiction();
    benzeneVert.straighten_depiction();
    assert(benzeneHoriz.get_molblock() !== benzeneHorizCopy.get_molblock());
    assert(benzeneVert.get_molblock() === benzeneVertCopy.get_molblock());
    benzeneHoriz = benzeneHorizCopy;
    benzeneVert = benzeneVertCopy;
    benzeneHorizCopy = RDKitModule.get_mol_copy(benzeneHoriz);
    benzeneVertCopy = RDKitModule.get_mol_copy(benzeneVert);
    benzeneHoriz.straighten_depiction(true);
    benzeneVert.straighten_depiction(true);
    assert(benzeneHoriz.get_molblock() === benzeneHorizCopy.get_molblock());
    assert(benzeneVert.get_molblock() === benzeneVertCopy.get_molblock());
}

function test_has_coords() {
    var mol = RDKitModule.get_mol('CC');
    assert(!mol.has_coords());
    var mol2 = RDKitModule.get_mol(mol.get_new_coords());
    assert.equal(mol2.has_coords(), 2);
    assert(!mol.has_coords());
    mol.set_new_coords();
    assert.equal(mol.has_coords(), 2);
    var mol3 = RDKitModule.get_mol(`
     RDKit          3D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.5909   -0.6086    0.0018 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2459    0.8185   -0.0936 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8271   -0.1760    0.1017 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1425   -0.9648    0.9113 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8767   -1.2011   -0.8967 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5849    1.4596    0.7584 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2393    1.3618   -1.0674 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4853   -0.4161   -0.7761 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.3678   -0.2732    1.0608 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  1  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
  2  7  1  0
  3  8  1  0
  3  9  1  0
M  END
`);
    assert.equal(mol3.has_coords(), 3);
}

function test_kekulize() {
    const badAromaticSmiles = 'c1cccc1';
    var mol = null;
    try {
        mol = RDKitModule.get_mol(badAromaticSmiles);
        assert(mol === null);
    } catch {
        assert(mol === null);
    }
    mol = RDKitModule.get_mol(badAromaticSmiles, JSON.stringify({ kekulize: false }));
    assert(mol !== null);
}

function test_sanitize() {
    const badValenceSmiles = 'N(C)(C)(C)C';
    var mol = null;
    try {
        mol = RDKitModule.get_mol(badValenceSmiles);
        assert(mol === null);
    } catch {
        assert(mol === null);
    }
    try {
        mol = RDKitModule.get_mol(badValenceSmiles, JSON.stringify({ kekulize: false }));
        assert(mol === null);
    } catch {
        assert(mol === null);
    }
    mol = RDKitModule.get_mol(badValenceSmiles, JSON.stringify({ sanitize: false }));
    assert(mol !== null);
}

function test_removehs() {
    const badValenceSmiles = 'N1C=CC(=O)c2ccc(N(C)(C)(C)(C)C)cc12';
    mol = RDKitModule.get_mol(badValenceSmiles, JSON.stringify({ sanitize: false, removeHs: false }));
    assert(mol !== null);
}

function test_flexicanvas() {
    var mol = RDKitModule.get_mol("CCCC");
    assert(mol !== null);

    var svg = mol.get_svg(-1,-1);
    assert(svg.search("svg")>0);
    assert(svg.search("width='87px'")>0);
    assert(svg.search("height='19px'")>0);
}

function test_run_reaction() {
    let rxn1;
    let molList;
    try {
        rxn1 = RDKitModule.get_rxn('[#6:1][O:2]>>[#6:1]=[O:2]');
        molList = molListFromSmiArray(['CC(C)O',]);
        let products;
        try {
            products = rxn1.run_reactants(molList, 10000);
            for (let i = 0; i < products.size(); i++) {
                let element;
                try {
                    element = products.get(i);
                    let mol;
                    try {
                        mol = element.next();
                        assert(mol && mol.get_smiles() === "CC(C)=O");
                    } finally {
                        if (mol) {
                            mol.delete();
                        }
                    }
                } finally {
                    if (element) {
                        element.delete();
                    }
                }
            }
        } finally {
            if (products) {
                products.delete();
            }
        }
    } finally {
        if (rxn1) {
            rxn1.delete();
        }
        if (molList) {
            molList.delete();
        }
    }
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
    var origSetting;
    try {
        origSetting = RDKitModule.use_legacy_stereo_perception(true);
        var mol = RDKitModule.get_mol("O[C@@]1(C)C/C(/C1)=C(/C)\\CC");
        assert(mol !== null);
        assert.equal(mol.get_smiles(),"CCC(C)=C1CC(C)(O)C1");

        RDKitModule.use_legacy_stereo_perception(false);
        mol = RDKitModule.get_mol("O[C@@]1(C)C/C(/C1)=C(/C)\\CC");
        assert(mol !== null);
        assert.equal(mol.get_smiles(),"CC/C(C)=C1\\C[C@](C)(O)C1");
    } finally {
        RDKitModule.use_legacy_stereo_perception(origSetting);
    }
}

function test_allow_non_tetrahedral_chirality() {
    var ctab = `
  Mrv2108 09132105183D          

  5  4  0  0  0  0            999 V2000
   -1.2500    1.4518    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
   -1.2500    2.2768    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4250    1.4518    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -2.0750    1.4518    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2500    0.6268    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
`;
    var origSetting;
    try {
        origSetting = RDKitModule.allow_non_tetrahedral_chirality(true);
        var mol = RDKitModule.get_mol(ctab);
        assert(mol !== null);
        assert.equal(mol.get_smiles(), "[F][Pt@SP3]([F])([Cl])[Cl]");
        RDKitModule.allow_non_tetrahedral_chirality(false);
        var mol = RDKitModule.get_mol(ctab);
        assert(mol !== null);
        assert.equal(mol.get_smiles(), "[F][Pt]([F])([Cl])[Cl]");
    } finally {
        RDKitModule.allow_non_tetrahedral_chirality(origSetting);
    }
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
    assert.equal(mol.clear_prop("test3"), false);
    assert.equal(mol.has_prop("test2"), true);
    assert.equal(mol.clear_prop("test2"), true);
    assert.equal(mol.has_prop("test2"), false);
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

function test_add_chiral_hs() {
    var mol = RDKitModule.get_mol(`
  MJ201100                      

 18 21  0  0  1  0  0  0  0  0999 V2000
   -0.8540   -1.4441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3019   -0.8310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5185   -0.9172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8540   -0.1635    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6825    0.6434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1379    0.7296    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5504    1.4441    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4734   -0.0239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2409    0.3885    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6609   -1.2726    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2130   -1.8857    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9580   -2.6703    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1511   -2.8419    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5990   -2.2287    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0201   -1.7143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5720   -2.3275    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3171   -3.1121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5100   -3.2835    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  4  3  1  0  0  0  0
  4  5  1  0  0  0  0
  6  5  1  0  0  0  0
  6  7  1  1  0  0  0
  6  8  1  0  0  0  0
  8  9  1  1  0  0  0
  8  2  1  0  0  0  0
  4  9  1  1  0  0  0
  2  1  1  1  0  0  0
 10 11  1  0  0  0  0
 11 12  2  0  0  0  0
 12 13  1  0  0  0  0
 13 14  2  0  0  0  0
  1 10  2  0  0  0  0
  1 14  1  0  0  0  0
 15 16  2  0  0  0  0
 16 17  1  0  0  0  0
 11 15  1  0  0  0  0
 17 18  2  0  0  0  0
 12 18  1  0  0  0  0
M  END
`);
    var quinoline_scaffold = RDKitModule.get_mol(`
  MJ201100                      

 10 11  0  0  1  0  0  0  0  0999 V2000
   -8.1001    2.8219    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -8.8145    2.4094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.8145    1.5843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1001    1.1718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3856    1.5843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3856    2.4094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6711    1.1718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9566    1.5842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9566    2.4092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6711    2.8218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
  5  7  1  0  0  0  0
 10  6  1  0  0  0  0
  1  2  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
`);
    var svg1 = mol.get_svg_with_highlights(JSON.stringify({width: 350, height: 300}));
    assert(svg1.includes("width='350px'"));
    assert(svg1.includes("height='300px'"));
    assert(svg1.includes("</svg>"));
    assert(svg1.includes("atom-17"));
    assert(svg1.includes("atom-18"));
    assert(svg1.includes("atom-19"));
    var svg2 = mol.get_svg_with_highlights(JSON.stringify({
        width: 350, height: 300, useMolBlockWedging: true, wedgeBonds: false, addChiralHs: false
    }));
    assert(svg2.includes("width='350px'"));
    assert(svg2.includes("height='300px'"));
    assert(svg2.includes("</svg>"));
    assert(svg2.includes("atom-17"));
    assert(!svg2.includes("atom-18"));
    assert(!svg2.includes("atom-19"));
    assert(mol.get_molblock().includes("4  3  1  6"));
    var molblock = mol.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(!molblock.includes("4  3  1  6"));
    assert(molblock.includes("6  7  1  1"));
    // Here we want to test that the original molblock wedging is preserved and inverted
    // as the coordinates are rigid-body rotated
    var molCopy;
    molCopy = RDKitModule.get_mol_copy(mol);
    assert(JSON.parse(molCopy.generate_aligned_coords(quinoline_scaffold, JSON.stringify({ acceptFailure: false, alignOnly: true }))));
    molblock = molCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(molblock.split('\n').some(line => line.match(/^ [1 ]\d [1 ]\d  [12]  6 *$/)));
    assert(!molblock.split('\n').some(line => line.match(/^ [1 ]\d [1 ]\d  [12]  1 *$/)));
    assert(!molblock.includes("4  3  1  6"));
    assert(molblock.includes("6  7  1  6"));
    molCopy.delete();
    // Here we want to test that the original molblock wedging gets cleared
    // and hence wedging is recomputed as the coordinates are re-generated
    molCopy = RDKitModule.get_mol_copy(mol);
    assert(JSON.parse(molCopy.generate_aligned_coords(quinoline_scaffold, JSON.stringify({ acceptFailure: false }))));
    molblock = molCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(molblock.split('\n').some(line => line.match(/^ [1 ]\d [1 ]\d  [12]  6 *$/)));
    assert(molblock.split('\n').some(line => line.match(/^ [1 ]\d [1 ]\d  [12]  1 *$/)));
    molCopy.delete();
    mol.delete();
    quinoline_scaffold.delete();
}

function getWedgedMolAndInvertedWedges() {
    const wedgedMol = RDKitModule.get_mol(`
     RDKit          2D

 29 34  0  0  1  0  0  0  0  0999 V2000
    1.3719    5.1304    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5985    3.7907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9482    3.7907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7216    5.1304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2685    5.1304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8994    3.5835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5597    4.3569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5597    5.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8994    6.6771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2389    5.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5784    6.6771    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2389    4.3569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3719    2.4510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5985    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3719   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9188   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6921    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9188    2.4510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2389    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0124   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2389   -1.5673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6921   -1.5673    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.8996   -5.0201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2391   -4.2467    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5777   -6.5331    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9909   -5.9040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0124   -2.9070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3306   -6.6772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5784   -5.0201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  3  4  1  0
  5  4  1  6
  5  6  1  0
  6  7  1  0
  7  8  1  0
  9  8  1  1
  5  9  1  0
  9 10  1  0
 10 11  1  1
 10 12  1  0
  6 12  1  1
  2 13  1  0
 13 14  2  0
 14 15  1  0
 15 16  2  0
 16 17  1  0
 17 18  2  0
 13 18  1  0
 17 19  1  0
 19 20  1  0
 20 21  1  0
 21 22  1  0
 16 22  1  0
 23 24  1  0
 23 25  1  0
 25 26  1  0
 24 27  1  0
 27 26  1  0
 26 28  1  0
 24 29  1  0
 28 29  1  0
 21 27  1  0
M  END
`);

    const invertedWedges = `  2  1  1  6
  2  3  1  0
  3  4  1  0
  5  4  1  1
  5  6  1  0
  6  7  1  0
  7  8  1  0
  9  8  1  6
  5  9  1  0
  9 10  1  0
 10 11  1  6
 10 12  1  0
  6 12  1  6
  2 13  1  0
 13 14  2  0
 14 15  1  0
 15 16  2  0
 16 17  1  0
 17 18  2  0
 13 18  1  0
 17 19  1  0
 19 20  1  0
 20 21  1  0
 21 22  1  0
 16 22  1  0
 23 24  1  0
 23 25  1  0
 25 26  1  0
 24 27  1  0
 27 26  1  0
 26 28  1  0
 24 29  1  0
 28 29  1  0
 21 27  1  0
`;
    return { wedgedMol, invertedWedges };
}

function test_wedging_all_within_scaffold() {
    const { wedgedMol, invertedWedges } = getWedgedMolAndInvertedWedges();
    const scaffold = RDKitModule.get_mol(`
     RDKit          2D

 13 14  0  0  1  0  0  0  0  0999 V2000
   -1.6549    2.5755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8814    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6653    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4385    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9854    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6161    1.0286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2766    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2766    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6161    4.1222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9558    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2953    4.1222    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    4.9558    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6549   -0.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  2  3  1  0
  3  4  1  0
  5  4  1  1
  5  6  1  0
  6  7  1  6
  7  8  1  0
  9  8  1  6
  5  9  1  0
  9 10  1  0
 10 11  1  6
 10 12  1  0
  6 12  1  0
  2 13  1  6
M  END
`);
    // the "alignOnly" alignment should succeed and preserve molblock wedging
    // (inverted with respect to the original molecule)
    // it should feature a narrow angle between the bridge bonds
    // as the original geometry of the bridge is preserved
    let wedgedMolCopy;
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(JSON.parse(wedgedMolCopy.generate_aligned_coords(scaffold, JSON.stringify({ acceptFailure: false, alignOnly: true }))));
    const mbAlignOnly = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    const svgAlignOnly = wedgedMolCopy.get_svg_with_highlights(JSON.stringify({
        width: 350, height: 300, useMolBlockWedging: true, wedgeBonds: false, addChiralHs: false
    }));
    {
        const [xy23, xy26] = extractBondCoords(svgAlignOnly, 'atom-23 atom-26');
        const [_, xy25] = extractBondCoords(svgAlignOnly, 'atom-26 atom-25');
        const v1 = [xy23[0] - xy26[0], xy23[1] - xy26[1]];
        const v2 = [xy25[0] - xy26[0], xy25[1] - xy26[1]];
        const v1v2Theta = angleDegBetweenVectors(v1, v2);
        assert(v1v2Theta > 10 && v1v2Theta < 15);
    }
    assert(mbAlignOnly.includes(invertedWedges));
    // the "rebuild" alignment should succeed and preserve molblock wedging
    // (inverted with respect to the original molecule)
    // it should feature a much wider angle between the bridge bonds as the
    // bridged system is entirely rebuilt since it is not part of the scaffold
    wedgedMolCopy.delete();
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(JSON.parse(wedgedMolCopy.generate_aligned_coords(scaffold, JSON.stringify({ acceptFailure: false }))));
    const mbRebuild = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    const svgRebuild = wedgedMolCopy.get_svg_with_highlights(JSON.stringify({
        width: 350, height: 300, useMolBlockWedging: true, wedgeBonds: false, addChiralHs: false
    }));
    {
        const [xy23, xy26] = extractBondCoords(svgRebuild, 'atom-23 atom-26');
        const [_, xy25] = extractBondCoords(svgRebuild, 'atom-26 atom-25');
        const v1 = [xy23[0] - xy26[0], xy23[1] - xy26[1]];
        const v2 = [xy25[0] - xy26[0], xy25[1] - xy26[1]];
        const v1v2Theta = angleDegBetweenVectors(v1, v2);
        assert(v1v2Theta > 105 && v1v2Theta < 110);
    }
    assert(mbRebuild.includes(invertedWedges));
    // the "rebuildCoordGen" alignment should succeed and clear original wedging
    // it should feature an even wider angle between the bridge bonds as CoordGen
    // has a template for the bridged system.
    // Additionally, CoordGen also rebuilds the scaffold, therefore original wedging
    // should be cleared
    wedgedMolCopy.delete();
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(JSON.parse(wedgedMolCopy.generate_aligned_coords(scaffold, JSON.stringify({ acceptFailure: false, useCoordGen: true }))));
    const mbRebuildCoordGen = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    const svgRebuildCoordGen = wedgedMolCopy.get_svg_with_highlights(JSON.stringify({
        width: 350, height: 300, useMolBlockWedging: true, wedgeBonds: false, addChiralHs: false
    }));
    {
        const [xy23, xy26] = extractBondCoords(svgRebuildCoordGen, 'atom-23 atom-26');
        const [_, xy25] = extractBondCoords(svgRebuildCoordGen, 'atom-26 atom-25');
        const v1 = [xy23[0] - xy26[0], xy23[1] - xy26[1]];
        const v2 = [xy25[0] - xy26[0], xy25[1] - xy26[1]];
        const v1v2Theta = angleDegBetweenVectors(v1, v2);
        assert(v1v2Theta > 145 && v1v2Theta < 150);
    }
    assert(!mbRebuildCoordGen.includes(invertedWedges));
    wedgedMolCopy.delete();
}

function test_wedging_outside_scaffold() {
    const { wedgedMol, invertedWedges } = getWedgedMolAndInvertedWedges();
    const scaffold = RDKitModule.get_mol(`
     RDKit          2D

  9 10  0  0  1  0  0  0  0  0999 V2000
   -0.8816    0.5663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6651    0.5663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2958   -0.9804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0435   -0.2072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0435    1.3395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2958    2.1129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6355    1.3395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9750    2.1129    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    2.6355   -0.2072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  3  4  1  6
  4  5  1  0
  6  5  1  6
  2  6  1  0
  6  7  1  0
  7  8  1  6
  7  9  1  0
  3  9  1  0
M  END
`);
    let wedgedMolCopy;
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    // the "alignOnly" alignment should succeed and preserve molblock wedging
    // (inverted with respect to the original molecule)
    // it should feature a narrow angle between the bridge bonds
    // as the original geometry of the bridge is preserved
    assert(JSON.parse(wedgedMolCopy.generate_aligned_coords(scaffold, JSON.stringify({ acceptFailure: false, alignOnly: true }))));
    const mbAlignOnly = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    const svgAlignOnly = wedgedMolCopy.get_svg_with_highlights(JSON.stringify({
        width: 350, height: 300, useMolBlockWedging: true, wedgeBonds: false, addChiralHs: false
    }));
    {
        const [xy23, xy26] = extractBondCoords(svgAlignOnly, 'atom-23 atom-26');
        const [_, xy25] = extractBondCoords(svgAlignOnly, 'atom-26 atom-25');
        const v1 = [xy23[0] - xy26[0], xy23[1] - xy26[1]];
        const v2 = [xy25[0] - xy26[0], xy25[1] - xy26[1]];
        const v1v2Theta = angleDegBetweenVectors(v1, v2);
        assert(v1v2Theta > 10 && v1v2Theta < 15);
    }
    assert(mbAlignOnly.includes(invertedWedges));
    // the "rebuild" alignment should succeed and clear original wedging
    // it should feature a much wider angle between the bridge bonds as the
    // bridged system is entirely rebuilt since it is not part of the scaffold
    wedgedMolCopy.delete();
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(JSON.parse(wedgedMolCopy.generate_aligned_coords(scaffold, JSON.stringify({ acceptFailure: false }))));
    const mbRebuild = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    const svgRebuild = wedgedMolCopy.get_svg_with_highlights(JSON.stringify({
        width: 350, height: 300, useMolBlockWedging: true, wedgeBonds: false, addChiralHs: false
    }));
    {
        const [xy23, xy26] = extractBondCoords(svgRebuild, 'atom-23 atom-26');
        const [_, xy25] = extractBondCoords(svgRebuild, 'atom-26 atom-25');
        const v1 = [xy23[0] - xy26[0], xy23[1] - xy26[1]];
        const v2 = [xy25[0] - xy26[0], xy25[1] - xy26[1]];
        const v1v2Theta = angleDegBetweenVectors(v1, v2);
        assert(v1v2Theta > 105 && v1v2Theta < 110);
    }
    assert(!mbRebuild.includes(invertedWedges));
    // the "rebuildCoordGen" alignment should succeed and clear original wedging
    // it should feature an even wider angle between the bridge bonds as CoordGen
    // has a template for the bridged system.
    // Additionally, CoordGen also rebuilds the scaffold, therefore original wedging
    // should be cleared
    wedgedMolCopy.delete();
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(JSON.parse(wedgedMolCopy.generate_aligned_coords(scaffold, JSON.stringify({ acceptFailure: false, useCoordGen: true }))));
    const mbRebuildCoordGen = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    const svgRebuildCoordGen = wedgedMolCopy.get_svg_with_highlights(JSON.stringify({
        width: 350, height: 300, useMolBlockWedging: true, wedgeBonds: false, addChiralHs: false
    }));
    {
        const [xy23, xy26] = extractBondCoords(svgRebuildCoordGen, 'atom-23 atom-26');
        const [_, xy25] = extractBondCoords(svgRebuildCoordGen, 'atom-26 atom-25');
        const v1 = [xy23[0] - xy26[0], xy23[1] - xy26[1]];
        const v2 = [xy25[0] - xy26[0], xy25[1] - xy26[1]];
        const v1v2Theta = angleDegBetweenVectors(v1, v2);
        assert(v1v2Theta > 145 && v1v2Theta < 150);
    }
    assert(!mbRebuildCoordGen.includes(invertedWedges));
    wedgedMolCopy.delete();
}

function test_wedging_if_no_match() {
    const { wedgedMol, invertedWedges } = getWedgedMolAndInvertedWedges();
    const scaffoldNoMatch = RDKitModule.get_mol(`
     RDKit          2D

 13 14  0  0  1  0  0  0  0  0999 V2000
   -1.6549    2.5755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8814    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6653    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4385    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9854    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6161    1.0286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2766    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2766    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6161    4.1222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9558    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2953    4.1222    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    4.9558    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6549   -0.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  2  3  1  0
  3  4  1  0
  5  4  1  1
  5  6  1  0
  6  7  1  6
  7  8  1  0
  9  8  1  6
  5  9  1  0
  9 10  1  0
 10 11  1  6
 10 12  1  0
  6 12  1  0
  2 13  1  6
M  END
`);
    let wedgedMolCopy;
    let mb;
    const origMolBlock = wedgedMol.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    // the "alignOnly" alignment should return "" if acceptFailure is false
    // and preserve the original coordinates
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(wedgedMolCopy.generate_aligned_coords(scaffoldNoMatch, JSON.stringify({ acceptFailure: false, alignOnly: true })) === "");
    mb = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(mb === origMolBlock);
    assert(!mb.includes(invertedWedges));
    wedgedMolCopy.delete();
    // the "alignOnly" alignment should return "{}" if acceptFailure is true
    // and generate new coordinates, hence wedging should be cleared
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(wedgedMolCopy.generate_aligned_coords(scaffoldNoMatch, JSON.stringify({ acceptFailure: true, alignOnly: true })) === "{}");
    mb = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(mb !== origMolBlock);
    assert(!mb.includes(invertedWedges));
    wedgedMolCopy.delete();
    // the "rebuild" alignment should return "" if acceptFailure is false
    // and preserve the original coordinates
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(wedgedMolCopy.generate_aligned_coords(scaffoldNoMatch, JSON.stringify({ acceptFailure: false })) === "");
    mb = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(mb === origMolBlock);
    assert(!mb.includes(invertedWedges));
    wedgedMolCopy.delete();
    // the "rebuild" alignment should return "{}" if acceptFailure is true
    // and generate new coordinates, hence wedging should be cleared
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(wedgedMolCopy.generate_aligned_coords(scaffoldNoMatch, JSON.stringify({ acceptFailure: true })) === "{}");
    mb = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(mb !== origMolBlock);
    assert(!mb.includes(invertedWedges));
    wedgedMolCopy.delete();
    // the "rebuildCoordGen" alignment should return "" if acceptFailure is false
    // and preserve the original coordinates
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(wedgedMolCopy.generate_aligned_coords(scaffoldNoMatch, JSON.stringify({ acceptFailure: false, useCoordGen: true })) === "");
    mb = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(mb === origMolBlock);
    assert(!mb.includes(invertedWedges));
    wedgedMolCopy.delete();
    // the "rebuildCoordGen" alignment should return "{}" if acceptFailure is true
    // and generate new coordinates, hence wedging should be cleared
    wedgedMolCopy = RDKitModule.get_mol_copy(wedgedMol);
    assert(wedgedMolCopy.generate_aligned_coords(scaffoldNoMatch, JSON.stringify({ acceptFailure: true, useCoordGen: true })) === "{}");
    mb = wedgedMolCopy.get_molblock(JSON.stringify({ useMolBlockWedging: true }));
    assert(mb !== origMolBlock);
    assert(!mb.includes(invertedWedges));
    wedgedMolCopy.delete();
}

function test_get_frags() {
    {
        var mol = RDKitModule.get_mol("n1ccccc1.CC(C)C.OCCCN");
        var expectedFragSmiles = ["c1ccncc1", "CC(C)C", "NCCCO"];
        var expectedFragSmilesNonSanitized = ["CN(C)(C)C", "c1ccc1"];
        var expectedMappings = {
            frags: [0,0,0,0,0,0,1,1,1,1,2,2,2,2,2],
            fragsMolAtomMapping: [[0,1,2,3,4,5],[6,7,8,9],[10,11,12,13,14]],
        };
        var { molList, mappings } = mol.get_frags();
        assert.equal(molList.size(), 3);
        assert(JSON.stringify(JSON.parse(mappings)) === JSON.stringify(expectedMappings));
        var i = 0;
        while (!molList.at_end()) {
            var mol = molList.next();
            assert(mol.get_smiles() === expectedFragSmiles[i++]);
            mol.delete();
        }
        assert(!molList.next());
        molList.delete();
    }
    {
        var mol = RDKitModule.get_mol("N(C)(C)(C)C.c1ccc1", JSON.stringify({sanitize: false}));
        var exceptionThrown = false;
        try {
            mol.get_frags();
        } catch (e) {
            exceptionThrown = true;
        }
        assert(exceptionThrown);
        var { molList, mappings } = mol.get_frags(JSON.stringify({sanitizeFrags: false}));
        assert.equal(molList.size(), 2);
        var i = 0;
        while (!molList.at_end()) {
            var mol = molList.next();
            assert(mol.get_smiles() === expectedFragSmilesNonSanitized[i++]);
            mol.delete();
        }
        assert(!molList.next());
        molList.delete();
    }
}

function test_get_mmpa_frags() {
    {
        var mol = RDKitModule.get_mol("CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3");
        var expectedCores = ["O=C1CN=C(C2CCCCC2)c2ccccc2N1CCC([*:1])[*:2]", "CC([*:1])[*:2]", "CC(C[*:2])[*:1]", "CC(CC[*:2])[*:1]",
        "CC(CCN1C(=O)CN=C([*:2])c2ccccc21)[*:1]", "C([*:1])[*:2]", "C(C[*:2])[*:1]", "O=C1CN=C([*:2])c2ccccc2N1CC[*:1]",
        "C([*:1])[*:2]", "O=C1CN=C([*:1])c2ccccc2N1C[*:2]", "O=C1CN=C([*:1])c2ccccc2N1[*:2]"];
        var expectedSidechains = ["C[*:1].C[*:2]", "C[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1CC[*:2]", "C[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1C[*:2]",
        "C[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1[*:2]", "C1CCC([*:2])CC1.C[*:1]", "CC(C)[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1C[*:2]",
        "CC(C)[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1[*:2]", "C1CCC([*:2])CC1.CC(C)[*:1]", "CC(C)C[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1[*:2]",
        "C1CCC([*:1])CC1.CC(C)C[*:2]", "C1CCC([*:1])CC1.CC(C)CC[*:2]"];
        var pairs = mol.get_mmpa_frags(2, 2, 20);
        assert(pairs.cores);
        assert.equal(pairs.cores.size(), 11);
        assert(pairs.sidechains);
        assert.equal(pairs.sidechains.size(), 11);
        var i = 0;
        while (!pairs.cores.at_end()) {
            var m = pairs.cores.next();
            assert(m.get_smiles() === expectedCores[i++]);
            m.delete();
        }
        i = 0;
        while (!pairs.sidechains.at_end()) {
            var m = pairs.sidechains.next();
            assert(m.get_smiles() === expectedSidechains[i++]);
            m.delete();
        }
        assert(!pairs.cores.next());
        assert(!pairs.sidechains.next());
        pairs.cores.delete();
        pairs.sidechains.delete();
        mol.delete();
    }
    {
        var mol = RDKitModule.get_mol("CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3");
        var expectedSidechains = ["CC(CCN1C(=O)CN=C(C2CCCCC2)c2ccccc21)[*:1].C[*:1]",
            "CC(C)[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1CC[*:1]", "CC(C)C[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1C[*:1]",
            "CC(C)CC[*:1].O=C1CN=C(C2CCCCC2)c2ccccc2N1[*:1]", "C1CCC([*:1])CC1.CC(C)CCN1C(=O)CN=C([*:1])c2ccccc21"];

        var pairs = mol.get_mmpa_frags(1, 1, 20);
        assert(pairs.cores);
        assert.equal(pairs.cores.size(), 5);
        assert(pairs.sidechains);
        assert.equal(pairs.sidechains.size(), 5);
        while (!pairs.cores.at_end()) {
            var m = pairs.cores.next();
            assert(m === null);
        }
        var i = 0;
        while (!pairs.sidechains.at_end()) {
            var m = pairs.sidechains.next();
            assert(m.get_smiles() === expectedSidechains[i++]);
            m.delete();
        }
        assert(!pairs.cores.next());
        assert(!pairs.sidechains.next());
        var numCores = pairs.cores.size();
        for (i = 0; i < numCores; ++i) {
            assert(pairs.cores.at(i) === null);
        }
        for (i = 0; i < numCores; ++i) {
            assert(pairs.cores.pop(0) === null);
        }
        assert.equal(pairs.cores.size(), 0);
        assert(pairs.cores.next() === null);
        pairs.cores.delete();
        pairs.sidechains.delete();
        mol.delete();
    }
}

function test_molzip() {
    {
        var mol1 = RDKitModule.get_mol("F/C=C/[*:1]");
        assert(mol1);
        var mol2 = RDKitModule.get_mol("[*:1]F");
        assert(mol2);
        var expectedLinkage = "F/C=C/F";
        var mol = RDKitModule.molzip(mol1, mol2);
        assert(mol);
        assert(mol.get_smiles() === expectedLinkage);
        mol1.delete();
        mol2.delete();
        mol.delete();
    }
    {
        var mol1 = RDKitModule.get_mol("[C@H]([Xe])(F)([V])");
        assert(mol1);
        var mol2 = RDKitModule.get_mol("[Xe]N.[V]I");
        assert(mol2);
        var expectedLinkage = "N[C@@H](F)I";
        var details = JSON.stringify({ Label: 'AtomType', AtomSymbols: ['Xe', 'V'] });
        var mol = RDKitModule.molzip(mol1, mol2, details);
        assert(mol);
        assert(mol.get_smiles() === expectedLinkage);
        mol1.delete();
        mol2.delete();
        mol.delete();
    }
    if (RDKitModule.RGroupDecomposition) {
        const smis = ['C1CN[C@H]1F', 'C1CN[C@]1(O)F', 'C1CN[C@@H]1F', 'C1CN[CH]1F'];
        const core = RDKitModule.get_qmol('C1CNC1[*:1]');
        const params = {
            rgroupLabelling: 'Isotope',
            allowMultipleRGroupsOnUnlabelled: true,
        };
        const rgd = RDKitModule.get_rgd(core, JSON.stringify(params));
        const smisCanon = smis.map((smi, i) => {
            const mol = RDKitModule.get_mol(smi);
            assert(mol);
            try {
                const res = rgd.add(mol);
                assert(res === i);
                return mol.get_smiles();
            } finally {
                mol.delete();
            }
        });
        rgd.process();
        const rows = rgd.get_rgroups_as_rows();
        const expectedRowMapping = [
            { Core: '[1*][C@@]1([2*])CCN1', R1: '[1*]F', R2: '[2*][H]' },
            { Core: '[1*][C@]1([2*])CCN1', R1: '[1*]F', R2: '[2*]O' },
            { Core: '[1*][C@]1([2*])CCN1', R1: '[1*]F', R2: '[2*][H]' },
            { Core: '[1*]C1([2*])CCN1', R1: '[1*]F', R2: '[2*][H]' },
        ];
        const details = JSON.stringify({ Label: 'Isotope' });
        rows.forEach((row, i) => {
            const foundMapping = getFoundRgdRowAsMap(row, true);
            assert(Object.keys(foundMapping).length === Object.keys(expectedRowMapping[i]).length);
            Object.entries(foundMapping).forEach(([rlabel, smi]) => {
                assert(expectedRowMapping[i][rlabel] && expectedRowMapping[i][rlabel] === smi);
            });
            let mol;
            try {
                mol = RDKitModule.molzip(row, details);
                assert(mol);
                mol.remove_hs_in_place();
                assert(mol.get_smiles() === smisCanon[i]);
            } finally {
                if (mol) {
                    mol.delete();
                }
                Object.values(row).forEach((rgroup) => {
                    if (rgroup) {
                        rgroup.delete();
                    }
                });
            }
    });
    }
}

function test_hs_in_place() {
    {
        var mol = RDKitModule.get_mol("CC");
        assert(!mol.has_coords());
        var descNoH = JSON.parse(mol.get_descriptors());
        assert(`${descNoH.chi0v}` === '2');
        assert(`${descNoH.chi1v}` === '1');
        mol.add_hs_in_place();
        assert(!mol.has_coords());
        assert(mol.get_smiles() === '[H]C([H])([H])C([H])([H])[H]');
        var descH = JSON.parse(mol.get_descriptors());
        assert(`${descH.chi0v}` === '1');
        assert(`${descH.chi1v}` === '0.25');
        mol.delete();
    }
    {
        var mol = RDKitModule.get_mol("C([H])([H])([H])C([H])([H])[H]", JSON.stringify({ removeHs: false }));
        assert(mol.get_smiles() === '[H]C([H])([H])C([H])([H])[H]');
        assert(!mol.has_coords());
        var descH = JSON.parse(mol.get_descriptors());
        assert(`${descH.chi0v}` === '1');
        assert(`${descH.chi1v}` === '0.25');
        // '' will use removeAllHs, while '{}' will use removeHs with default parameters
        ['', '{}'].forEach((details) => {
            var molCopy = mol.copy();
            molCopy.remove_hs_in_place(details);
            assert(!molCopy.has_coords());
            assert(molCopy.get_smiles() === 'CC');
            var descNoH = JSON.parse(molCopy.get_descriptors());
            assert(`${descNoH.chi0v}` === '2');
            assert(`${descNoH.chi1v}` === '1');
            molCopy.delete();
        });
        mol.delete();
    }
    {
        var mol = RDKitModule.get_mol(`
  MJ201100                      

  2  1  0  0  0  0  0  0  0  0999 V2000
  -13.9955    5.1116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.2810    5.5241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
`);
        assert(mol.has_coords());
        assert(mol.get_molblock() === `
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
  -13.9955    5.1116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.2810    5.5241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
`);
        var descNoH = JSON.parse(mol.get_descriptors());
        assert(`${descNoH.chi0v}` === '2');
        assert(`${descNoH.chi1v}` === '1');
        mol.add_hs_in_place();
        assert(mol.has_coords());
        assert(mol.get_smiles() === '[H]C([H])([H])C([H])([H])[H]');
        var descH = JSON.parse(mol.get_descriptors());
        assert(`${descH.chi0v}` === '1');
        assert(`${descH.chi1v}` === '0.25');
        assert(mol.get_molblock().includes(`
     RDKit          2D

  8  7  0  0  0  0  0  0  0  0999 V2000
  -13.9955    5.1116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.2810    5.5241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
`));
        assert(mol.get_molblock().includes(`  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
  2  7  1  0
  2  8  1  0
M  END
`));
        mol.delete();
    }
    {
        var mol = RDKitModule.get_mol(`
  MJ201100                      

  8  7  0  0  0  0  0  0  0  0999 V2000
  -13.9955    5.1116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.2810    5.5241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.4080    5.8260    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -14.7100    4.6991    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -13.5830    4.3971    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -12.8685    4.8096    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -12.5665    5.9366    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -13.6935    6.2385    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
M  END
`, JSON.stringify({ removeHs: false }));
        assert(mol.has_coords());
        assert(mol.get_molblock() === `
     RDKit          2D

  8  7  0  0  0  0  0  0  0  0999 V2000
  -13.9955    5.1116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.2810    5.5241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.4080    5.8260    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -14.7100    4.6991    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -13.5830    4.3971    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -12.8685    4.8096    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -12.5665    5.9366    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -13.6935    6.2385    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
  2  7  1  0
  2  8  1  0
M  END
`);
        var descH = JSON.parse(mol.get_descriptors());
        assert(`${descH.chi0v}` === '1');
        assert(`${descH.chi1v}` === '0.25');
        mol.remove_hs_in_place();
        assert(mol.has_coords());
        assert(mol.get_smiles() === 'CC');
        var descNoH = JSON.parse(mol.get_descriptors());
        assert(`${descNoH.chi0v}` === '2');
        assert(`${descNoH.chi1v}` === '1');
        assert(mol.get_molblock() === `
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
  -13.9955    5.1116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.2810    5.5241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
`);
        mol.delete();
    }
}

function test_query_colour() {
    var mol;
    try {
        mol = RDKitModule.get_qmol('c1ccc2nc([*:1])nc([*:2])c2c1');
        var svg1 = mol.get_svg_with_highlights(JSON.stringify({width: 350, height: 300}));
        assert(svg1.includes("width='350px'"));
        assert(svg1.includes("height='300px'"));
        assert(svg1.includes("</svg>"));
        assert(svg1.includes("#7F7F7F"));
        var svg2 = mol.get_svg_with_highlights(JSON.stringify({width: 350, height: 300, queryColour: [0.0, 0.0, 0.0]}));
        assert(svg2.includes("width='350px'"));
        assert(svg2.includes("height='300px'"));
        assert(svg2.includes("</svg>"));
        assert(!svg2.includes("#7F7F7F"));
    } finally {
        if (mol) {
            mol.delete();
        }
    }
}

function test_alignment_r_groups_aromatic_ring() {
    var mol;
    var scaffold;
    try {
        mol = RDKitModule.get_mol('c1ccc2nccnc2c1');
        assert(mol !== null);
        scaffold = RDKitModule.get_mol(`
  MJ201100                      

  8  8  0  0  0  0  0  0  0  0999 V2000
   -1.0263   -0.3133    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -2.4553    0.5116    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -1.7408   -0.7258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7408   -1.5509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4553   -1.9633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1698   -1.5509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1698   -0.7258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4553   -0.3133    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  3  1  1  0  0  0  0
  8  2  1  0  0  0  0
  4  3  2  0  0  0  0
  5  4  1  0  0  0  0
  6  5  2  0  0  0  0
  7  6  1  0  0  0  0
  8  3  1  0  0  0  0
  8  7  2  0  0  0  0
M  RGP  2   1   2   2   1
M  END`);
        assert(scaffold !== null);
        var res = mol.generate_aligned_coords(scaffold, JSON.stringify({useCoordGen: true, allowRGroups: true}));
        assert(res);
        assert.equal(JSON.parse(res).atoms.length, 8);
        assert.equal(JSON.parse(res).bonds.length, 8);
    } finally {
        if (mol) {
            mol.delete();
        }
    }
    try {
        mol = RDKitModule.get_mol(`
  MJ201100                      

 10 11  0  0  0  0  0  0  0  0999 V2000
    3.6937    2.5671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8687    2.5671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4561    1.8526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8687    1.1382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6937    1.1381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1062    1.8526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9313    1.8527    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.3438    2.5671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9313    3.2816    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1062    3.2816    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  5  6  1  0  0  0  0
  4  5  2  0  0  0  0
  3  4  1  0  0  0  0
  2  3  2  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10  1  1  0  0  0  0
  7  8  2  0  0  0  0
  6  7  1  0  0  0  0
M  END`);
        var res = mol.generate_aligned_coords(scaffold, JSON.stringify({allowRGroups: true, alignOnly: true}));
        assert(res);
        assert.equal(JSON.parse(res).atoms.length, 8);
        assert.equal(JSON.parse(res).bonds.length, 8);
    } finally {
        if (mol) {
            mol.delete();
        }
        if (scaffold) {
            scaffold.delete();
        }
    }
}

function test_is_valid_deprecated() {
    var mol = RDKitModule.get_mol('C');
    assert(mol !== null);
    assert(mol.is_valid());
    var mol;
    try {
        mol = RDKitModule.get_mol('CN(C)(C)C');
    } catch (e) {
        // in case MinimalLib was built without exception support
        mol = null;
    }
    assert(mol === null);
    if (RDKitModule.get_rxn)  {
        var rxn;
        rxn = RDKitModule.get_rxn('C>>N');
        assert(rxn !== null);
        try {
            rxn = RDKitModule.get_rxn('Z>>C');
        } catch (e) {
            // in case MinimalLib was built without exception support
            rxn = null;
        }
        assert(rxn === null);
    }
}

function molListFromSmiArray(smiArray) {
    const molList = new RDKitModule.MolList();
    assert(molList);
    smiArray.forEach((smiName) => {
        const [smi, name] = smiName.split(' ');
        let mol;
        try {
            mol = RDKitModule.get_mol(smi);
            assert(mol);
            molList.append(mol);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
    });
    return molList;
}

function test_mol_list() {
    const smiArray = [ 'C1CC1', 'C1CCCC1' ];
    let molList;
    let mol;
    try {
        molList = molListFromSmiArray(smiArray);
        assert(molList);
        assert.equal(molList.size(), 2);
        assert(!molList.at_end());
        try {
            mol = molList.next();
            assert(mol);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        try {
            mol = molList.next();
            assert(mol);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        mol = molList.next();
        assert(!mol);
        assert(molList.at_end());
        try {
            mol = molList.at(0);
            assert.equal(mol.get_num_atoms(), 3);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        try {
            mol = molList.at(1);
            assert.equal(mol.get_num_atoms(), 5);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        assert(molList.at_end());
        try {
            mol = RDKitModule.get_mol('C1CCC1');
            assert(mol);
            molList.insert(1, mol);
            assert.equal(molList.size(), 3);
            assert(!molList.at_end());
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        molList.reset();
        assert(!molList.at_end());
        try {
            mol = molList.at(1);
            assert.equal(mol.get_num_atoms(), 4);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        try {
            mol = molList.pop(0);
            assert.equal(mol.get_num_atoms(), 3);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        assert.equal(molList.size(), 2);
        // Modifications to molecules in the molList should persist
        // as the underlying C++ object is a shared_ptr
        [0, 1].forEach((loopIdx) => {
            let i = 0;
            molList.reset();
            while (!molList.at_end()) {
                try {
                    mol = molList.next();
                    assert(mol);
                    if (loopIdx == 0) {
                        assert(!mol.has_prop('molIdx'));
                        mol.set_prop('molIdx', `${++i}`);
                    } else {
                        assert(mol.has_prop('molIdx'));
                        i = parseInt(mol.get_prop('molIdx'));
                    }
                } finally {
                    if (mol) {
                        mol.delete();
                    }
                }
            }
            assert.equal(i, 2);
            assert(molList.at_end());
        });
        try {
            mol = molList.pop(0);
            assert.equal(mol.get_num_atoms(), 4);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        assert.equal(molList.size(), 1);
        assert(molList.at_end());
        try {
            mol = molList.pop(0);
            assert.equal(mol.get_num_atoms(), 5);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        assert.equal(molList.size(), 0);
        assert(molList.at_end());
        assert(!molList.pop(0));
        try {
            mol = RDKitModule.get_mol('C1CCCCC1');
            assert(mol);
            molList.append(mol);
            assert.equal(molList.size(), 1);
            assert(!molList.at_end());
        } finally {
            if (mol) {
                mol.delete();
            }
        }
        assert(!molList.at(1));
        try {
            mol = molList.at(0);
            assert(mol);
            assert(mol.get_num_atoms(), 6);
        } finally {
            if (mol) {
                mol.delete();
            }
        }
    } finally {
        if (molList) {
            molList.delete();
        }
    }
}

function test_mcs() {
    {
        const smiArray = [
            "COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)oc2cc1  CHEMBL1479679",
            "COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1333382",
            "Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)nc2cc1  CHEMBL1437584",
            "COc1c(NC(=O)CSc2ccc(Cl)cc2)cc(-c2nc3ccccc3o2)cc1  CHEMBL1601350",
            "Cc1cc2nc(-c3cccc(NC(=O)CSc4ccc(Cl)cc4)c3)oc2cc1C  CHEMBL1398008",
            "Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)nc2cc1  CHEMBL1612903",
            "COc1cc2nc(-c3cc(NC(=O)Cc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1316483",
            "Cc1c(NC(=O)CSc2ccc(Cl)cc2)cccc1-c1nc2cc(Cl)ccc2o1  CHEMBL1568754",
            "COc1ccc2oc(-c3ccc(C)c(NC(=O)COc4cc(C)cc(C)c4)c3)nc2c1  CHEMBL1436972",
            "Cc1ccc(SCC(=O)Nc2cc(-c3nc4cc(C)ccc4o3)c(O)cc2)cc1  CHEMBL1611932",
        ];
        let molList;
        let mcsSmarts;
        try {
            molList = molListFromSmiArray(smiArray);
            let mcsMol;
            try {
                mcsMol = RDKitModule.get_mcs_as_mol(molList);
                assert(mcsMol);
                assert.equal(mcsMol.get_smarts(), '[#6]1:[#6]:[#6]2:[#8]:[#6](:[#7]:[#6]:2:[#6]:[#6]:1)-[#6]1:[#6]:[#6](-[#7]-[#6](=[#8])-[#6]):[#6]:[#6]:[#6]:1');
                mcsSmarts = RDKitModule.get_mcs_as_smarts(molList);
            } finally {
                if (mcsMol) {
                    mcsMol.delete();
                }
            }
        } finally {
            if (molList) {
                molList.delete();
            }
        }
        assert(mcsSmarts);
        assert.equal(mcsSmarts, '[#6]1:[#6]:[#6]2:[#8]:[#6](:[#7]:[#6]:2:[#6]:[#6]:1)-[#6]1:[#6]:[#6](-[#7]-[#6](=[#8])-[#6]):[#6]:[#6]:[#6]:1');
    }
    {
        const smiArray = ["C1CC1N2CC2", "C1CC1N"];
        let molList;
        try {
            molList = molListFromSmiArray(smiArray);
            let mcsMol;
            try {
                mcsMol = RDKitModule.get_mcs_as_mol(molList);
                assert(mcsMol);
                assert.equal(mcsMol.get_num_atoms(), 4);
                assert.equal(mcsMol.get_num_bonds(), 4);
            } finally {
                if (mcsMol) {
                    mcsMol.delete();
                }
            }
            try {
                mcsMol = RDKitModule.get_mcs_as_mol(molList, JSON.stringify({ RingMatchesRingOnly: true }));
                assert(mcsMol);
                assert.equal(mcsMol.get_num_atoms(), 3);
                assert.equal(mcsMol.get_num_bonds(), 3);
            } finally {
                if (mcsMol) {
                    mcsMol.delete();
                }
            }
        } finally {
            if (molList) {
                molList.delete();
            }
        }
    }
    {
        const smiArray = ["NC1CCCCC1", "Cc1ccccc1", "Cc1cnccc1", "CC1CCCCN1"];
        let molList;
        const res = new Set();
        try {
            molList = molListFromSmiArray(smiArray);
            ["Elements", "Any"].forEach((AtomCompare) => {
                ["Order", "OrderExact"].forEach((BondCompare) => {
                    const details = {
                        AtomCompare,
                        BondCompare
                    };
                    let mcsSmarts;
                    mcsSmarts = RDKitModule.get_mcs_as_smarts(molList, JSON.stringify(details));
                    assert(mcsSmarts);
                    res.add(mcsSmarts);
                });
            });
        } finally {
            if (molList) {
                molList.delete();
            }
        }
        assert.equal(res.size, 4);
    }
    {
        const smiArray = [
            "c1cccc(c12)ccc(c2)-c3n(CCC[NH3+])c(nn3)SCCc4c[nH]c(c45)cccc5",
            "c1cccc(c12)sc(c2)-c3n(CCCC[NH3+])c(nn3)SCCc4c[nH]c(c45)cccc5",
        ];
        let molList;
        let mcs;
        try {
            molList = molListFromSmiArray(smiArray);
            mcs = RDKitModule.get_mcs_as_json(molList, JSON.stringify({
                RingMatchesRingOnly: true,
                CompleteRingsOnly: true,
                BondCompare: 'Any',
                AtomCompare: 'Any',
                Timeout: 1,
            }));
        } finally {
            if (molList) {
                molList.delete();
            }
        }
        assert(mcs);
        mcs = JSON.parse(mcs);
        assert(mcs.canceled);
    }
    {
        const smiArray = [
            "Nc1ccc(cc1)C-Cc1c(N)cccc1",
            "Nc1ccc(cc1)C=Cc1c(N)cccc1",
        ];
        let molList;
        let mcs;
        try {
            molList = molListFromSmiArray(smiArray);
            mcs = RDKitModule.get_mcs_as_json(molList);
        } finally {
            if (molList) {
                molList.delete();
            }
        }
        assert(mcs);
        mcs = JSON.parse(mcs);
        assert(!mcs.canceled);
        assert(!Array.isArray(mcs.smarts));
        assert.equal(mcs.numAtoms, 8);
        assert.equal(mcs.numBonds, 8);
        assert(mcs.smarts === '[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]');
        try {
            molList = molListFromSmiArray(smiArray);
            mcs = RDKitModule.get_mcs_as_json(molList, JSON.stringify({StoreAll: true}));
        } finally {
            if (molList) {
                molList.delete();
            }
        }
        assert(mcs);
        mcs = JSON.parse(mcs);
        assert(!mcs.canceled);
        assert(Array.isArray(mcs.smarts));
        assert.equal(mcs.numAtoms, 8);
        assert.equal(mcs.numBonds, 8);
        assert.equal(mcs.smarts.length, 2);
        assert(mcs.smarts.includes('[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7]'));
        assert(mcs.smarts.includes('[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]'));
    }
    {
        const smiArray = [
            "C1CC1",
            "c1ccccc1",
        ];
        let molList;
        let mcs;
        try {
            molList = molListFromSmiArray(smiArray);
            mcs = RDKitModule.get_mcs_as_json(molList, JSON.stringify({ CompleteRingsOnly: true }));
        } finally {
            if (molList) {
                molList.delete();
            }
        }
        assert(mcs);
        mcs = JSON.parse(mcs);
        assert(!mcs.canceled);
        assert(!mcs.numAtoms);
        assert(!mcs.numBonds);
        assert(!mcs.smarts);
        try {
            molList = molListFromSmiArray(smiArray);
            mcs = RDKitModule.get_mcs_as_json(molList, JSON.stringify({ CompleteRingsOnly: true, StoreAll: true }));
        } finally {
            if (molList) {
                molList.delete();
            }
        }
        assert(mcs);
        mcs = JSON.parse(mcs);
        assert(!mcs.canceled);
        assert(!mcs.numAtoms);
        assert(!mcs.numBonds);
        assert(!mcs.smarts.length);
    }
}

function test_get_num_atoms_bonds() {
    var mol = RDKitModule.get_mol('CCCC');
    var molH = RDKitModule.get_mol_copy(mol);
    molH.add_hs_in_place();
    assert.equal(mol.get_num_atoms(), molH.get_num_atoms(true));
    assert.equal(mol.get_num_atoms(), 4);
    assert.equal(molH.get_num_atoms(), 14);
    assert.equal(molH.get_num_atoms(false), 14);
    assert.equal(mol.get_num_bonds(), 3);
    assert.equal(molH.get_num_bonds(), 13);
}

function test_sanitize_no_kekulize_no_setaromaticity() {
    var molblock1 = `
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END`;
    var molblock2 = `
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  1  2  0
M  END`;
    var mol1 = RDKitModule.get_mol(molblock1, JSON.stringify({kekulize: false, setAromaticity: false}));
    var mol2 = RDKitModule.get_mol(molblock2, JSON.stringify({kekulize: false, setAromaticity: false}));
    assert.notEqual(mol1.get_molblock(), mol2.get_molblock());
    mol1.delete();
    mol2.delete();
    mol1 = RDKitModule.get_mol(molblock1);
    mol2 = RDKitModule.get_mol(molblock2);
    assert.equal(mol1.get_molblock(), mol2.get_molblock());
    mol1.delete();
    mol2.delete();
}

function test_partial_sanitization() {
    var mol1 = RDKitModule.get_mol('C1CCC2CCCC2C1', JSON.stringify({
        sanitize: false, removeHs: false, assignStereo: false,
    }));
    var fp1 = mol1.get_morgan_fp(JSON.stringify({radius: 2, nBits: 32}));
    assert.equal(fp1, '00001001010101000000100000110010');
    mol1.delete();
    var mol2 = RDKitModule.get_mol('C1CCC2CCCC2C1', JSON.stringify({
        sanitize: false, removeHs: false, assignStereo: false, fastFindRings: false
    }));
    var exceptionThrown = false;
    try {
        mol2.get_morgan_fp(JSON.stringify({radius: 2, nBits: 32}));
    } catch (e) {
        exceptionThrown = true;
    }
    assert(exceptionThrown);
    mol2.delete();
    let mb;
    mol1 = RDKitModule.get_mol('c1ccccc1N(=O)=O', JSON.stringify({sanitize: false}));
    mb = mol1.get_molblock(JSON.stringify({kekulize: false}));
    assert(mb.includes('  1  2  4  0'));
    assert(mb.includes('  7  8  2  0') && mb.includes('  7  9  2  0'));
    assert(!mb.includes('M  CHG'));
    mol1.delete();
    mol1 = RDKitModule.get_mol('c1ccccc1N(=O)=O', JSON.stringify({sanitize: {SANITIZE_CLEANUP: true}}));
    mb = mol1.get_molblock(JSON.stringify({kekulize: false}));
    assert(mb.includes('  1  2  4  0'));
    assert((mb.includes('  7  8  1  0') && mb.includes('  7  9  2  0'))
        || (mb.includes('  7  8  2  0') && mb.includes('  7  9  1  0')));
    assert(mb.includes('M  CHG  2'));
    mol1.delete();
    mol1 = RDKitModule.get_mol('c1ccccc1N(=O)=O', JSON.stringify({sanitize: {
        SANITIZE_CLEANUP: true, SANITIZE_KEKULIZE: true
    }}));
    mb = mol1.get_molblock(JSON.stringify({kekulize: false}));
    assert(!mb.includes('  1  2  4  0'));
    assert((mb.includes('  7  8  1  0') && mb.includes('  7  9  2  0'))
        || (mb.includes('  7  8  2  0') && mb.includes('  7  9  1  0')));
    assert(mb.includes('M  CHG  2'));
    mol1.delete();
}

function captureStdoutStderr(stdoutCallback, optStderrCallback) {
    if (!stdoutCallback) {
        return null;
    }
    const stderrCallback = optStderrCallback || stdoutCallback;
    const origStdoutWrite = process.stdout.write;
    const origStderrWrite = process.stderr.write;
    process.stdout.write = (chunk) => stdoutCallback(chunk);
    process.stderr.write = (chunk) => stderrCallback(chunk);
    return () => {
        process.stdout.write = origStdoutWrite;
        process.stderr.write = origStderrWrite;
    };
}

function captureStdoutStderr(stdoutCallback, optStderrCallback) {
    if (!stdoutCallback) {
        return null;
    }
    const stderrCallback = optStderrCallback || stdoutCallback;
    const origStdoutWrite = process.stdout.write;
    const origStderrWrite = process.stderr.write;
    process.stdout.write = (chunk) => stdoutCallback(chunk);
    process.stderr.write = (chunk) => stderrCallback(chunk);
    return () => {
        process.stdout.write = origStdoutWrite;
        process.stderr.write = origStderrWrite;
    };
}

function test_capture_logs() {
    const PENTAVALENT_CARBON = 'CC(C)(C)(C)C';
    const PENTAVALENT_CARBON_VALENCE_ERROR = 'Explicit valence for atom # 1 C, 5, is greater than permitted';
    const TETRAVALENT_NITROGEN = 'CN(C)(C)C';
    const TETRAVALENT_NITROGEN_VALENCE_ERROR = 'Explicit valence for atom # 1 N, 4, is greater than permitted';
    RDKitModule.disable_logging();
    assert(!RDKitModule.enable_logging('dummy'));
    assert(RDKitModule.enable_logging('rdApp.info'));
    var restoreStreams;
    var captureHasHappened;
    // Should see no warning on pentavalent carbon below
    captureHasHappened = false;
    restoreStreams = captureStdoutStderr((capturedStreams) => {
        captureHasHappened = true;
    });
    assert(restoreStreams);
    try {
        var mol = RDKitModule.get_mol(PENTAVALENT_CARBON);
        assert(!mol);
    } finally {
        restoreStreams();
    }
    assert(!captureHasHappened);
    RDKitModule.enable_logging("rdApp.error");
    // Should see warning on pentavalent carbon below
    restoreStreams = captureStdoutStderr((capturedStreams) => {
        captureHasHappened = true;
        assert(capturedStreams.includes(PENTAVALENT_CARBON_VALENCE_ERROR));
    });
    assert(restoreStreams);
    try {
        var mol = RDKitModule.get_mol(PENTAVALENT_CARBON);
        assert(!mol);
    } finally {
        restoreStreams();
    }
    assert(captureHasHappened);
    captureHasHappened = false;
    RDKitModule.disable_logging();
    // Should again see no warning on pentavalent carbon below
    restoreStreams = captureStdoutStderr((capturedStreams) => {
        captureHasHappened = true;
    });
    assert(restoreStreams);
    try {
        var mol = RDKitModule.get_mol(PENTAVALENT_CARBON);
        assert(!mol);
    } finally {
        restoreStreams();
    }
    assert(!captureHasHappened);
    var logHandle1;
    var logHandle2;
    var logBuffer1;
    var mol;
    ['set_log_tee', 'set_log_capture'].forEach((func, i) => {
        logHandle1 = RDKitModule[func]('rdApp.*');
        assert(logHandle1);
        logBuffer1 = logHandle1.get_buffer();
        assert(!logBuffer1);
        restoreStreams = captureStdoutStderr((capturedStreams) => {
            captureHasHappened = true;
            assert(capturedStreams.includes(TETRAVALENT_NITROGEN_VALENCE_ERROR));
        });
        assert(restoreStreams);
        try {
            mol = RDKitModule.get_mol(TETRAVALENT_NITROGEN);
            assert(!mol);
        } finally {
            restoreStreams();
        }
        assert(!((func === 'set_log_tee') ^ captureHasHappened));
        captureHasHappened = false;
        logBuffer1 = logHandle1.get_buffer();
        assert(logBuffer1);
        assert(logBuffer1.includes(TETRAVALENT_NITROGEN_VALENCE_ERROR));
        logHandle1.clear_buffer();
        assert(!logHandle1.get_buffer());
        logHandle2 = RDKitModule[func]('rdApp.*');
        assert(!logHandle2);
        restoreStreams = captureStdoutStderr((capturedStreams) => {
            captureHasHappened = true;
            assert(capturedStreams.includes(PENTAVALENT_CARBON_VALENCE_ERROR));
        });
        assert(restoreStreams);
        try {
            mol = RDKitModule.get_mol(PENTAVALENT_CARBON);
            assert(!mol);
        } finally {
            restoreStreams();
        }
        assert(!((func === 'set_log_tee') ^ captureHasHappened));
        captureHasHappened = false;
        logBuffer1 = logHandle1.get_buffer();
        assert(logBuffer1);
        assert(logBuffer1.includes(PENTAVALENT_CARBON_VALENCE_ERROR));
        assert(!logBuffer1.includes(TETRAVALENT_NITROGEN_VALENCE_ERROR));
        logHandle1.delete();
        logHandle2 = RDKitModule[func]('rdApp.*');
        assert(logHandle2);
        logHandle2.delete();
    })
    // Should again see no warning on pentavalent carbon below
    restoreStreams = captureStdoutStderr((capturedStreams) => {
        captureHasHappened = true;
    });
    assert(restoreStreams);
    try {
        var mol = RDKitModule.get_mol(PENTAVALENT_CARBON);
        assert(!mol);
    } finally {
        restoreStreams();
    }
    assert(!captureHasHappened);
}

function test_rgroup_match_heavy_hydro_none_charged() {
    var templateRef = RDKitModule.get_mol(`
  MJ201100                      

  7  7  0  0  0  0  0  0  0  0999 V2000
   -0.5804    1.2045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2948    0.7920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2948   -0.0330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5804   -0.4455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1340   -0.0330    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    0.1340    0.7920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8485   -0.4455    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  6  1  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  5  7  1  0  0  0  0
M  RGP  1   7   1
M  END
`);
    assert(templateRef);
    var mol;
    var allowRGroups = true;
    [true, false].forEach((alignOnly) => {
        mol = RDKitModule.get_mol("Cc1ccccc1");
        assert(mol);
        assert.equal(mol.get_num_atoms(), 7);
        assert.equal(JSON.parse(mol.generate_aligned_coords(templateRef, JSON.stringify({ allowRGroups, alignOnly }))).atoms.length, 7);
        mol.delete();
        mol = RDKitModule.get_mol("c1ccccc1");
        assert(mol);
        assert.equal(mol.get_num_atoms(), 6);
        assert.equal(JSON.parse(mol.generate_aligned_coords(templateRef, JSON.stringify({ allowRGroups, alignOnly }))).atoms.length, 6);
        mol.delete();
        mol = RDKitModule.get_mol("[H]c1ccccc1", JSON.stringify({removeHs: false}));
        assert(mol);
        assert.equal(mol.get_num_atoms(), 7);
        assert.equal(JSON.parse(mol.generate_aligned_coords(templateRef, JSON.stringify({ allowRGroups, alignOnly }))).atoms.length, 7);
        mol.delete();
        mol = RDKitModule.get_mol("n1ccccc1");
        assert(mol);
        assert.equal(mol.get_num_atoms(), 6);
        assert.equal(JSON.parse(mol.generate_aligned_coords(templateRef, JSON.stringify({ allowRGroups, alignOnly }))).atoms.length, 6);
        mol.delete();
        mol = RDKitModule.get_mol("C[n+]1ccccc1");
        assert(mol);
        assert.equal(mol.get_num_atoms(), 7);
        assert.equal(JSON.parse(mol.generate_aligned_coords(templateRef, JSON.stringify({ allowRGroups, alignOnly }))).atoms.length, 7);
        mol.delete();
    });
    templateRef.delete();
}

function test_get_sss_json() {
    var fentanylScaffold = RDKitModule.get_mol(`
  MJ201100                      

 25 27  0  0  0  0  0  0  0  0999 V2000
   -0.3910   -1.2720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2160   -1.2721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6285   -1.9866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2161   -2.7010    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3911   -2.7011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0214   -1.9865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0213   -0.5575    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3911    0.1568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0213    0.8713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3911    1.5859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2161    1.5858    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6287    0.8714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2161    0.1568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6286    2.3002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8463   -0.5575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2588    0.1569    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4536    2.3002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8661    1.5858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6910    1.5858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1036    0.8712    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6912    0.1568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8662    0.1567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4536    0.8713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2589   -1.2720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0839   -1.2719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  9 10  1  0  0  0  0
 12 13  1  0  0  0  0
  8  9  1  0  0  0  0
  8 13  1  0  0  0  0
 10 11  1  0  0  0  0
 11 12  1  0  0  0  0
 11 14  1  0  0  0  0
  1  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7 15  1  0  0  0  0
 14 17  1  0  0  0  0
 17 18  1  0  0  0  0
 19 20  1  0  0  0  0
 20 21  2  0  0  0  0
 21 22  1  0  0  0  0
 22 23  2  0  0  0  0
 18 19  2  0  0  0  0
 23 18  1  0  0  0  0
 15 24  1  0  0  0  0
 24 25  1  0  0  0  0
 15 16  2  0  0  0  0
M  END
`);
    var bilastine = RDKitModule.get_mol("O=C(O)C(c1ccc(cc1)CCN4CCC(c2nc3ccccc3n2CCOCC)CC4)(C)C");
    var referenceSmarts = "c1ccccc1CCN1CCCCC1";
    var res = bilastine.generate_aligned_coords(fentanylScaffold, JSON.stringify({useCoordGen: true, referenceSmarts}));
    assert(res);
    var resExpected = "{\"atoms\":[6,5,4,9,8,7,10,11,12,13,14,15,30,31],\"bonds\":[13,30,14,29,12,36,11,10,9,5,4,33,8,6,7]}";
    assert.equal(res, resExpected);
    assert.equal(JSON.parse(res).bonds.length, 15);
    bilastine.delete();
    fentanylScaffold.delete();
}

function test_relabel_mapped_dummies() {
    var core = RDKitModule.get_mol("c1cc([4*:2])c([3*:1])cn1");
    assert.equal(core.get_cxsmiles(), "c1cc([4*:2])c([3*:1])cn1");
    core.delete();
    core = RDKitModule.get_mol("c1cc([4*:2])c([3*:1])cn1", JSON.stringify({mappedDummiesAreRGroups: true}));
    assert.equal(core.get_cxsmiles(), "*c1ccncc1* |atomProp:0.dummyLabel.R2:7.dummyLabel.R1|");
    core.delete();
}

function test_assign_cip_labels() {
    var origSetting;
    const getTextSection = (svg) => (
      svg.split('\n').map((line) => line.replace(/^<text.+>([^<]*)<\/text>$/, '$1')).join('')
    );
    try {
        origSetting = RDKitModule.use_legacy_stereo_perception(true);
        {
            var mol = RDKitModule.get_mol('C/C=C/c1ccccc1[S@@](C)=O');
            var svg = mol.get_svg_with_highlights(JSON.stringify({noFreetype: true, addStereoAnnotation: true}));
            assert(getTextSection(svg).includes("(S)"));
            assert(!getTextSection(svg).includes("(R)"));
            mol.delete();
        }
        RDKitModule.use_legacy_stereo_perception(false);
        {
            var mol = RDKitModule.get_mol('C/C=C/c1ccccc1[S@@](C)=O');
            var svg = mol.get_svg_with_highlights(JSON.stringify({noFreetype: true, addStereoAnnotation: true}));
            assert(!getTextSection(svg).includes("(S)"));
            assert(!getTextSection(svg).includes("(R)"));
            mol.delete();
        }
        {
            var mol = RDKitModule.get_mol('C/C=C/c1ccccc1[S@@](C)=O', JSON.stringify({assignCIPLabels: true}));
            var svg = mol.get_svg_with_highlights(JSON.stringify({noFreetype: true, addStereoAnnotation: true}));
            assert(!getTextSection(svg).includes("(S)"));
            assert(getTextSection(svg).includes("(R)"));
            mol.delete();
        }
    } finally {
        RDKitModule.use_legacy_stereo_perception(origSetting);
    }
}

function test_smiles_smarts_params() {
    {
        const amoxicillinPubChem = 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)[C@@H](C3=CC=C(C=C3)O)N)C(=O)O)C';
        const mol = RDKitModule.get_mol(amoxicillinPubChem);
        {
            const canonicalSmiles = mol.get_smiles();
            assert(canonicalSmiles === 'CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O');
        }
        ['{}', ''].forEach((emptyJson) => {
            const canonicalSmiles = mol.get_smiles(emptyJson);
            assert(canonicalSmiles === 'CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O');
        });
        const nonCanonicalSmiles = mol.get_smiles(JSON.stringify({canonical: false}));
        assert(nonCanonicalSmiles === 'CC1(C)[C@H](C(=O)O)N2[C@H](S1)[C@H](NC(=O)[C@@H](c1ccc(O)cc1)N)C2=O');
        const canonicalSmilesNoStereo = mol.get_smiles(JSON.stringify({doIsomericSmiles: false}));
        assert(canonicalSmilesNoStereo === 'CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O');
        const nonCanonicalSmilesNoStereo = mol.get_smiles(JSON.stringify({doIsomericSmiles: false, canonical: false}));
        assert(nonCanonicalSmilesNoStereo === 'CC1(C)C(C(=O)O)N2C(S1)C(NC(=O)C(c1ccc(O)cc1)N)C2=O');
        mol.delete();
    }
    {
        const bicyclo221heptane = `
     RDKit          2D

  9 10  0  0  1  0  0  0  0  0999 V2000
   -2.8237   -1.3088    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5723   -0.3996    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5723    1.1473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1011    1.6253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3701    1.1474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3701   -0.3995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6217   -1.3087    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1009   -0.8775    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8083    0.3739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  4  3  1  0
  4  5  1  0
  6  5  1  0
  6  7  1  1
  6  8  1  0
  8  9  1  1
  8  2  1  0
  4  9  1  1
M  END
`;
        const mol = RDKitModule.get_mol(bicyclo221heptane);
        {
            const canonicalCXSmiles = mol.get_cxsmiles();
            const [_, canonicalSmiles, wedging] = canonicalCXSmiles.match(/^(\S+) \|\([^\)]+\),([^\|]+)\|$/);
            assert(canonicalSmiles === 'N[C@@H]1C[C@@H]2C[C@H]1[C@@H](O)C2');
            assert(wedging === 'wD:3.9,wU:1.0,5.4,6.7');
        }
        ['{}', ''].forEach((emptyJson) => {
            const canonicalCXSmiles = mol.get_cxsmiles(emptyJson);
            const [_, canonicalSmiles, wedging] = canonicalCXSmiles.match(/^(\S+) \|\([^\)]+\),([^\|]+)\|$/);
            assert(canonicalSmiles === 'N[C@@H]1C[C@@H]2C[C@H]1[C@@H](O)C2');
            assert(wedging === 'wD:3.9,wU:1.0,5.4,6.7');
        });
        {
            const canonicalCXSmiles = mol.get_cxsmiles(JSON.stringify({restoreBondDirOption: 'RestoreBondDirOptionTrue'}));
            const [_, canonicalSmiles, wedging] = canonicalCXSmiles.match(/^(\S+) \|\([^\)]+\),([^\|]+)\|$/);
            assert(canonicalSmiles === 'N[C@@H]1C[C@@H]2C[C@H]1[C@@H](O)C2');
            assert(wedging === 'wU:1.0,3.3,5.4,6.7');
        }
        {
            const nonCanonicalCXSmilesNoStereo = mol.get_cxsmiles(JSON.stringify({doIsomericSmiles: false, canonical: false, CX_ALL_BUT_COORDS: true}));
            assert(nonCanonicalCXSmilesNoStereo === 'OC1CC2CC(N)C1C2');
            const nonCanonicalCXSmilesNoStereoAtomProp = `${nonCanonicalCXSmilesNoStereo} |atomProp:1.atomProp.1&#46;234|`;
            const molWithAtomProp = RDKitModule.get_mol(nonCanonicalCXSmilesNoStereoAtomProp);
            assert(molWithAtomProp);
            const cxSmilesWithAtomProp = molWithAtomProp.get_cxsmiles(JSON.stringify({CX_ALL_BUT_COORDS: true}));
            assert(cxSmilesWithAtomProp === 'NC1CC2CC(O)C1C2 |atomProp:5.atomProp.1&#46;234|');
            molWithAtomProp.delete();
        }
        mol.delete();
    }
    {
        const chiralQuery = RDKitModule.get_qmol('N-[C@H](-C(-O)=O)-C(-C)-C');
        assert(chiralQuery.get_smarts() === 'N-[C@&H1](-C(-O)=O)-C(-C)-C');
        ['', '{}'].forEach((emptyJson) => {
            assert(chiralQuery.get_smarts(emptyJson) === 'N-[C@&H1](-C(-O)=O)-C(-C)-C');
        });
        assert(chiralQuery.get_smarts(JSON.stringify({doIsomericSmiles: false})) === 'N-[C&H1](-C(-O)=O)-C(-C)-C');
    }
    {
        const chiralQuery = RDKitModule.get_qmol('N-[C@H](-C(-O)=O)-C(-C)-C |atomProp:1.atomProp.1&#46;234|');
        assert(chiralQuery.get_cxsmarts() === 'N-[C@&H1](-C(-O)=O)-C(-C)-C |atomProp:1.atomProp.1&#46;234|');
        ['', '{}'].forEach((emptyJson) => {
            assert(chiralQuery.get_cxsmarts(emptyJson) === 'N-[C@&H1](-C(-O)=O)-C(-C)-C |atomProp:1.atomProp.1&#46;234|');
        });
        assert(chiralQuery.get_cxsmarts(JSON.stringify({doIsomericSmiles: false})) === 'N-[C&H1](-C(-O)=O)-C(-C)-C |atomProp:1.atomProp.1&#46;234|');
    }
}

function test_wedged_bond_atropisomer() {
    var atropisomer = `
  Mrv2311 05242408162D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.0006 -1.54 0 0
M  V30 2 N 2.0006 -3.08 0 0
M  V30 3 C 0.6669 -3.85 0 0
M  V30 4 C -0.6668 -3.08 0 0
M  V30 5 C -0.6668 -1.54 0 0
M  V30 6 C -2.0006 -0.77 0 0
M  V30 7 C 0.6669 -0.77 0 0
M  V30 8 C 0.6669 0.77 0 0
M  V30 9 C -0.6668 1.54 0 0
M  V30 10 C -2.0006 0.77 0 0
M  V30 11 C -0.6668 3.08 0 0
M  V30 12 C 0.6669 3.85 0 0
M  V30 13 C 2.0006 3.08 0 0
M  V30 14 C 2.0006 1.54 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 1 7 5 CFG=3
M  V30 7 2 7 1
M  V30 8 1 7 8
M  V30 9 2 8 9
M  V30 10 1 9 10
M  V30 11 1 9 11
M  V30 12 2 11 12
M  V30 13 1 12 13
M  V30 14 2 13 14
M  V30 15 1 8 14
M  V30 END BOND
M  V30 END CTAB
M  END
`;
    var mol = RDKitModule.get_mol(atropisomer);
    assert(mol);
    var svg = mol.get_svg_with_highlights(JSON.stringify({
        useMolBlockWedging: true, noFreetype: true
    }));
    assert(svg.match(/<path class='bond-5 atom-6 atom-4'.*style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1\.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' \/>/g).length > 6);
    assert(!svg.match(/<path class='bond-5 atom-6 atom-4'.*style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2\.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' \/>/g));
    mol.delete();
}

function test_get_molblock_use_molblock_wedging() {
    var mb = `
     RDKit          2D

  9 10  0  0  1  0  0  0  0  0999 V2000
    1.4885   -4.5513    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0405   -3.9382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8610   -4.0244    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1965   -3.2707    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0250   -2.4637    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2045   -2.3775    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7920   -1.6630    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.8690   -3.1311    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5834   -2.7186    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  4  3  1  0
  4  5  1  0
  6  5  1  0
  6  7  1  1
  6  8  1  0
  8  9  1  1
  8  2  1  0
  4  9  1  1
M  END
`;
    var mol = RDKitModule.get_mol(mb);
    assert(mol);
    var molCopy = RDKitModule.get_mol_copy(mol);
    assert(molCopy);
    var mbRDKitWedging = mol.get_molblock();
    assert(mbRDKitWedging !== mb);
    var mbOrigWedging = mol.get_molblock(JSON.stringify({useMolBlockWedging: true}));
    assert(mb === mbOrigWedging);
    var mbRDKitWedgingPostOrig = mol.get_molblock();
    assert(mb !== mbRDKitWedgingPostOrig);
    assert(mbRDKitWedging === mbRDKitWedgingPostOrig);
    mol.delete();
}

function test_assign_chiral_tags_from_mol_parity() {
    let mol;
    const artemisininCTAB = `68827
  -OEChem-03262404452D

 20 23  0     1  0  0  0  0  0999 V2000
    4.3177    0.4203    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7899    1.1100    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4870   -0.3207    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5402    1.3953    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4004   -1.8275    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6664   -0.2988    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    3.7655    0.1351    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    4.7603   -1.3362    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    2.8959   -0.4383    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    5.5674    0.1351    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    3.9042   -1.9296    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5430    1.1100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9657   -1.4776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6389   -1.8668    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    5.1664    1.8919    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    4.1664    1.8919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5237   -1.3465    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6330   -2.8668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3890    2.8668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0  0  0  0
  6  1  1  0  0  0  0
  2 10  1  0  0  0  0
  2 15  1  0  0  0  0
  3 10  1  0  0  0  0
  3 18  1  0  0  0  0
  4 15  1  0  0  0  0
  5 18  2  0  0  0  0
  6  7  1  0  0  0  0
  6  8  1  0  0  0  0
  6 10  1  0  0  0  0
  7  9  1  0  0  0  0
  7 12  1  0  0  0  0
  8 11  1  0  0  0  0
  8 14  1  0  0  0  0
  9 13  1  0  0  0  0
  9 17  1  0  0  0  0
 11 13  1  0  0  0  0
 12 16  1  0  0  0  0
 14 18  1  0  0  0  0
 14 19  1  0  0  0  0
 15 16  1  0  0  0  0
 15 20  1  0  0  0  0
M  END
`;
    try {
        mol = RDKitModule.get_mol(artemisininCTAB);
        assert(mol);
        assert(mol.get_smiles() === 'CC1CCC2C(C)C(=O)OC3OC4(C)CCC1C32OO4');
    } finally {
        if (mol) {
            mol.delete();
        }
    }
    try {
        mol = RDKitModule.get_mol(artemisininCTAB, JSON.stringify({ assignChiralTypesFromMolParity: true }));
        assert(mol);
        assert(mol.get_smiles() === 'C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4');
    } finally {
        if (mol) {
            mol.delete();
        }
    }
}

function test_make_dummies_queries() {
    let mol;
    let query;
    let match;
    try {
        mol = RDKitModule.get_mol('CN');
        assert(mol);
        query = RDKitModule.get_mol('*N');
        assert(query);
        match = JSON.parse(mol.get_substruct_match(query));
        assert(!match.atoms);
        query.delete();
        query = RDKitModule.get_mol('*N', JSON.stringify(({ makeDummiesQueries: true })));
        assert(query);
        match = JSON.parse(mol.get_substruct_match(query));
        assert(match.atoms.toString() === [0, 1].toString());
    } finally {
        if (mol) {
            mol.delete();
        }
        if (query) {
            query.delete();
        }
    }
}

function test_get_mol_copy() {
    let mol;
    let mol_copy1;
    let mol_copy2;
    try {
        mol = RDKitModule.get_mol('c1ccccn1');
        assert(mol);
        mol_copy1 = RDKitModule.get_mol_copy(mol);
        assert(mol_copy1);
        mol_copy2 = mol.copy();
        assert(mol_copy2);
        assert(mol.get_molblock() === mol_copy1.get_molblock());
        assert(mol.get_molblock() === mol_copy2.get_molblock());
    } finally {
        if (mol) {
            mol.delete();
        }
        if (mol_copy1) {
            mol_copy1.delete();
        }
        if (mol_copy2) {
            mol_copy2.delete();
        }
    }
}

function test_multi_highlights() {
    const mol = RDKitModule.get_mol('[H]c1cc2c(-c3ccnc(Nc4ccc(F)c(F)c4)n3)c(-c3cccc(C(F)(F)F)c3)nn2nc1C', JSON.stringify({removeHs: false}));
    const details = '{"width":250,"height":200,"highlightAtomMultipleColors":{"15":[[0.941,0.894,0.259]],"17":[[0,0.62,0.451]],"21":[[0.902,0.624,0]],"22":[[0.902,0.624,0]],"23":[[0.902,0.624,0]],"24":[[0.902,0.624,0]],"25":[[0.902,0.624,0]],"26":[[0.902,0.624,0]],"27":[[0.902,0.624,0]],"28":[[0.902,0.624,0]],"29":[[0.902,0.624,0]],"30":[[0.902,0.624,0]],"35":[[0.337,0.706,0.914]]},"highlightBondMultipleColors":{"14":[[0.941,0.894,0.259]],"16":[[0,0.62,0.451]],"20":[[0.902,0.624,0]],"21":[[0.902,0.624,0]],"22":[[0.902,0.624,0]],"23":[[0.902,0.624,0]],"24":[[0.902,0.624,0]],"25":[[0.902,0.624,0]],"26":[[0.902,0.624,0]],"27":[[0.902,0.624,0]],"28":[[0.902,0.624,0]],"29":[[0.902,0.624,0]],"34":[[0.337,0.706,0.914]],"38":[[0.902,0.624,0]]},"highlightAtomRadii":{"15":0.4,"17":0.4,"21":0.4,"22":0.4,"23":0.4,"24":0.4,"25":0.4,"26":0.4,"27":0.4,"28":0.4,"29":0.4,"30":0.4,"35":0.4},"highlightLineWidthMultipliers":{"14":2,"16":2,"20":2,"21":2,"22":2,"23":2,"24":2,"25":2,"26":2,"27":2,"28":2,"29":2,"34":2,"38":2}}';
    const svgWithDetails = mol.get_svg_with_highlights(details);
    assert(svgWithDetails.includes('ellipse'));
    const COLORS = ['#009E73', '#55B4E9', '#E69F00', '#EFE342'];
    assert(COLORS.every((color) => svgWithDetails.includes(color)));
    const svgWithOutDetails = mol.get_svg_with_highlights('');
    assert(!svgWithOutDetails.includes('ellipse'));
    assert(!COLORS.some((color) => svgWithOutDetails.includes(color)));
    mol.delete();
}

const getFoundRgdRowAsMap = (row, keep) => Object.fromEntries(Object.entries(row).map(([rlabel, mol]) => {
    try {
        assert(mol);
        assert(mol instanceof RDKitModule.Mol);
        const smi = mol.get_smiles();
        return [rlabel, smi];
    } finally {
        if (!keep && mol) {
            mol.delete();
        }
    }
}));

const getExpectedRgdRowAsMap = (row) => Object.fromEntries(row.split(' ').map((rgroup) => {
    const match = rgroup.match(/^([^:]+):(.+)$/);
    assert(match);
    return match.slice(1);
}));

const getExpectedRgdAsCols = (rowArray) => {
    const res = {};
    rowArray.forEach((row, i) => {
        const rgroupToSmiMap = getExpectedRgdRowAsMap(row);
        Object.entries(rgroupToSmiMap).forEach(([rlabel, smi]) => {
            if (!i) {
                res[rlabel] = [];
            }
            const arr = res[rlabel];
            assert(Array.isArray(arr));
            arr.push(smi);
        });
    });
    return res;
};

function test_singlecore_rgd() {
    const ringData3 = ['c1cocc1CCl', 'c1c[nH]cc1CI', 'c1cscc1CF'];
    const expectedRingData3Rgd = ['Core:c1cc([*:1])co1 R1:ClC[*:1]',
                                  'Core:c1cc([*:1])c[nH]1 R1:IC[*:1]',
                                  'Core:c1cc([*:1])cs1 R1:FC[*:1]'];
    const expectedRingData3RgdAsCols = getExpectedRgdAsCols(expectedRingData3Rgd);
    const core = RDKitModule.get_qmol('*1***[*:1]1');
    assert(core);
    try {
        const scoreMethods = ['Match', 'FingerprintVariance'];
        scoreMethods.forEach((scoreMethod) => {
            const params = {
                scoreMethod,
                allowNonTerminalRGroups: true,
            };
            const rgd = RDKitModule.get_rgd(core, JSON.stringify(params));
            ringData3.forEach((ringData, i) => {
                const mol = RDKitModule.get_mol(ringData);
                assert(mol);
                try {
                    const res = rgd.add(mol);
                    assert(res === i);
                } finally {
                    mol.delete();
                }
            });
            rgd.process();
            const rows = rgd.get_rgroups_as_rows();
            assert(Array.isArray(rows));
            assert(rows.length === ringData3.length);
            rows.forEach((row, i) => {
                const expectedRowMapping = getExpectedRgdRowAsMap(expectedRingData3Rgd[i]);
                const foundMapping = getFoundRgdRowAsMap(row);
                assert(Object.keys(foundMapping).length === Object.keys(expectedRowMapping).length);
                Object.entries(foundMapping).forEach(([rlabel, smi]) => {
                    assert(expectedRowMapping[rlabel] && expectedRowMapping[rlabel] === smi);
                });
            })
            const cols = rgd.get_rgroups_as_columns();
            assert(typeof cols === 'object');
            assert(Object.keys(cols).length === Object.keys(expectedRingData3RgdAsCols).length);
            Object.keys(cols).forEach((rlabel) => {
                const expectedRGroupsAsSmiles = expectedRingData3RgdAsCols[rlabel];
                assert(Array.isArray(expectedRGroupsAsSmiles));
                const rgroupsAsMolList = cols[rlabel];
                assert(rgroupsAsMolList);
                assert(rgroupsAsMolList instanceof RDKitModule.MolList);
                try {
                    assert(expectedRGroupsAsSmiles.length === rgroupsAsMolList.size());
                    let i = 0;
                    while (!rgroupsAsMolList.at_end()) {
                        const mol = rgroupsAsMolList.next();
                        assert(mol);
                        try {
                            assert(mol.get_smiles() === expectedRGroupsAsSmiles[i++]);
                        } finally {
                            mol.delete();
                        }
                    }
                } finally {
                    rgroupsAsMolList.delete();
                }
            });
        });
    } finally {
        core.delete();
    }
}

function test_multicore_rgd() {
    const smiArray = [
        'C1CCNC(Cl)CC1', 'C1CC(Cl)NCCC1', 'C1CCNC(I)CC1', 'C1CC(I)NCCC1',
        'C1CCSC(Cl)CC1', 'C1CC(Cl)SCCC1', 'C1CCSC(I)CC1', 'C1CC(I)SCCC1',
        'C1CCOC(Cl)CC1', 'C1CC(Cl)OCCC1', 'C1CCOC(I)CC1', 'C1CC(I)OCCC1'
    ];
    const expectedRgd = [
        'Core:C1CCNC([*:1])CC1 R1:Cl[*:1]', 'Core:C1CCNC([*:1])CC1 R1:Cl[*:1]',
        'Core:C1CCNC([*:1])CC1 R1:I[*:1]',  'Core:C1CCNC([*:1])CC1 R1:I[*:1]',
        'Core:C1CCSC([*:1])CC1 R1:Cl[*:1]', 'Core:C1CCSC([*:1])CC1 R1:Cl[*:1]',
        'Core:C1CCSC([*:1])CC1 R1:I[*:1]',  'Core:C1CCSC([*:1])CC1 R1:I[*:1]',
        'Core:C1CCOC([*:1])CC1 R1:Cl[*:1]', 'Core:C1CCOC([*:1])CC1 R1:Cl[*:1]',
        'Core:C1CCOC([*:1])CC1 R1:I[*:1]',  'Core:C1CCOC([*:1])CC1 R1:I[*:1]'
    ];
    const expectedRgdAsCols = getExpectedRgdAsCols(expectedRgd);

    const cores = molListFromSmiArray(['C1CCNCCC1', 'C1CCOCCC1', 'C1CCSCCC1']);
    try {
        const rgd = RDKitModule.get_rgd(cores);
        smiArray.forEach((smi, i) => {
            const mol = RDKitModule.get_mol(smi);
            assert(mol);
            try {
                assert(rgd.add(mol) === i);
            } finally {
                mol.delete();
            }
        });
        assert(rgd.process());
        const rows = rgd.get_rgroups_as_rows();
        assert(Array.isArray(rows));
        assert(rows.length === smiArray.length);
        rows.forEach((row, i) => {
            const expectedRowMapping = getExpectedRgdRowAsMap(expectedRgd[i]);
            const foundMapping = getFoundRgdRowAsMap(row);
            assert(Object.keys(foundMapping).length === Object.keys(expectedRowMapping).length);
            Object.entries(foundMapping).forEach(([rlabel, smi]) => {
                assert(expectedRowMapping[rlabel] && expectedRowMapping[rlabel] === smi);
            });
        })
        const cols = rgd.get_rgroups_as_columns();
        assert(typeof cols === 'object');
        assert(Object.keys(cols).length === Object.keys(expectedRgdAsCols).length);
        Object.keys(cols).forEach((rlabel) => {
            const expectedRGroupsAsSmiles = expectedRgdAsCols[rlabel];
            assert(Array.isArray(expectedRGroupsAsSmiles));
            const rgroupsAsMolList = cols[rlabel];
            assert(rgroupsAsMolList);
            assert(rgroupsAsMolList instanceof RDKitModule.MolList);
            try {
                assert(expectedRGroupsAsSmiles.length === rgroupsAsMolList.size());
                let i = 0;
                while (!rgroupsAsMolList.at_end()) {
                    const mol = rgroupsAsMolList.next();
                    assert(mol);
                    try {
                        assert(mol.get_smiles() === expectedRGroupsAsSmiles[i++]);
                    } finally {
                        mol.delete();
                    }
                }
            } finally {
                rgroupsAsMolList.delete();
            }
        });
    } finally {
        cores.delete();
    }
}

function test_multi_highlights() {
    const mol = RDKitModule.get_mol('[H]c1cc2c(-c3ccnc(Nc4ccc(F)c(F)c4)n3)c(-c3cccc(C(F)(F)F)c3)nn2nc1C', JSON.stringify({removeHs: false}));
    const details = '{"width":250,"height":200,"highlightAtomMultipleColors":{"15":[[0.941,0.894,0.259]],"17":[[0,0.62,0.451]],"21":[[0.902,0.624,0]],"22":[[0.902,0.624,0]],"23":[[0.902,0.624,0]],"24":[[0.902,0.624,0]],"25":[[0.902,0.624,0]],"26":[[0.902,0.624,0]],"27":[[0.902,0.624,0]],"28":[[0.902,0.624,0]],"29":[[0.902,0.624,0]],"30":[[0.902,0.624,0]],"35":[[0.337,0.706,0.914]]},"highlightBondMultipleColors":{"14":[[0.941,0.894,0.259]],"16":[[0,0.62,0.451]],"20":[[0.902,0.624,0]],"21":[[0.902,0.624,0]],"22":[[0.902,0.624,0]],"23":[[0.902,0.624,0]],"24":[[0.902,0.624,0]],"25":[[0.902,0.624,0]],"26":[[0.902,0.624,0]],"27":[[0.902,0.624,0]],"28":[[0.902,0.624,0]],"29":[[0.902,0.624,0]],"34":[[0.337,0.706,0.914]],"38":[[0.902,0.624,0]]},"highlightAtomRadii":{"15":0.4,"17":0.4,"21":0.4,"22":0.4,"23":0.4,"24":0.4,"25":0.4,"26":0.4,"27":0.4,"28":0.4,"29":0.4,"30":0.4,"35":0.4},"highlightLineWidthMultipliers":{"14":2,"16":2,"20":2,"21":2,"22":2,"23":2,"24":2,"25":2,"26":2,"27":2,"28":2,"29":2,"34":2,"38":2}}';
    const svgWithDetails = mol.get_svg_with_highlights(details);
    assert(svgWithDetails.includes('ellipse'));
    const COLORS = ['#009E73', '#55B4E9', '#E69F00', '#EFE342'];
    assert(COLORS.every((color) => svgWithDetails.includes(color)));
    const svgWithOutDetails = mol.get_svg_with_highlights('');
    assert(!svgWithOutDetails.includes('ellipse'));
    assert(!COLORS.some((color) => svgWithOutDetails.includes(color)));
    mol.delete();
}

function test_bw_palette() {
    const mol = RDKitModule.get_mol('N');
    assert(mol);
    const details = '{"atomColourPalette":"bw"}';
    const svgWithDetails = mol.get_svg_with_highlights(details);
    assert(svgWithDetails.includes('#000000'));
    assert(!svgWithDetails.includes('#0000FF'));
    const svgWithoutDetails = mol.get_svg();
    assert(!svgWithoutDetails.includes('#000000'));
    assert(svgWithoutDetails.includes('#0000FF'));
    mol.delete();
}

function test_custom_palette() {
    const mol = RDKitModule.get_mol('N');
    assert(mol);
    const details = '{"atomColourPalette\":{"7":[1.0,0.0,0.0]}}';
    const svgWithDetails = mol.get_svg_with_highlights(details);
    assert(svgWithDetails.includes('#FF0000'));
    assert(!svgWithDetails.includes('#0000FF'));
    const svgWithoutDetails = mol.get_svg();
    assert(!svgWithoutDetails.includes('#FF0000'));
    assert(svgWithoutDetails.includes('#0000FF'));
    mol.delete();
}

function test_pickle() {
    const mol = RDKitModule.get_mol('c1ccccn1');
    assert(mol);
    mol.set_prop('a', '1');
    assert(mol.get_prop('a') === '1');
    const pklWithProps = mol.get_as_uint8array();
    assert(pklWithProps);
    let molFromPkl = RDKitModule.get_mol_from_uint8array(pklWithProps);
    assert(molFromPkl);
    assert(molFromPkl.has_prop('a'));
    assert(molFromPkl.get_prop('a') === '1');
    molFromPkl.delete();
    const pklWithoutProps = mol.get_as_uint8array(JSON.stringify({propertyFlags: { NoProps: true } }));
    molFromPkl = RDKitModule.get_mol_from_uint8array(pklWithoutProps);
    assert(molFromPkl);
    assert(!molFromPkl.has_prop('a'));
    molFromPkl.delete();
}

function test_remove_hs_details() {
    let mol;
    let smi;
    mol = RDKitModule.get_mol('C([H])([H])([H])C([2H])([H])C([H])([H])[H]', JSON.stringify({removeHs: false}));
    smi = mol.get_smiles();
    assert(smi === '[H]C([H])([H])C([H])([2H])C([H])([H])[H]');
    assert(mol.remove_hs_in_place(JSON.stringify({removeIsotopes: false})));
    smi = mol.get_smiles();
    assert(smi === '[2H]C(C)C');
    assert(mol.remove_hs_in_place(JSON.stringify({removeIsotopes: true})));
    smi = mol.get_smiles();
    assert(smi === 'CCC');
    mol.delete();
}

function test_get_mol_remove_hs() {
    const mbIn = `
  MJ240300                      

  8  8  0  0  0  0  0  0  0  0999 V2000
   -1.4955    1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2099    0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2099   -0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4955   -0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7810   -0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7810    0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0666    1.1152    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9244   -0.5348    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  0  0
  3  8  1  0  0  0  0
M  ISO  1   7   2
M  END
`;
    const noDetails = {};
    const removeHsTrue = { removeHs: true };
    const removeHsFalse = { removeHs: false };
    const deuteriumCoords = '  -0.0666    1.1152    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0';
    const hydrogenCoords = '  -2.9244   -0.5348    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0';
    let mol;
    let mbOut;
    let smiOut;
    let jb;
    let smi;
    [noDetails, removeHsFalse].forEach((details) => {
        mol = RDKitModule.get_mol(mbIn, JSON.stringify(details));
        assert(mol);
        mbOut = mol.get_molblock();
        assert(mbOut);
        assert(mbOut.includes(deuteriumCoords));
        assert(mbOut.includes(hydrogenCoords));
        if (!jb) {
            jb = mol.get_json();
            assert(jb);
        }
        if (!smi) {
            smi = mol.get_smiles();
            assert(smi);
        }
        mol.delete();
    });
    mol = RDKitModule.get_mol(mbIn, JSON.stringify(removeHsTrue));
    assert(mol);
    mbOut = mol.get_molblock();
    assert(mbOut);
    assert(mbOut.includes(deuteriumCoords));
    assert(!mbOut.includes(hydrogenCoords));
    mol.delete();
    [noDetails, removeHsFalse].forEach((details) => {
        mol = RDKitModule.get_mol(mbIn, JSON.stringify(details));
        assert(mol);
        mbOut = mol.get_molblock();
        assert(mbOut);
        assert(mbOut.includes(deuteriumCoords));
        assert(mbOut.includes(hydrogenCoords));
        mol.delete();
    });
    mol = RDKitModule.get_mol(jb, JSON.stringify(removeHsTrue));
    assert(mol);
    mbOut = mol.get_molblock();
    assert(mbOut);
    assert(mbOut.includes(deuteriumCoords));
    assert(!mbOut.includes(hydrogenCoords));
    mol.delete();
    [noDetails, removeHsTrue].forEach((details) => {
        mol = RDKitModule.get_mol(smi, JSON.stringify(details));
        assert(mol);
        smiOut = mol.get_smiles();
        assert(smiOut);
        assert(smiOut.includes('[2H]'));
        assert(!smiOut.includes('[H]'));
        mol.delete();
    });
    mol = RDKitModule.get_mol(smi, JSON.stringify(removeHsFalse));
    assert(mol);
    smiOut = mol.get_smiles();
    assert(smiOut);
    assert(smiOut.includes('[2H]'));
    assert(smiOut.includes('[H]'));
    mol.delete();
}

function test_png_metadata() {
    const PNG_COLCHICINE_NO_METADATA = "/../../GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    const PNG_COLCHICINE_WITH_METADATA = "/../../GraphMol/FileParsers/test_data/colchicine.png";
    const PNG_PENICILLIN_METADATA = "penicillin_metadata.png";
    const PNG_PENICILLIN_AMOXICILLIN_METADATA = "penicillin_amoxicillin_metadata.png";
    const PNG_COLCHICINE_AMOXICILLIN_METADATA = "colchicine_amoxicillin_metadata.png";
    const BENZYLPENICILLIN_SMI = "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)Cc3ccccc3)C(=O)O)C";
    const BENZYLPENICILLIN_CAN_SMI = "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O";
    const AMOXICILLIN_SMI = "O=C(O)[C@@H]2N3C(=O)[C@@H](NC(=O)[C@@H](c1ccc(O)cc1)N)[C@H]3SC2(C)C";
    const AMOXICILLIN_CAN_SMI = "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O";
    let mol;
    let mols;
    let molSan;
    const png_no_metadata_buf = fs.readFileSync(__dirname + PNG_COLCHICINE_NO_METADATA);
    const png_no_metadata_blob = new Uint8Array(png_no_metadata_buf.length);
    const png_no_metadata_blob2 = new Uint8Array(png_no_metadata_buf.length);
    png_no_metadata_buf.copy(png_no_metadata_blob);
    png_no_metadata_buf.copy(png_no_metadata_blob2);
    const png_with_metadata_buf = fs.readFileSync(__dirname + PNG_COLCHICINE_WITH_METADATA);
    const png_with_metadata_blob = new Uint8Array(png_with_metadata_buf.length);
    png_with_metadata_buf.copy(png_with_metadata_blob);
    assert(!RDKitModule.get_mol_from_png_blob(png_no_metadata_blob));
    assert(!RDKitModule.get_mols_from_png_blob(png_no_metadata_blob));
    let penicillin = RDKitModule.get_mol(BENZYLPENICILLIN_SMI);
    assert(penicillin);
    assert(penicillin.set_new_coords());
    let png_penicillin_metadata_blob = penicillin.add_to_png_blob(png_no_metadata_blob,
        "{\"includePkl\":false,\"includeSmiles\":true,\"includeMol\":true}");
    penicillin.delete();
    assert(png_penicillin_metadata_blob);
    fs.writeFileSync(PNG_PENICILLIN_METADATA, png_penicillin_metadata_blob);
    penicillin = RDKitModule.get_mol_from_png_blob(png_penicillin_metadata_blob);
    assert(penicillin);
    assert.equal(penicillin.has_coords(), 2);
    penicillin.delete();
    RDKitModule.enable_logging()
    mol = RDKitModule.get_mol_from_png_blob(png_with_metadata_blob);
    assert(mol);
    assert.equal(mol.has_coords(), 2);
    mol.delete();
    mols = RDKitModule.get_mols_from_png_blob(png_with_metadata_blob);
    assert(!mols);
    mols = RDKitModule.get_mols_from_png_blob(png_with_metadata_blob,
        JSON.stringify({ includePkl: true, includeSmiles: true }));
    assert(!mols);
    mols = RDKitModule.get_mols_from_png_blob(
        png_with_metadata_blob, JSON.stringify({ includeSmiles: true }));
    assert(mols);
    assert.equal(mols.size(), 1);
    mol = mols.next();
    assert.equal(mol.has_coords(), 2);
    mol.delete();
    mols.delete();
    assert(!RDKitModule.get_mol_from_png_blob(png_penicillin_metadata_blob,
        "{\"includePkl\":true,\"includeSmiles\":false,\"includeMol\":false}"));
    penicillin = RDKitModule.get_mol_from_png_blob(png_penicillin_metadata_blob,
        "{\"includePkl\":false,\"includeSmiles\":false,\"includeMol\":true,\"sanitize\":false,\"removeHs\":false,\"assignStereo\":false,\"fastFindRings\":false}");
    assert(penicillin);
    assert.equal(penicillin.has_coords(), 2);
    smi = penicillin.get_smiles();
    assert.notEqual(smi, BENZYLPENICILLIN_CAN_SMI);
    molSan = RDKitModule.get_mol(smi)
    assert.equal(molSan.get_smiles(), BENZYLPENICILLIN_CAN_SMI);
    molSan.delete();
    penicillin.delete();
    assert(!RDKitModule.get_mol_from_png_blob(png_no_metadata_blob2));
    penicillin = RDKitModule.get_mol(BENZYLPENICILLIN_SMI);
    assert(penicillin);
    let amoxicillin = RDKitModule.get_mol(AMOXICILLIN_SMI);
    assert(amoxicillin);
    assert(amoxicillin.set_new_coords());
    png_penicillin_metadata_blob = penicillin.add_to_png_blob(png_no_metadata_blob2,
        "{\"includePkl\":false,\"includeMol\":true,\"CX_ALL_BUT_COORDS\":true}");
    assert(png_penicillin_metadata_blob);
    let png_penicillin_amoxicillin_metadata_blob = amoxicillin.add_to_png_blob(png_penicillin_metadata_blob,
        "{\"includePkl\":false,\"includeMol\":true,\"CX_ALL_BUT_COORDS\":true}");
    assert(png_penicillin_amoxicillin_metadata_blob);
    fs.writeFileSync(PNG_PENICILLIN_AMOXICILLIN_METADATA, png_penicillin_amoxicillin_metadata_blob);
    mol = RDKitModule.get_mol_from_png_blob(png_penicillin_amoxicillin_metadata_blob,
        "{\"sanitize\":false,\"removeHs\":false,\"assignStereo\":false,\"fastFindRings\":false}");
    assert(mol);
    assert(!mol.has_coords());
    smi = mol.get_smiles();
    assert.equal(smi, BENZYLPENICILLIN_CAN_SMI);
    mol.delete();
    assert(!RDKitModule.get_mol_from_png_blob(png_penicillin_amoxicillin_metadata_blob,
        "{\"includePkl\":false,\"includeSmiles\":false,\"includeMol\":false}"));
    mol = RDKitModule.get_mol_from_png_blob(png_penicillin_amoxicillin_metadata_blob,
        "{\"includePkl\":true,\"includeSmiles\":false,\"includeMol\":true}");
    assert(mol);
    assert.equal(mol.has_coords(), 2);
    smi = mol.get_smiles();
    assert.equal(smi, BENZYLPENICILLIN_CAN_SMI);
    mol.delete();
    mols = RDKitModule.get_mols_from_png_blob(png_penicillin_amoxicillin_metadata_blob);
    assert(!mols);
    mols = RDKitModule.get_mols_from_png_blob(
        png_penicillin_amoxicillin_metadata_blob, JSON.stringify({ includeSmiles: true }));
    assert(mols);
    assert.equal(mols.size(), 2);
    mol = mols.at(0);
    assert(!mol.has_coords());
    smi = mol.get_smiles();
    assert.equal(smi, BENZYLPENICILLIN_CAN_SMI);
    assert(mol.get_morgan_fp());
    mol.delete();
    mol = mols.at(1);
    assert(!mol.has_coords());
    smi = mol.get_smiles();
    assert.equal(smi, AMOXICILLIN_CAN_SMI);
    assert(mol.get_morgan_fp());
    mol.delete();
    mols.delete();
    mols = RDKitModule.get_mols_from_png_blob(png_penicillin_amoxicillin_metadata_blob, JSON.stringify({
        includePkl: false, includeMol: true, sanitize: false, removeHs: false, assignStereo: false, fastFindRings: false
    }));
    assert(mols);
    assert.equal(mols.size(), 2);
    mol = mols.at(0);
    assert(mol.has_coords());
    smi = mol.get_smiles();
    assert.notEqual(smi, BENZYLPENICILLIN_CAN_SMI);
    molSan = RDKitModule.get_mol(smi)
    assert.equal(molSan.get_smiles(), BENZYLPENICILLIN_CAN_SMI);
    molSan.delete();
    mol.delete();
    mol = mols.at(1);
    assert(mol.has_coords());
    smi = mol.get_smiles();
    assert.notEqual(smi, AMOXICILLIN_CAN_SMI);
    molSan = RDKitModule.get_mol(smi)
    assert.equal(molSan.get_smiles(), AMOXICILLIN_CAN_SMI);
    molSan.delete();
    mol.delete();
    mols.delete();
    let png_colchicine_amoxicillin_metadata_blob = amoxicillin.add_to_png_blob(
        png_with_metadata_blob, JSON.stringify({ includeMol: true, CX_ALL_BUT_COORDS:true }));
    assert(png_colchicine_amoxicillin_metadata_blob);
    fs.writeFileSync(PNG_COLCHICINE_AMOXICILLIN_METADATA, png_colchicine_amoxicillin_metadata_blob);
    mols = RDKitModule.get_mols_from_png_blob(png_colchicine_amoxicillin_metadata_blob,
        "{\"includePkl\":false,\"includeMol\":true}");
    assert.equal(mols.size(), 1);
    mol = mols.at(0);
    assert.equal(mol.has_coords(), 2);
    mol.delete();
    mols.delete();
    mols = RDKitModule.get_mols_from_png_blob(png_colchicine_amoxicillin_metadata_blob,
        "{\"includeSmiles\":true}");
    assert.equal(mols.size(), 2);
    assert.equal(mols.at(0).has_coords(), 2);
    assert(!mols.at(1).has_coords());
    mols.delete();
    penicillin.delete();
    amoxicillin.delete();

    const colchicine = RDKitModule.get_mol('COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1[C@@H](NC(C)=O)CC2');
    assert(colchicine);
    png_no_metadata_buf.copy(png_no_metadata_blob);
    let png_colchicine_metadata_blob;
    png_colchicine_metadata_blob = colchicine.add_to_png_blob(png_no_metadata_blob);
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert(!mol.has_coords());
    mol.delete();
    // use SMILES
    png_colchicine_metadata_blob = colchicine.add_to_png_blob(png_no_metadata_blob,
        JSON.stringify({ includePkl: false }));
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert(!mol.has_coords());
    mol.delete();
    // use MOL
    png_colchicine_metadata_blob = colchicine.add_to_png_blob(png_no_metadata_blob,
        JSON.stringify({ includePkl: false, includeSmiles: false, includeMol: true }));
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert.equal(mol.has_coords(), 2);
    mol.delete();
    // use PKL
    colchicine.set_new_coords();
    assert.equal(colchicine.has_coords(), 2);
    const PROPERTY_NAME = 'property';
    const PROPERTY_VALUE = 'value';
    colchicine.set_prop(PROPERTY_NAME, PROPERTY_VALUE);
    png_colchicine_metadata_blob = colchicine.add_to_png_blob(png_no_metadata_blob, JSON.stringify({
        includePkl: true, includeSmiles: false, includeMol: false, propertyFlags: { NoProps: true }
    }));
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert.equal(mol.has_coords(), 2);
    assert(!mol.has_prop(PROPERTY_NAME));
    mol.delete();
    png_colchicine_metadata_blob = colchicine.add_to_png_blob(png_no_metadata_blob, JSON.stringify({
        includePkl: true, includeSmiles: false, includeMol: false, propertyFlags: { AllProps: true }
    }));
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert.equal(mol.has_coords(), 2);
    assert(mol.has_prop(PROPERTY_NAME));
    assert(mol.get_prop(PROPERTY_NAME) == PROPERTY_VALUE);
    mol.delete();
    png_colchicine_metadata_blob = colchicine.add_to_png_blob(png_no_metadata_blob, JSON.stringify({
        includePkl: false, includeSmiles: true, includeMol: false,
        propertyFlags: { NoProps: true }, CX_ALL_BUT_COORDS: true
    }));
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert(!mol.has_coords());
    assert(!mol.has_prop(PROPERTY_NAME));
    mol.delete();
    // use original wedging
    const colchicineUnusualWedging = RDKitModule.get_mol(`
     RDKit          2D

 29 31  0  0  0  0  0  0  0  0999 V2000
    6.4602    1.0300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3062    1.9883    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8993    1.4680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7453    2.4262    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3384    1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0856    0.4273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2396   -0.5309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9868   -2.0094    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1408   -2.9677    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6465   -0.0106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8005   -0.9688    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5477   -2.4474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2280   -0.2968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1857   -1.7387    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6836   -2.9611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1813   -3.0436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7569   -4.4288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2442   -4.6230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1797   -1.9240    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6215   -2.3378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9268   -0.4455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6132    0.2787    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0269    1.7205    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5055    1.9733    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0258    3.3802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5043    3.6330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0675    4.5342    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1576    2.9429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3401    3.0254    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  1  0
  8  9  1  0
  7 10  2  0
 10 11  1  0
 11 12  1  0
  6 13  1  0
 13 14  2  0
 14 15  1  0
 15 16  2  0
 16 17  1  0
 17 18  1  0
 16 19  1  0
 19 20  2  0
 19 21  1  0
 21 22  2  0
 23 22  1  1
 23 24  1  0
 24 25  1  0
 25 26  1  0
 25 27  2  0
 23 28  1  0
 28 29  1  0
 10  3  1  0
 22 13  1  0
 29  5  1  0
M  END
`);
    assert(colchicineUnusualWedging);
    assert.equal(colchicineUnusualWedging.has_coords(), 2);
    assert(colchicineUnusualWedging.get_cxsmiles().includes('wU:22.24|'));
    assert(colchicineUnusualWedging.get_cxsmiles(JSON.stringify({ CX_ALL: true, restoreBondDirOption: 'RestoreBondDirOptionTrue' })).includes('wU:22.23|'));
    png_colchicine_metadata_blob = colchicineUnusualWedging.add_to_png_blob(png_no_metadata_blob, JSON.stringify({
        includePkl: true, includeSmiles: true, includeMol: true,
        propertyFlags: { AtomProps: true, BondProps: true }
    }));
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    // the mol is restored from pickle, so it will retain
    // original molblock wedging
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert.equal(mol.has_coords(), 2);
    assert(mol.get_cxsmiles(JSON.stringify({ CX_ALL: true, restoreBondDirOption: 'RestoreBondDirOptionTrue' })).includes('wU:22.23|'));
    assert(mol.get_cxsmiles(JSON.stringify({ CX_ALL: true, restoreBondDirOption: 'RestoreBondDirOptionClear' })).includes('wU:22.24|'));
    assert(mol.get_molblock().includes(' 23 24  1  1'));
    assert(mol.get_molblock(JSON.stringify({ useMolBlockWedging: true })).includes(' 23 22  1  1'));
    mol.delete();
    // the mol is restored from CXSMILES, so it will not retain
    // original molblock wedging, as it was not stored in the CXSMILES string
    png_colchicine_metadata_blob = colchicineUnusualWedging.add_to_png_blob(png_no_metadata_blob, JSON.stringify({
        includePkl: false, includeSmiles: true, includeMol: true,
        propertyFlags: { AtomProps: true, BondProps: true }
    }));
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert.equal(mol.has_coords(), 2);
    assert(mol.get_cxsmiles(JSON.stringify({ CX_ALL: true, restoreBondDirOption: 'RestoreBondDirOptionTrue' })).includes('wU:22.24|'));
    assert(mol.get_cxsmiles(JSON.stringify({ CX_ALL: true, restoreBondDirOption: 'RestoreBondDirOptionClear' })).includes('wU:22.24|'));
    assert(mol.get_molblock().includes(' 23 24  1  1'));
    assert(mol.get_molblock(JSON.stringify({ useMolBlockWedging: true })).includes(' 23 24  1  1'));
    mol.delete();
    // the mol is restored from CXSMILES, but restoreBondDirOption was set
    // to 'RestoreBondDirOptionTrue', so it will not retain
    // original molblock wedging, as it was stored in the CXSMILES string
    png_colchicine_metadata_blob = colchicineUnusualWedging.add_to_png_blob(png_no_metadata_blob, JSON.stringify({
        includePkl: false, includeSmiles: true, includeMol: true,
        propertyFlags: { AtomProps: true, BondProps: true },
        restoreBondDirOption: 'RestoreBondDirOptionTrue'
    }));
    mol = RDKitModule.get_mol_from_png_blob(png_colchicine_metadata_blob);
    assert(mol);
    assert.equal(mol.get_num_atoms(), 29);
    assert.equal(mol.has_coords(), 2);
    assert(mol.get_cxsmiles(JSON.stringify({ CX_ALL: true, restoreBondDirOption: 'RestoreBondDirOptionTrue' })).includes('wU:22.23|'));
    assert(mol.get_cxsmiles(JSON.stringify({ CX_ALL: true, restoreBondDirOption: 'RestoreBondDirOptionClear' })).includes('wU:22.24|'));
    assert(mol.get_molblock().includes(' 23 22  1  1'));
    assert(mol.get_molblock(JSON.stringify({ useMolBlockWedging: true })).includes(' 23 22  1  1'));
    mol.delete();
    colchicineUnusualWedging.delete();
}

function test_combine_with() {
    {
        var mol = RDKitModule.get_mol("CC");
        assert(mol);
        var other = RDKitModule.get_mol("NCC");
        assert(other);
        assert(!mol.combine_with(other));
        assert.equal(mol.get_num_atoms(), 5);
        assert(mol.get_smiles() === "CC.CCN");
        mol.delete();
        other.delete();
    }
    {
        var mol = RDKitModule.get_mol("C1CC1 |(0.866025,0,;-0.433013,0.75,;-0.433013,-0.75,)|");
        var mol_copy = RDKitModule.get_mol_copy(mol);
        assert(mol);
        var other = RDKitModule.get_mol("C1CNC1 |(-1.06066,0,;0,-1.06066,;1.06066,0,;0,1.06066,)|");
        assert(other);
        assert(!mol.combine_with(other, JSON.stringify({offset: [4.0, 0.0, 0.0]})));
        assert.equal(mol.get_num_atoms(), 7);
        assert(mol.get_smiles() === "C1CC1.C1CNC1");
        assert(mol.get_molblock().includes("    4.0000"));
        mol.delete();
        assert(!mol_copy.combine_with(other, JSON.stringify({offset: [9.0, 0.0, 0.0]})));
        assert.equal(mol_copy.get_num_atoms(), 7);
        assert(mol_copy.get_smiles() === "C1CC1.C1CNC1");
        assert(mol_copy.get_molblock().includes("    9.0000"));
        mol_copy.delete();
        other.delete();
    }
}

function test_get_coords() {
    {
        var mol = RDKitModule.get_mol("C1CC1 |(0.866025,0,;-0.433013,0.75,;-0.433013,-0.75,)|");
        assert(mol);
        assert.equal(mol.has_coords(), 2);
        var pos = mol.get_coords();
        assert(Array.isArray(pos));
        assert(pos.length === mol.get_num_atoms());
        assert(pos.every((xyz) => Array.isArray(xyz) && xyz.length === 3 && xyz.every((c) => typeof c === "number")));
        assert(JSON.stringify(pos) === "[[0.866025,0,0],[-0.433013,0.75,0],[-0.433013,-0.75,0]]");
        mol.delete();
    }
    {
        var mol = RDKitModule.get_mol("C1CC1");
        assert(mol);
        assert(!mol.has_coords());
        var pos = mol.get_coords();
        assert(Array.isArray(pos));
        assert(!pos.length);
        mol.delete();
    }
}

function test_get_v3K_v2K_molblock() {
    var mol = RDKitModule.get_mol('c1cc(O)ccc1');
    assert(mol);
    var molblock = mol.get_v3Kmolblock();
    var mol2 = RDKitModule.get_mol(molblock);
    assert(mol2);
    var smiles = mol2.get_smiles();
    assert(smiles === 'Oc1ccccc1');
    mol2.delete();
    molblock = mol.get_v3Kmolblock(JSON.stringify({ kekulize: false }));
    assert(molblock.includes('M  V30 1 4 1 2'));
    molblock = mol.get_molblock(JSON.stringify({ forceMDLVersion: 'V3000' }));
    assert(molblock.includes('V3000'));
    molblock = mol.get_v3Kmolblock(JSON.stringify({ forceMDLVersion: 'V2000' }));
    assert(molblock.includes('V3000'));
    mol.delete();
    var mol3 = RDKitModule.get_mol('N->[Pt+2](Cl)(Cl)<-N');
    assert(mol3);
    molblock = mol3.get_molblock();
    assert(molblock.includes('V3000'));
    assert(molblock.includes('M  V30 4 9 5 2'));
    molblock = mol3.get_molblock(JSON.stringify({forceMDLVersion: 'V2000'}));
    assert(molblock.includes('V2000'));
    assert(molblock.includes('  5  2  9  0'));
    molblock = mol3.get_v3Kmolblock(JSON.stringify({forceMDLVersion: 'V2000'}));
    assert(molblock.includes('V3000'));
    assert(molblock.includes('M  V30 4 9 5 2'));
    molblock = mol3.get_v2Kmolblock(JSON.stringify({forceMDLVersion: 'V3000'}));
    assert(molblock.includes('V2000'));
    assert(molblock.includes('  5  2  9  0'));
    mol3.delete();
}

function test_return_draw_coords() {
    var mb = `
     RDKit          2D

  3  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.8930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7734   -0.4465    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7734   -0.4465    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  3  1  0
M  END`;
    var reference;
    var highlight;
    var aboveTol;
    var mol = RDKitModule.get_mol(mb);
    assert(mol);
    res = mol.get_svg_with_highlights(JSON.stringify({
        width: 300,
        height: 200,
        padding: 0.2,
        returnDrawCoords: true
    }));
    assert(res);
    res = JSON.parse(res);
    assert(res.drawCoords);
    assert(res.svg);
    reference = res.drawCoords.flat();
    res = mol.get_svg_with_highlights(JSON.stringify({
        width: 300,
        height: 200,
        padding: 0.2,
        atoms: [0],
        returnDrawCoords: true
    }));
    assert(res);
    res = JSON.parse(res);
    assert(res.drawCoords);
    assert(res.svg);
    highlight = res.drawCoords.flat();
    aboveTol = reference.some((ref, i) => (Math.abs(ref - highlight[i]) > 0.1));
    assert(aboveTol);
    res = mol.get_svg_with_highlights(JSON.stringify({
        width: 300,
        height: 200,
        padding: 0.2,
        atoms: [0],
        returnDrawCoords: true,
        drawingExtentsInclude: {
            ALL: true,
            HIGHLIGHTS: false
        }
    }));
    assert(res);
    res = JSON.parse(res);
    assert(res.drawCoords);
    assert(res.svg);
    highlight = res.drawCoords.flat();
    aboveTol = reference.some((ref, i) => (Math.abs(ref - highlight[i]) > 0.1));
    assert(!aboveTol);
    mol.delete();
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
    if (RDKitModule.SubstructLibrary) {
        test_substruct_library(done);
        test_substruct_library_merge_hs();
        test_substruct_library_empty_mols();
        test_substruct_library_empty_lib();
        test_substruct_library_empty_query();
    }
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
    test_removehs();
    test_normalize_depiction();
    test_straighten_depiction();
    test_flexicanvas();
    if (RDKitModule.get_rxn)  {
        test_rxn_drawing();
        test_run_reaction();
    }
    test_legacy_stereochem();
    test_allow_non_tetrahedral_chirality();
    test_prop();
    test_highlights();
    test_add_chiral_hs();
    test_wedging_all_within_scaffold();
    test_wedging_outside_scaffold();
    test_wedging_if_no_match();
    test_get_frags();
    if (RDKitModule.Mol.prototype.get_mmpa_frags) {
        test_get_mmpa_frags();
    }
    test_hs_in_place();
    test_query_colour();
    test_alignment_r_groups_aromatic_ring();
    test_is_valid_deprecated();
    test_mol_list();
    test_get_num_atoms_bonds();
    if (RDKitModule.get_mcs_as_mol)  {
        test_mcs();
    }
    test_sanitize_no_kekulize_no_setaromaticity();
    test_partial_sanitization();
    test_capture_logs();
    test_rgroup_match_heavy_hydro_none_charged();
    test_get_sss_json();
    test_relabel_mapped_dummies();
    test_assign_cip_labels();
    test_smiles_smarts_params();
    test_wedged_bond_atropisomer();
    test_get_molblock_use_molblock_wedging();
    test_assign_chiral_tags_from_mol_parity();
    test_make_dummies_queries();
    test_get_mol_copy();
    test_multi_highlights();
    if (RDKitModule.RGroupDecomposition)  {
        test_singlecore_rgd();
        test_multicore_rgd();
    }
    test_bw_palette();
    test_custom_palette();
    if (RDKitModule.molzip)  {
        test_molzip();
    }
    test_pickle();
    test_remove_hs_details();
    test_get_mol_remove_hs();
    test_png_metadata();
    test_combine_with();
    test_get_coords();
    test_get_v3K_v2K_molblock();
    test_return_draw_coords();

    waitAllTestsFinished().then(() =>
        console.log("Tests finished successfully")
    );
});
