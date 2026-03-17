// Test RDKit WASM modules in Pyodide (Node.js)
const { loadPyodide } = require("pyodide");
const fs = require("fs");
const path = require("path");

async function main() {
    console.log("Loading Pyodide...");
    const pyodide = await loadPyodide();

    const rdkitWasmDir = path.resolve(__dirname, "rdkit-wasm");

    // Copy rdkit package into Pyodide's filesystem
    function copyDirToFS(srcDir, destDir) {
        const entries = fs.readdirSync(srcDir, { withFileTypes: true });
        for (const entry of entries) {
            const srcPath = path.join(srcDir, entry.name);
            const destPath = destDir + "/" + entry.name;
            if (entry.isDirectory()) {
                try { pyodide.FS.mkdirTree(destPath); } catch(e) {}
                copyDirToFS(srcPath, destPath);
            } else if (entry.name.endsWith(".py") || entry.name.endsWith(".so") ||
                       entry.name.endsWith(".txt") || entry.name.endsWith(".fdef")) {
                const data = fs.readFileSync(srcPath);
                pyodide.FS.writeFile(destPath, data);
            }
        }
    }

    console.log("Copying RDKit files to Pyodide filesystem...");
    const sitePackages = "/lib/python3.13/site-packages";
    pyodide.FS.mkdirTree(sitePackages + "/rdkit");
    copyDirToFS(rdkitWasmDir, sitePackages + "/rdkit");

    console.log("Installing numpy...");
    await pyodide.loadPackage("numpy");

    console.log("Loading librdkit_core.so...");
    try {
        await pyodide.runPythonAsync(`
import sys
sys.path.insert(0, '/lib/python3.13/site-packages')

from pyodide_js._module import loadDynamicLibrary
import js

# Load core library with global symbols (must be loaded before any wrapper)
core_path = "/lib/python3.13/site-packages/rdkit/librdkit_core.so"
result = loadDynamicLibrary(core_path, js.JSON.parse('{"global": true, "loadAsync": true}'))
if hasattr(result, 'then'):
    await result
print("Core library loaded!")

# Basic imports
import rdkit
print(f"RDKit version: {rdkit.__version__}")

from rdkit import Chem
print("Chem imported!")

# === Test 1: SMILES parsing ===
mol = Chem.MolFromSmiles("CCO")
assert mol is not None
assert mol.GetNumAtoms() == 3
assert mol.GetNumBonds() == 2
print(f"[PASS] SMILES parsing: ethanol = {mol.GetNumAtoms()} atoms")

# === Test 2: Canonical SMILES ===
smi = Chem.MolToSmiles(Chem.MolFromSmiles("OCC"))
assert smi == "CCO"
print(f"[PASS] Canonical SMILES: OCC -> {smi}")

# === Test 3: Aromatic molecules ===
benzene = Chem.MolFromSmiles("c1ccccc1")
assert benzene is not None
assert benzene.GetNumAtoms() == 6
print(f"[PASS] Aromatic: benzene = {benzene.GetNumAtoms()} atoms")

# === Test 4: Substructure search ===
pat = Chem.MolFromSmarts("[OH]")
assert mol.HasSubstructMatch(pat)
assert not benzene.HasSubstructMatch(pat)
print("[PASS] Substructure search (SMARTS)")

# === Test 5: InChI ===
from rdkit.Chem import inchi
inchi_str = inchi.MolToInchi(mol)
assert inchi_str.startswith("InChI=")
print(f"[PASS] InChI: {inchi_str}")

# Roundtrip InChI -> mol
mol_from_inchi = inchi.MolFromInchi(inchi_str)
assert mol_from_inchi is not None
assert Chem.MolToSmiles(mol_from_inchi) == "CCO"
print("[PASS] InChI roundtrip")

# === Test 6: MolBlock (SDF) ===
aspirin = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
assert aspirin is not None
molblock = Chem.MolToMolBlock(aspirin)
assert "V2000" in molblock
aspirin2 = Chem.MolFromMolBlock(molblock)
assert aspirin2 is not None
assert Chem.MolToSmiles(aspirin) == Chem.MolToSmiles(aspirin2)
print(f"[PASS] MolBlock roundtrip: aspirin ({aspirin.GetNumAtoms()} atoms)")

# === Test 7: Ring info ===
ri = benzene.GetRingInfo()
assert ri.NumRings() == 1
print(f"[PASS] Ring info: benzene has {ri.NumRings()} ring(s)")

# === Test 8: Atom/bond properties ===
atom = mol.GetAtomWithIdx(2)  # oxygen
assert atom.GetSymbol() == "O"
assert atom.GetAtomicNum() == 8
bond = mol.GetBondWithIdx(1)  # C-O bond
assert bond.GetBondTypeAsDouble() == 1.0
print(f"[PASS] Atom/bond properties: atom 2 = {atom.GetSymbol()} (Z={atom.GetAtomicNum()})")

# === Test 9: Molecular formula ===
from rdkit.Chem import rdMolDescriptors
formula = rdMolDescriptors.CalcMolFormula(aspirin)
assert formula == "C9H8O4"
print(f"[PASS] Molecular formula: aspirin = {formula}")

# === Test 10: Chemical reactions ===
from rdkit.Chem import AllChem
rxn = AllChem.ReactionFromSmarts("[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]")
assert rxn is not None
print("[PASS] Reaction SMARTS parsing")

# === Test 11: 2D coordinate generation ===
AllChem.Compute2DCoords(benzene)
conf = benzene.GetConformer()
pos = conf.GetAtomPosition(0)
print(f"[PASS] 2D coords: atom 0 at ({pos.x:.2f}, {pos.y:.2f})")

# === Test 12: SVG drawing ===
from rdkit.Chem.Draw import rdMolDraw2D
drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
drawer.DrawMolecule(benzene)
drawer.FinishDrawing()
svg = drawer.GetDrawingText()
assert "<svg" in svg and "</svg>" in svg
print(f"[PASS] SVG drawing: {len(svg)} chars")

# === Test 13: SMILES of drug-like molecules ===
drugs = {
    "caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "penicillin_G": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
}
for name, smi in drugs.items():
    m = Chem.MolFromSmiles(smi)
    assert m is not None, f"Failed to parse {name}"
    can = Chem.MolToSmiles(m)
    m2 = Chem.MolFromSmiles(can)
    assert Chem.MolToSmiles(m2) == can
    print(f"[PASS] {name}: {m.GetNumAtoms()} atoms")

# === Test 14: Boost serialization (ToBinary/FromBinary) ===
aspirin_bin = aspirin.ToBinary()
assert len(aspirin_bin) > 0
aspirin3 = Chem.Mol(aspirin_bin)
assert Chem.MolToSmiles(aspirin3) == Chem.MolToSmiles(aspirin)
print(f"[PASS] Binary serialization roundtrip: {len(aspirin_bin)} bytes")

# === Test 15: Fingerprint generators with numpy ===
import numpy as np
from rdkit.Chem import rdFingerprintGenerator
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
fp_np = fpgen.GetFingerprintAsNumPy(mol)
assert isinstance(fp_np, np.ndarray)
assert fp_np.shape[0] > 0
print(f"[PASS] Fingerprint as NumPy array: shape={fp_np.shape}, dtype={fp_np.dtype}")

# === Test 16: 3D descriptors ===
from rdkit.Chem import Descriptors3D, AllChem
mol3d = Chem.MolFromSmiles("CCCC")
mol3d = Chem.AddHs(mol3d)
AllChem.EmbedMolecule(mol3d, randomSeed=42)
asphericity = Descriptors3D.Asphericity(mol3d)
print(f"[PASS] 3D descriptors: Asphericity(butane) = {asphericity:.4f}")

# === Test 17: 3D embedding ===
embed_mol = Chem.MolFromSmiles("c1ccc(O)cc1")
embed_mol = Chem.AddHs(embed_mol)
res = AllChem.EmbedMolecule(embed_mol, randomSeed=42)
assert res == 0, f"EmbedMolecule failed with code {res}"
conf = embed_mol.GetConformer()
assert conf.Is3D()
pos0 = conf.GetAtomPosition(0)
assert not (pos0.x == 0.0 and pos0.y == 0.0 and pos0.z == 0.0)
# Optimize with UFF
res_opt = AllChem.UFFOptimizeMolecule(embed_mol, maxIters=200)
assert res_opt == 0, f"UFF optimization failed with code {res_opt}"
print(f"[PASS] 3D embedding + UFF optimization: {embed_mol.GetNumAtoms()} atoms, 3D={conf.Is3D()}")

# === Test 18: Add/remove hydrogens ===
phenol = Chem.MolFromSmiles("c1ccc(O)cc1")
assert phenol.GetNumAtoms() == 7  # heavy atoms only
phenol_h = Chem.AddHs(phenol)
assert phenol_h.GetNumAtoms() == 13  # 7 heavy + 6 H
phenol_noh = Chem.RemoveHs(phenol_h)
assert phenol_noh.GetNumAtoms() == 7
# Verify SMILES roundtrip is preserved
assert Chem.MolToSmiles(phenol) == Chem.MolToSmiles(phenol_noh)
print(f"[PASS] Add/RemoveHs: {phenol.GetNumAtoms()} -> {phenol_h.GetNumAtoms()} -> {phenol_noh.GetNumAtoms()} atoms")

# === Test 19: Kekulization ===
arom = Chem.MolFromSmiles("c1ccccc1")
# Check aromatic bonds
bond = arom.GetBondWithIdx(0)
assert bond.GetIsAromatic()
assert bond.GetBondType() == Chem.rdchem.BondType.AROMATIC
# Kekulize
Chem.Kekulize(arom, clearAromaticFlags=True)
bond_k = arom.GetBondWithIdx(0)
assert not bond_k.GetIsAromatic()
assert bond_k.GetBondType() in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE)
kek_smi = Chem.MolToSmiles(arom, kekuleSmiles=True)
assert "c" not in kek_smi  # no lowercase = no aromatic
print(f"[PASS] Kekulization: c1ccccc1 -> {kek_smi}")

# === Test 20: MMFF optimization ===
mmff_mol = Chem.MolFromSmiles("c1ccc(O)cc1")
mmff_mol = Chem.AddHs(mmff_mol)
AllChem.EmbedMolecule(mmff_mol, randomSeed=42)
props = AllChem.MMFFGetMoleculeProperties(mmff_mol)
assert props is not None
ff = AllChem.MMFFGetMoleculeForceField(mmff_mol, props)
assert ff is not None
e_before = ff.CalcEnergy()
res_mmff = AllChem.MMFFOptimizeMolecule(mmff_mol, maxIters=200)
assert res_mmff == 0
ff2 = AllChem.MMFFGetMoleculeForceField(mmff_mol, AllChem.MMFFGetMoleculeProperties(mmff_mol))
e_after = ff2.CalcEnergy()
print(f"[PASS] MMFF optimization: energy {e_before:.1f} -> {e_after:.1f}")

# === Test 21: Multiple conformers (ETKDGv3) ===
conf_mol = Chem.MolFromSmiles("CCCCCCC")
conf_mol = Chem.AddHs(conf_mol)
params = AllChem.ETKDGv3()
params.randomSeed = 42
params.numThreads = 1
cids = AllChem.EmbedMultipleConfs(conf_mol, numConfs=5, params=params)
assert len(cids) == 5
assert conf_mol.GetNumConformers() == 5
print(f"[PASS] EmbedMultipleConfs (ETKDGv3): {len(cids)} conformers generated")

# === Test 22: Molecular alignment (rdMolAlign) ===
from rdkit.Chem import rdMolAlign
# Align conformers of the same molecule
rmsds = []
rdMolAlign.AlignMolConformers(conf_mol, RMSlist=rmsds)
assert len(rmsds) > 0
# Align two different molecules
ref = Chem.MolFromSmiles("c1ccccc1O")
ref = Chem.AddHs(ref)
AllChem.EmbedMolecule(ref, randomSeed=42)
probe = Chem.MolFromSmiles("c1ccccc1O")
probe = Chem.AddHs(probe)
AllChem.EmbedMolecule(probe, randomSeed=123)
rmsd = rdMolAlign.AlignMol(probe, ref)
print(f"[PASS] Molecular alignment: conformer RMSDs={[f'{r:.2f}' for r in rmsds[:3]]}, AlignMol RMSD={rmsd:.4f}")

# === Test 23: Maximum Common Substructure (rdFMCS) ===
from rdkit.Chem import rdFMCS
mol_a = Chem.MolFromSmiles("c1ccccc1CCO")
mol_b = Chem.MolFromSmiles("c1ccccc1CCCO")
mcs = rdFMCS.FindMCS([mol_a, mol_b])
assert mcs.numAtoms > 0
assert mcs.numBonds > 0
print(f"[PASS] MCS: {mcs.numAtoms} atoms, {mcs.numBonds} bonds, SMARTS={mcs.smartsString}")

# === Test 24: Tautomer enumeration ===
from rdkit.Chem.MolStandardize import rdMolStandardize
taut_enum = rdMolStandardize.TautomerEnumerator()
keto = Chem.MolFromSmiles("OC1=CC=CC=C1")  # phenol with explicit keto form
canonical = taut_enum.Canonicalize(keto)
assert canonical is not None
tautomers = list(taut_enum.Enumerate(keto))
assert len(tautomers) >= 1
print(f"[PASS] Tautomer enumeration: {len(tautomers)} tautomer(s), canonical={Chem.MolToSmiles(canonical)}")

# === Test 25: Salt removal (using default Data/Salts.txt) ===
from rdkit.Chem.SaltRemover import SaltRemover
remover = SaltRemover()  # uses RDConfig.RDDataDir/Salts.txt
salt_mol = Chem.MolFromSmiles("[Na+].OC1=CC=CC=C1")  # sodium phenolate
stripped = remover.StripMol(salt_mol)
assert stripped is not None
assert stripped.GetNumAtoms() < salt_mol.GetNumAtoms()
print(f"[PASS] Salt removal: {salt_mol.GetNumAtoms()} -> {stripped.GetNumAtoms()} atoms, SMILES={Chem.MolToSmiles(stripped)}")

# === Test 26: Stereochemistry assignment ===
chiral = Chem.MolFromSmiles("C[C@@H](O)F")
Chem.AssignStereochemistry(chiral, cleanIt=True, force=True)
stereo_atom = chiral.GetAtomWithIdx(1)
assert stereo_atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED
print(f"[PASS] Stereochemistry: chiral tag = {stereo_atom.GetChiralTag()}")

# AssignStereochemistryFrom3D
mol_3d = Chem.MolFromSmiles("C[C@@H](O)F")
mol_3d = Chem.AddHs(mol_3d)
AllChem.EmbedMolecule(mol_3d, randomSeed=42)
Chem.AssignStereochemistryFrom3D(mol_3d)
print(f"[PASS] Stereochemistry from 3D: SMILES={Chem.MolToSmiles(Chem.RemoveHs(mol_3d))}")

# === Test 27: RWMol (editable molecule) and CombineMols ===
rwmol = Chem.RWMol(Chem.MolFromSmiles("C"))
idx = rwmol.AddAtom(Chem.Atom(8))  # add oxygen
rwmol.AddBond(0, idx, Chem.rdchem.BondType.SINGLE)
Chem.SanitizeMol(rwmol)
assert Chem.MolToSmiles(rwmol) == "CO"
# CombineMols
m1 = Chem.MolFromSmiles("C")
m2 = Chem.MolFromSmiles("O")
combined = Chem.CombineMols(m1, m2)
assert combined.GetNumAtoms() == 2
print(f"[PASS] RWMol edit: C + O -> {Chem.MolToSmiles(rwmol)}, CombineMols: {combined.GetNumAtoms()} atoms")

# === Test 28: Descriptors.descList ===
from rdkit.Chem import Descriptors
assert len(Descriptors.descList) > 0
# Calculate a few descriptors
mw = Descriptors.MolWt(aspirin)
logp = Descriptors.MolLogP(aspirin)
hbd = Descriptors.NumHDonors(aspirin)
hba = Descriptors.NumHAcceptors(aspirin)
print(f"[PASS] Descriptors: MW={mw:.1f}, LogP={logp:.2f}, HBD={hbd}, HBA={hba}, total={len(Descriptors.descList)} descriptors")

# === Test 29: MACCS fingerprints ===
from rdkit.Chem import MACCSkeys
maccs = MACCSkeys.GenMACCSKeys(aspirin)
assert maccs is not None
assert maccs.GetNumOnBits() > 0
print(f"[PASS] MACCS fingerprints: {maccs.GetNumOnBits()} on bits out of {maccs.GetNumBits()}")

# === Test 30: SDF reading/writing ===
import io
# Write to SDF
sdf_out = Chem.SDWriter("/tmp/test.sdf")
sdf_out.write(aspirin)
sdf_out.write(benzene)
sdf_out.close()
# Read back
suppl = Chem.SDMolSupplier("/tmp/test.sdf")
mols_read = [m for m in suppl if m is not None]
assert len(mols_read) == 2
assert Chem.MolToSmiles(mols_read[0]) == Chem.MolToSmiles(aspirin)
print(f"[PASS] SDF read/write: wrote 2 mols, read back {len(mols_read)}")

# === Test 31: RDConfig.RDDataDir and ChemicalFeatures ===
from rdkit import RDConfig
import os
assert os.path.isdir(RDConfig.RDDataDir), f"RDDataDir not found: {RDConfig.RDDataDir}"
fdef_path = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
assert os.path.isfile(fdef_path), f"BaseFeatures.fdef not found: {fdef_path}"
from rdkit.Chem import ChemicalFeatures
feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
assert feat_factory is not None
aspirin_h = Chem.AddHs(aspirin)
AllChem.EmbedMolecule(aspirin_h, randomSeed=42)
feats = feat_factory.GetFeaturesForMol(aspirin_h)
assert len(feats) > 0
feat_types = set(f.GetFamily() for f in feats)
print(f"[PASS] ChemicalFeatures: {len(feats)} features, families={feat_types}")

print(f"\\n=== ALL 31 TESTS PASSED === RDKit {rdkit.__version__} works in Pyodide!")
`);
    } catch (e) {
        console.error("FAILED:", e.message);
        process.exit(1);
    }
}

main().catch(console.error);
