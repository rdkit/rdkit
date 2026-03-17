"""Test RDKit in Pyodide - run with: node test_pyodide.js"""
import sys


async def main():
    try:
        from pyodide_js._module import loadDynamicLibrary
        import js

        # Load core library with global symbols (must happen before any wrapper import)
        core_path = "/lib/python3.13/site-packages/rdkit/librdkit_core.so"
        result = loadDynamicLibrary(
            core_path, js.JSON.parse('{"global": true, "loadAsync": true}')
        )
        if hasattr(result, "then"):
            await result
        print("Core library loaded!")

        import rdkit

        print(f"RDKit version: {rdkit.__version__}")

        from rdkit import Chem

        # Test basic SMILES parsing
        mol = Chem.MolFromSmiles("CCO")
        if mol is None:
            print("FAIL: Could not parse CCO")
            return 1

        print(f"Ethanol atoms: {mol.GetNumAtoms()}")
        print(f"Canonical SMILES: {Chem.MolToSmiles(mol)}")

        # Test aromatic molecule
        mol2 = Chem.MolFromSmiles("c1ccccc1")
        if mol2:
            print(f"Benzene atoms: {mol2.GetNumAtoms()}")

        print("SUCCESS: Basic RDKit tests passed!")
        return 0

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    import asyncio

    sys.exit(asyncio.run(main()))
