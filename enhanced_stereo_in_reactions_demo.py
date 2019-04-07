from rdkit import Chem
from rdkit.Chem import AllChem

def get_product_cxsmiles(product):
    # Clear properties. Mapping properties show up in CXSMILES
    for a in product.GetAtoms():
        for k in a.GetPropsAsDict():
            a.ClearProp(k)
    return Chem.MolToCXSmiles(product)

def print_rxn(rxn_smarts, smiles, comment):
    print(f"\n{comment}:\n  {rxn_smarts}")
    print('   ', ' + '.join(smiles))
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    for prods in rxn.RunReactants(mols):
        print('      >>')
        print('     ', ' + '.join(map(get_product_cxsmiles, prods)))    


# StereoGroup atoms are in the reaction, but the reaction doesn't affect
# the chirality at the stereo centers
# -> preserve stereo group
print_rxn('[C@:1]>>[C@:1]', ['F[C@H](Cl)Br |o1:1|'], '0a. Reaction preserves stereo')
print_rxn('[C@:1]>>[C@:1]', ['F[C@@H](Cl)Br |&1:1|'], '0b. Reaction preserves stereo')
print_rxn('[C@:1]>>[C@:1]', ['FC(Cl)Br'],  '0c. Reaction preserves stereo')

# StereoGroup atoms are in the reaction, but the reaction doesn't specify the
# chirality at the stereo centers
# -> preserve stereo group
print_rxn('[C:1]>>[C:1]', ['F[C@H](Cl)Br |a:1|'], '1a. Reaction ignores stereo')
print_rxn('[C:1]>>[C:1]', ['F[C@@H](Cl)Br |&1:1|'], '1b. Reaction ignores stereo')
print_rxn('[C:1]>>[C:1]', ['FC(Cl)Br'],  '1c. Reaction ignores stereo')

# StereoGroup atoms are in the reaction, and the reaction inverts the specified
# chirality at the stereo centers
# -> preserve stereo group
print_rxn('[C@:1]>>[C@@:1]', ['F[C@H](Cl)Br |o1:1|'], '2a. Reaction inverts stereo')
print_rxn('[C@:1]>>[C@@:1]', ['F[C@@H](Cl)Br |&1:1|'], '2b. Reaction inverts stereo')
print_rxn('[C@:1]>>[C@@:1]', ['FC(Cl)Br'],  '2c. Reaction inverts stereo')

# StereoGroup atoms are in the reaction, and the reaction destroys the specified
# chirality at the stereo centers
# -> invalidate stereo group
print_rxn('[C@:1]>>[C:1]', ['F[C@H](Cl)Br |o1:1|'], '3a. Reaction destroys stereo')
print_rxn('[C@:1]>>[C:1]', ['F[C@@H](Cl)Br |&1:1|'], '3b. Reaction destroys stereo')
print_rxn('[C@:1]>>[C:1]', ['FC(Cl)Br'],  '3c. Reaction destroys stereo')
print_rxn('[C@:1]F>>[C:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |o1:1,&2:3|'], '3d. Reaction destroys stereo (but preserves unaffected group)' )
#### -> Should the group be invalidated or trimmed?
print_rxn('[C@:1]F>>[C:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |&1:1,3|'], '3e. Reaction destroys stereo' )

# StereoGroup atoms are in the reaction, and the reaction creates the specified
# chirality at the stereo centers
# -> invalidate stereo group
print_rxn('[C:1]>>[C@@:1]', ['F[C@H](Cl)Br |o1:1|'], '4a. Reaction creates stereo')
print_rxn('[C:1]>>[C@@:1]', ['F[C@@H](Cl)Br |&1:1|'], '4b. Reaction creates stereo')
print_rxn('[C:1]>>[C@@:1]', ['FC(Cl)Br'],  '4c. Reaction creates stereo')
print_rxn('[C:1]F>>[C@@:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |o1:1,&2:3|'], '4d. Reaction creates stereo (preserve unaffected group)' )
#### -> Should the group be invalidated or trimmed?
print_rxn('[C:1]F>>[C@@:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |o1:1,3|'], '4e. Reaction creates stereo' )


# StereoGroup atoms are not in the reaction
# -> preserve stereo group
print_rxn('[C@:1]F>>[C@:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |o1:3|'], '5a. Reaction preserves unrelated stereo' )
print_rxn('[C:1]F>>[C:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |o1:3|'], '5b. Reaction ignores unrelated stereo' )
print_rxn('[C@:1]F>>[C@@:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |o1:3|'], '5c. Reaction inverts unrelated stereo' )
print_rxn('[C@:1]F>>[C:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |o1:3|'], '5d. Reaction destroys unrelated stereo' )
print_rxn('[C:1]F>>[C@@:1]F', ['F[C@H](Cl)[C@@H](Cl)Br |o1:3|'], '5e. Reaction creates unrelated stereo' )

# StereoGroup atoms are split into two products by the reaction
#### -> Should the group be invalidated or trimmed?
print_rxn('[C:1]OO[C:2]>>[C:2]O.O[C:1]', ['F[C@H](Cl)OO[C@@H](Cl)Br |o1:1,5|'], '6e. Reaction splits StereoGroup atoms into two Mols' )

# Stereogroup atoms are in the reaction with multiple copies in the product
# -> All copies should be in the same product stereo group
print_rxn('[O:1].[C:2]=O>>[O:1][C:2][O:1]', ['Cl[C@@H](Br)C[C@H](Br)CCO |&1:1,4|', 'CC(=O)C'], '7. Add two copies')

# Stereogroup atoms are not in the reaction, but have multiple copies in the
# product.
# -> All copies should be in the same product stereo group
print_rxn('[O:1].[C:2]=O>>[O:1][C:2][O:1]', ['Cl[C@@H](Br)C[C@H](Br)CCO |&1:1,4|', 'CC(=O)C'], '8. Add two copies')

