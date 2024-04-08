from rdkit import Chem

from rdkit.Chem import rdChemReactions

smarts = "[CH3:1][CH:2]([CH3:3])[*:4].[OH:5][CH2:6][*:7]>O=C=O>[CH3:1][CH:2]([CH3:3])[CH2:6][OH:5] |$;;;_AP1;;;_AP1;;;;;;;;$|"
smarts = "[OH:5][CH2:6][*:7].[CH3:1][CH:2]([CH3:3])[*:4]>O=C=O>[CH3:1][CH:2]([CH3:3])[CH2:6][OH:5] |$;;;_AP1;;;_AP1;;;;;;;;$|"
smarts = "[Na].[Cl]>>[Na+].[Cl-]"
smarts = "[CH3:1][CH:2]([CH3:3])[*:4].[Fe:8][OH:5][CH2:6][*:7]>>[Fe:8][OH:5][CH2:6][CH2:1][CH:2]([CH3:3])[*:4] |$;;;_AP1;;;;_AP1;;;;;;;_AP1$,C:5.3,9.6,SgD:6:foo:bar::::,SgD:10:bar:baz::::|"
# [C&H3:1][C&H1:2]([C&H3:3])[*:4].[O&H1:5][C&H2:6][*:7]>O=C=O>[C&H3:1][C&H1:2]([C&H3:3])[C&H2:6][O&H1|$;;;_AP1;;;_AP1;;;;;;;;$|
# smarts = "CO.OC1CCC(F)C1>>COC1CC(O)CC1F |LN:3:1.3.4.8,13:2.5.12.15|"

smarts = "[CH3:6][O:5][CH:3](-*)[O:2]-*>>[CH3:6][NH:5][CH:3](-*)[O:2]-* |$;;;star_e;;star_e;;;;star_e;;star_e$,SgD:1,0:foo:bar::::,SgD:7,6:foo:baz::::,Sg:n:4,2,1,0::ht,Sg:n:10,8,7,6::ht,SgH:2:0,3:1|"
# smarts = "[C&H3:6][O:5][C&H1:3](-*)[O:2]-*>>[C&H3:6][N&H1:5             |$;;;star_e;;star_e;;;;star_e;;star_e$,SgD:1,0:foo:bar::::,SgD:7,6:foo:baz::::,Sg:n:4,2,1,0::ht,Sg:n:10,8,7,6::ht:::,SgH:3:1.1|"
rxn = rdChemReactions.ReactionFromSmarts(smarts)

for m in rxn.GetReactants():
    print("Reactant:")
    for at in m.GetAtoms():
        print("\t", at.GetPropsAsDict())
    print()

for m in rxn.GetAgents():
    print("Agent:")
    for at in m.GetAtoms():
        print("\t", at.GetPropsAsDict())
    print()

for m in rxn.GetProducts():
    print("Product:")
    for at in m.GetAtoms():
        print("\t", at.GetPropsAsDict())
    print()

# print(rdChemReactions.ReactionToSmarts(rxn))
# print(rdChemReactions.ReactionToSmarts(rxn, includeCX=True))
print(rdChemReactions.ReactionToSmarts(rxn, includeCX=True))