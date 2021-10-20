//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "Filters.h"
#include "FilterCatalog.h"

namespace RDKit {

/////////////////////////////////////////////////////////////////////////////////////////
// BRENK data
// # Reference: Brenk R et al. Lessons Learnt from Assembling Screening
// Libraries for Drug Discovery for Neglected Diseases. ChemMedChem 3 (2008)
// 435-444. doi:10.1002/cmdc.200700139.
// # Scope: unwanted functionality due to potential tox reasons or unfavourable
// pharmacokinetic properties
//

const FilterData_t BRENK[] = {
    {">_2_ester_groups", "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]", 0, ""},
    {"2-halo_pyridine", "n1c([F,Cl,Br,I])cccc1", 0, ""},
    {"acid_halide", "C(=O)[Cl,Br,I,F]", 0, ""},
    {"acyclic_C=C-O", "C=[C!r]O", 0, ""},
    {"acyl_cyanide", "N#CC(=O)", 0, ""},
    {"acyl_hydrazine", "C(=O)N[NH2]", 0, ""},
    {"aldehyde", "[CH1](=O)", 0, ""},
    {"Aliphatic_long_chain", "[R0;D2][R0;D2][R0;D2][R0;D2]", 0, ""},
    {"alkyl_halide", "[CX4][Cl,Br,I]", 0, ""},
    {"amidotetrazole", "c1nnnn1C=O", 0, ""},
    {"aniline", "c1cc([NH2])ccc1", 0, ""},
    {"azepane", "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1", 0, ""},
    {"Azido_group", "N=[N+]=[N-]", 0, ""},
    {"Azo_group", "N#N", 0, ""},
    {"azocane", "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1", 0, ""},
    {"benzidine",
     "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2](["
     "Nv3X3,Nv4X4])[cR2][cR2]2",
     0, ""},
    {"beta-keto/anhydride", "[C,c](=O)[CX4,CR0X3,O][C,c](=O)", 0, ""},
    {"biotin_analogue", "C12C(NC(N1)=O)CSC2", 0, ""},
    {"Carbo_cation/anion", "[C+,c+,C-,c-]", 0, ""},
    {"catechol", "c1c([OH])c([OH,NH2,NH])ccc1", 0, ""},
    {"charged_oxygen_or_sulfur_atoms", "[O+,o+,S+,s+]", 0, ""},
    {"chinone_1", "C1(=[O,N])C=CC(=[O,N])C=C1", 0, ""},
    {"chinone_2", "C1(=[O,N])C(=[O,N])C=CC=C1", 0, ""},
    {"conjugated_nitrile_group", "C=[C!r]C#N", 0, ""},
    {"crown_ether", "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]", 0,
     ""},
    {"cumarine", "c1ccc2c(c1)ccc(=O)o2", 0, ""},
    {"cyanamide", "N[CH2]C#N", 0, ""},
    {"cyanate_/aminonitrile_/thiocyanate", "[N,O,S]C#N", 0, ""},
    {"cyanohydrins", "N#CC[OH]", 0, ""},
    {"cycloheptane_1", "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2]1", 0, ""},
    {"cycloheptane_2", "[CR2]1[CR2][CR2]cc[CR2][CR2]1", 0, ""},
    {"cyclooctane_1", "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1", 0, ""},
    {"cyclooctane_2", "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1", 0, ""},
    {"diaminobenzene_1",
     "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1", 0, ""},
    {"diaminobenzene_2",
     "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1", 0, ""},
    {"diaminobenzene_3",
     "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])", 0, ""},
    {"diazo_group", "[N!R]=[N!R]", 0, ""},
    {"diketo_group", "[C,c](=O)[C,c](=O)", 0, ""},
    {"disulphide", "SS", 0, ""},
    {"enamine", "[CX2R0][NX3R0]", 0, ""},
    {"ester_of_HOBT", "C(=O)Onnn", 0, ""},
    {"four_member_lactones", "C1(=O)OCC1", 0, ""},
    {"halogenated_ring_1", "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]", 0,
     ""},
    {"halogenated_ring_2", "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]", 0,
     ""},
    {"heavy_metal", "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si]", 0, ""},
    {"het-C-het_not_in_ring",
     "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]", 0, ""},
    {"hydantoin", "C1NC(=O)NC(=O)1", 0, ""},
    {"hydrazine", "N[NH2]", 0, ""},
    {"hydroquinone", "[OH]c1ccc([OH,NH2,NH])cc1", 0, ""},
    {"hydroxamic_acid", "C(=O)N[OH]", 0, ""},
    {"imine_1", "C=[N!R]", 0, ""},
    {"imine_2", "N=[CR0][N,n,O,S]", 0, ""},
    {"iodine", "I", 0, ""},
    {"isocyanate", "N=C=O", 0, ""},
    {"isolated_alkene",
     "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4]"
     ")[CX4])]",
     0, ""},
    {"ketene", "C=C=O", 0, ""},
    {"methylidene-1,3-dithiole", "S1C=CSC1=S", 0, ""},
    {"Michael_acceptor_1", "C=!@CC=[O,S]", 0, ""},
    {"Michael_acceptor_2", "[$([CH]),$(CC)]#CC(=O)[C,c]", 0, ""},
    {"Michael_acceptor_3", "[$([CH]),$(CC)]#CS(=O)(=O)[C,c]", 0, ""},
    {"Michael_acceptor_4", "C=C(C=O)C=O", 0, ""},
    {"Michael_acceptor_5", "[$([CH]),$(CC)]#CC(=O)O[C,c]", 0, ""},
    {"N_oxide", "[NX2,nX3][OX1]", 0, ""},
    {"N-acyl-2-amino-5-mercapto-1,3,4-_thiadiazole", "s1c(S)nnc1NC=O", 0, ""},
    {"N-C-halo", "NC[F,Cl,Br,I]", 0, ""},
    {"N-halo", "[NX3,NX4][F,Cl,Br,I]", 0, ""},
    {"N-hydroxyl_pyridine", "n[OH]", 0, ""},
    {"nitro_group", "[N+](=O)[O-]", 0, ""},
    {"N-nitroso", "[#7]-N=O", 0, ""},
    {"oxime_1", "[C,c]=N[OH]", 0, ""},
    {"oxime_2", "[C,c]=NOC=O", 0, ""},
    {"Oxygen-nitrogen_single_bond", "[OR0,NR0][OR0,NR0]", 0, ""},
    {"Perfluorinated_chain", "[CX4](F)(F)[CX4](F)F", 0, ""},
    {"peroxide", "OO", 0, ""},
    {"phenol_ester", "c1ccccc1OC(=O)[#6]", 0, ""},
    {"phenyl_carbonate", "c1ccccc1OC(=O)O", 0, ""},
    {"phosphor", "P", 0, ""},
    {"phthalimide", "[cR,CR]~C(=O)NC(=O)~[cR,CR]", 0, ""},
    {"Polycyclic_aromatic_hydrocarbon_1", "a1aa2a3a(a1)A=AA=A3=AA=A2", 0, ""},
    {"Polycyclic_aromatic_hydrocarbon_2", "a21aa3a(aa1aaaa2)aaaa3", 0, ""},
    {"Polycyclic_aromatic_hydrocarbon_3", "a31a(a2a(aa1)aaaa2)aaaa3", 0, ""},
    {"polyene", "[CR0]=[CR0][CR0]=[CR0]", 0, ""},
    {"quaternary_nitrogen_1",
     "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]", 0, ""},
    {"quaternary_nitrogen_2",
     "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,"
     "C,n,N,o,O]",
     0, ""},
    {"quaternary_nitrogen_3", "[*]=[N+]=[*]", 0, ""},
    {"saponine_derivative", "O1CCCCC1OC2CCC3CCCCC3C2", 0, ""},
    {"silicon_halogen", "[Si][F,Cl,Br,I]", 0, ""},
    {"stilbene", "c1ccccc1C=Cc2ccccc2", 0, ""},
    {"sulfinic_acid", "[SX3](=O)[O-,OH]", 0, ""},
    {"Sulfonic_acid_1", "[C,c]S(=O)(=O)O[C,c]", 0, ""},
    {"Sulfonic_acid_2", "S(=O)(=O)[O-,OH]", 0, ""},
    {"sulfonyl_cyanide", "S(=O)(=O)C#N", 0, ""},
    {"sulfur_oxygen_single_bond", "[SX2]O", 0, ""},
    {"sulphate", "OS(=O)(=O)[O-]", 0, ""},
    {"sulphur_nitrogen_single_bond", "[SX2H0][N]", 0, ""},
    {"Thiobenzothiazole_1", "c12ccccc1(SC(S)=N2)", 0, ""},
    {"thiobenzothiazole_2", "c12ccccc1(SC(=S)N2)", 0, ""},
    {"Thiocarbonyl_group", "[C,c]=S", 0, ""},
    {"thioester", "SC=O", 0, ""},
    {"thiol_1", "[S-]", 0, ""},
    {"thiol_2", "[SH]", 0, ""},
    {"Three-membered_heterocycle", "*1[O,S,N]*1", 0, ""},
    {"triflate", "OS(=O)(=O)C(F)(F)F", 0, ""},
    {"triphenyl_methyl-silyl", "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)", 0,
     ""},
    {"triple_bond", "C#C", 0, ""}};
const unsigned int NUM_BRENK =
    static_cast<unsigned int>(sizeof(BRENK) / sizeof(FilterData_t));

const FilterProperty_t BRENK_PROPS[] = {
    {"FilterSet", "Brenk"},
    {"Reference",
     "Brenk R et al. Lessons Learnt from Assembling Screening Libraries for "
     "Drug Discovery for Neglected Diseases. ChemMedChem 3 (2008) 435-444. "
     "doi:10.1002/cmdc.200700139."},
    {"Scope",
     "unwanted functionality due to potential tox reasons or unfavourable "
     "pharmacokinetic properties"}};
const unsigned int NUM_BRENK_PROPS =
    static_cast<unsigned int>(sizeof(BRENK_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// NIH data
// # Scope: annotate compounds with problematic functional groups
// # Reference: Doveston R, et al. A Unified Lead-oriented Synthesis of over
// Fifty Molecular Scaffolds. Org Biomol Chem 13 (2014) 859D65.
// doi:10.1039/C4OB02287D.
// # Reference: Jadhav A, et al. Quantitative Analyses of Aggregation,
// Autofluorescence, and Reactivity Artifacts in a Screen for Inhibitors of a
// Thiol Protease. J Med Chem 53 (2009) 37D51. doi:10.1021/jm901070c.
//

const FilterData_t NIH[] = {
    {"2halo_pyrazine_3EWG",
     "[#7;R1]1[#6]([F,Cl,Br,I])[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)("
     "=O)),$([N+](=O)[O-]),$(C=O)])[#7][#6][#6]1",
     0, ""},
    {"2halo_pyrazine_5EWG",
     "[#7;R1]1[#6]([F,Cl,Br,I])[#6;!$(c-N)][#7][#6]([$(S(=O)(=O)),$(C(F)(F)(F))"
     ",$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C=O)])[#6;!$(c-N)]1",
     0, ""},
    {"2halo_pyridazine_3EWG",
     "[#7;R1]1[#6]([F,Cl,Br,I])[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)("
     "=O)),$([N+](=O)[O-]),$(C=O)])[#6][#6][#7]1",
     0, ""},
    {"2halo_pyridazine_5EWG",
     "[#7;R1]1[#6]([F,Cl,Br,I])[#6][#6][#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),"
     "$(N(=O)(=O)),$([N+](=O)[O-]),$(C=O)])[#7]1",
     0, ""},
    {"2halo_pyridine_3EWG",
     "[#7;R1]1[#6;!$(c=O)]([F,Cl,Br,I])[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),"
     "$(N(=O)(=O)),$([N+](=O)[O-]),$(C=O)])[#6;!$(c-N)][#6][#6;!$(c-N)]1",
     0, ""},
    {"2halo_pyridine_5EWG",
     "[#7;R1]1[#6;!$(c=O)]([F,Cl,Br,I])[#6][#6;!$(c-N)][#6]([$(S(=O)(=O)),$(C("
     "F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C=O)])[#6;!$(c=O);!$(c-N)"
     "]1",
     0, ""},
    {"2halo_pyrimidine_5EWG",
     "[#7;R1]1[#6]([F,Cl,Br,I])[#7][#6][#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),"
     "$(N(=O)(=O)),$([N+](=O)[O-]),$(C=O)])[#6]1",
     0, ""},
    {"3halo_pyridazine_2EWG",
     "[#7;R1]1[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)["
     "O-]),$(C=O)])[#6]([F,Cl,Br,I])[#6][#6][#7]1",
     0, ""},
    {"3halo_pyridazine_4EWG",
     "[#7;R1]1[#6][#6]([F,Cl,Br,I])[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N("
     "=O)(=O)),$([N+](=O)[O-]),$(C=O)])[#6][#7]1",
     0, ""},
    {"4_pyridone_3_5_EWG",
     "[#7,#8,#16]1~[#6;H]~[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),"
     "$([N+](=O)[O-]),$(C=O)])~[#6](=O)~[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N)"
     ",$(N(=O)(=O)),$([N+](=O)[O-]),$(C=O)])~[#6;H]1",
     0, ""},
    {"4halo_pyridine_3EWG",
     "[#7;R1]1[#6;!$(c=O);!$(c-N)][#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(="
     "O)(=O)),$([N+](=O)[O-]),$(C=O)])[#6]([F,Cl,Br,I])[#6][#6;!$(c=O);!$(c-N)]"
     "1",
     0, ""},
    {"4halo_pyrimidine_2_6EWG",
     "[#7]1[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-])"
     ",$(C=O)])[#7;R1][#6]([F,Cl,Br,I])[#6][#6]1([$(S(=O)(=O)),$(C(F)(F)(F)),$("
     "C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C=O)])",
     0, ""},
    {"4halo_pyrimidine_5EWG",
     "[#7]1[#6][#7;R1][#6]([F,Cl,Br,I])[#6]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),"
     "$(N(=O)(=O)),$([N+](=O)[O-]),$(C=O)])[#6]1",
     0, ""},
    {"CH2_S#O_3_ring", "[CH2]1[O,S]C1", 0, ""},
    {"HOBT_ester", "O=C(-[!N])O[$(nnn),$([#7]-[#7]=[#7])]", 0, ""},
    {"NO_phosphonate", "P(=O)ON", 0, ""},
    {"acrylate", "[CH2]=[C;!$(C-N);!$(C-O)]C(=O)", 0, ""},
    {"activated_4mem_ring",
     "[#6]1~[$(C(=O)),$(S(=O))]~[O,S,N]~[$(C(=O)),$(S(=O))]1", 0, ""},
    {"activated_S#O_3_ring", "C1~[O,S]~[C,N,O,S]1[a,N,O,S]", 0, ""},
    {"activated_acetylene",
     "[$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O))]"
     "C#[C;!$(C-N);!$(C-n)]",
     0, ""},
    {"activated_diazo",
     "[N;!R]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$("
     "C(=O))])=[N;!R]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](="
     "O)[O-]),$(C(=O))])",
     0, ""},
    {"activated_vinyl_ester",
     "O=COC=[$(C(S(=O)(=O))),$(C(C(F)(F)(F))),$(C(C#N)),$(C(N(=O)(=O))),$(C([N+"
     "](=O)[O-])),$(C(C(=O)));!$(C(N))]",
     0, ""},
    {"activated_vinyl_sulfonate",
     "O(-S(=O)(=O))C=[$(C(S(=O)(=O))),$(C(C(F)(F)(F))),$(C(C#N)),$(C(N(=O)(=O))"
     "),$(C([N+](=O)[O-])),$(C(C(=O)));!$(C(N))]",
     0, ""},
    {"acyclic_imide", "[C,c][C;!R](=O)[N;!R][C;!R](=O)[C,c]", 0, ""},
    {"acyl_123_triazole", "[#7;R1]1~[#7;R1]~[#7;R1](-C(=O))~[#6]~[#6]1", 0, ""},
    {"acyl_134_triazole", "[#7]1~[#7]~[#6]~[#7](-C(=O)[!N])~[#6]1", 0, ""},
    {"acyl_activated_NO", "O=C(-[!N])O[$([#7;+]),$(N(C=[O,S,N])(C=[O,S,N]))]",
     0, ""},
    {"acyl_cyanide", "C(=O)-C#N", 0, ""},
    {"acyl_imidazole",
     "[C;!$(C-N)](=O)[#7]1[#6;H1,$([#6]([*;!R]))][#7][#6;H1,$([#6]([*;!R]))][#"
     "6;H1,$([#6]([*;!R]))]1",
     0, ""},
    {"acyl_pyrazole",
     "[C;!$(C-N)](=O)[#7]1[#7][#6;H1,$([#6]([*;!R]))][#6;H1,$([#6]([*;!R]))][#"
     "6;H1,$([#6]([*;!R]))]1",
     0, ""},
    {"aldehyde", "[C,c][C;H1](=O)", 0, ""},
    {"alpha_dicarbonyl", "C(=O)!@C(=O)", 0, ""},
    {"alpha_halo_EWG",
     "[$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-])]-[CH,CH2]-[Cl,Br,I,$("
     "O(S(=O)(=O)))]",
     0, ""},
    {"alpha_halo_amine",
     "[F,Cl,Br,I,$(O(S(=O)(=O)))]-[CH,CH2;!$(C(F)(F))]-[N,n]", 0, "Edited"},
    {"alpha_halo_carbonyl", "C(=O)([CH,CH2][Cl,Br,I,$(O(S(=O)(=O)))])", 0, ""},
    {"alpha_halo_heteroatom",
     "[N,n,O,S;!$(S(=O)(=O))]-[CH,CH2;!$(C(F)(F))][F,Cl,Br,I,$(O(S(=O)(=O)))]",
     0, ""},
    {"alpha_halo_heteroatom_tert",
     "[N,n,O,S;!$(S(=O)(=O))]-C([Cl,Br,I,$(O(S(=O)(=O)))])(C)(C)", 0, ""},
    {"anhydride",
     "[$(C(=O)),$(C(=S))]-[O,S]-[$(C(=O)),$(C(=S)),$(C(=[N;!R])),$(C(=N(-[C;X4]"
     ")))]",
     0, ""},
    {"aryl_phosphonate", "P(=O)-[O;!R]-a", 0, ""},
    {"aryl_thiocarbonyl", "a-[S;X2;!R]-[C;!R](=O)", 0, ""},
    {"azide", "[$(N#[N+]-[N-]),$([N-]=[N+]=N)]", 0, ""},
    {"aziridine_diazirine", "[C,N]1~[C,N]~N~1", 0, ""},
    {"azo_amino", "[N]=[N;!R]-[N]", 0, ""},
    {"azo_aryl", "c[N;!R;!+]=[N;!R;!+]-c", 0, ""},
    {"azo_filter1", "[N;!R]=[N;!R]-[N]=[*]", 0, ""},
    {"azo_filter2",
     "[N;!$(N-S(=O)(=O));!$(N-C=O)]-[N;!r3;!$(N-S(=O)(=O));!$(N-C=O)]-[N;!$(N-"
     "S(=O)(=O));!$(N-C=O)]",
     0, ""},
    {"azo_filter3", "[N;!R]-[N;!R]-[N;!R]", 0, ""},
    {"azo_filter4", "a-N=N-[N;H2]", 0, ""},
    {"bad_boron", "[B-,BH2,BH3,$(B(F)(F))]", 0, ""},
    {"bad_cations", "[C+,F+,Cl+,Br+,I+,Se+]", 0, ""},
    {"benzidine_like", "c([N;!+])1ccc(c2ccc([N;!+])cc2)cc1", 0, ""},
    {"beta_lactone", "[#6,#15,#16]1(=O)~[#6]~[#6]~[#8,#16]1", 0, ""},
    {"betalactam", "C1(=O)~[#6]~[#6]N1", 0, ""},
    {"betalactam_EWG",
     "C1(=O)~[#6]~[#6]N1([$(S(=O)(=O)[C,c,O&D2]),$(C(F)(F)(F)),$(C#N),$(N(=O)(="
     "O)),$([N+](=O)[O-]),$(C(=O)[C,c,O&D2])])",
     0, ""},
    {"bis_activated_aryl_ester",
     "O=[C,S]Oc1aaa([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)["
     "O-]),$(C(=O)O),$(C(=O)N)])aa([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(="
     "O)),$([N+](=O)[O-]),$(C(=O)O),$(C(=O)N)])1",
     0, ""},
    {"bis_keto_olefin",
     "CC(=O)[$([C&H1]),$(C-F),$(C-Cl),$(C-Br),$(C-I)]=[$([C&H1]),$(C-F),$(C-Cl)"
     ",$(C-Br),$(C-I)]C(=O)C",
     0, ""},
    {"boron_warhead", "[C,c]~[#5]", 0, ""},
    {"branched_polycyclic_aromatic", "a1(a2aa(a3aaaaa3)aa(a4aaaaa4)a2)aaaaa1",
     0, ""},
    {"carbodiimide_iso#thio#cyanate", "N=C=[N,O,S]", 0, ""},
    {"carbonyl_halide", "O=C[F,Cl,Br,I]", 0, ""},
    {"contains_metal",
     "[$([Ru]),$([#45]),$([Se]),$([se]),$([Pd]),$([#21]),$([Bi]),$([Sb]),$([Ag]"
     "),$([Ti]),$([Al]),$([Cd]),$([V]),$([In]),$([#24]),$([#50]),$([Mn]),$([La]"
     "),$([Fe]),$([Er]),$([Tm]),$([Yb]),$([Lu]),$([Hf]),$([Ta]),$([W]),$([Re]),"
     "$([#27]),$([#76]),$([Ni]),$([Ir]),$([Cu]),$([Zn]),$([Ga]),$([Ge]),$([As])"
     ",$([Y]),$([Zr]),$([Nb]),$([Ce]),$([#59]),$([Nd]),$([Sm]),$([Eu]),$([Gd]),"
     "$([Tb]),$([Dy]),$([#67]),$([Pt]),$([Au]),$([Hg]),$([Tl]),$([Pb]),$([Ac]),"
     "$([Th]),$([Pa]),$([Mo]),$([U]),$([Tc]),$([Te]),$([#84]),$([At])]",
     0, "Edited"},
    {"crown_ether",
     "[$([O,S,#7;R1;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,"
     "r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,r13,r14,r15,r16,r17,"
     "r18][O,S,#7;R1;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,"
     "r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,r13,r14,r15,r16,r17,"
     "r18][O,S,#7;R1;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18]),$([O,S,#7;R1;r9,"
     "r10,r11,r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,r13,r14,r15,"
     "r16,r17,r18][CH,CH2;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,"
     "r10,r11,r12,r13,r14,r15,r16,r17,r18][O,S,#7;R1;r9,r10,r11,r12,r13,r14,"
     "r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18][CH,CH2;"
     "r9,r10,r11,r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,r13,r14,"
     "r15,r16,r17,r18][O,S,#7;R1;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18]),$([O,"
     "S,#7;R1;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,"
     "r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18]["
     "O,S,#7;R1;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,"
     "r13,r14,r15,r16,r17,r18][CH,CH2;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18]["
     "CH,CH2;r9,r10,r11,r12,r13,r14,r15,r16,r17,r18][O,S,#7;R1;r9,r10,r11,r12,"
     "r13,r14,r15,r16,r17,r18])]",
     0, ""},
    {"cyano_phosphonate", "P(O[A,a])(O[A,a])(=O)C#N", 0, ""},
    {"cyanohydrin", "[C;X4](-[OH,NH1,NH2,SH])(-C#N)", 0, ""},
    {"diamino_sulfide", "[N,n]~[S;!R;D2]~[N,n]", 0, ""},
    {"diazo_carbonyl", "[$(N=N=C~C=O),$(N#N-C~C=O)]", 0, ""},
    {"diazonium", "a[N+]#N", 0, ""},
    {"dicarbonyl_sulfonamide",
     "[$(N(-C(=O))(-C(=O))(-S(=O))),$(n([#6](=O))([#6](=O))([#16](=O)))]", 0,
     ""},
    {"disulfide_acyclic", "[S;!R;X2]-[S;!R;X2]", 0, ""},
    {"disulfonyliminoquinone", "S(=O)(=O)N=C1C=CC(=NS(=O)(=O))C=C1", 0, ""},
    {"double_trouble_warhead", "NC(C[S;D1])C([N;H1]([O;D1]))=O", 0, ""},
    {"flavanoid", "O=C2CC(a3aaaaa3)Oa1aaaaa12", 0, ""},
    {"four_nitriles", "C#N.C#N.C#N.C#N", 0, ""},
    {"gte_10_carbon_sb_chain",
     "[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]", 0,
     ""},
    {"gte_2_N_quats", "[N,n;H0;+;!$(N~O);!$(n~O)].[N,n;H0;+;!$(N~O);!$(n~O)]",
     0, ""},
    {"gte_2_free_phos", "P([O;D1])=O.P([O;D1])=O", 0, ""},
    {"gte_2_sulfonic_acid", "[C,c]S(=O)(=O)[O;D1].[C,c]S(=O)(=O)[O;D1]", 0, ""},
    {"gte_3_COOH", "C(=O)[O;D1].C(=O)[O;D1].C(=O)[O;D1]", 0, ""},
    {"gte_3_iodine", "[#53].[#53].[#53]", 0, ""},
    {"gte_4_basic_N",
     "[N;!$(N(=[N,O,S,C]));!$(N(S(=O)(=O)));!$(N(C(F)(F)(F)));!$(N(C#N));!$(N("
     "C(=O)));!$(N(C(=S)));!$(N(C(=N)));!$(N(#C));!$(N-c)].[N;!$(N(=[N,O,S,C]))"
     ";!$(N(S(=O)(=O)));!$(N(C(F)(F)(F)));!$(N(C#N));!$(N(C(=O)));!$(N(C(=S)));"
     "!$(N(C(=N)));!$(N(#C));!$(N-c)].[N;!$(N(=[N,O,S,C]));!$(N(S(=O)(=O)));!$("
     "N(C(F)(F)(F)));!$(N(C#N));!$(N(C(=O)));!$(N(C(=S)));!$(N(C(=N)));!$(N(#C)"
     ");!$(N-c)].[N;!$(N(=[N,O,S,C]));!$(N(S(=O)(=O)));!$(N(C(F)(F)(F)));!$(N("
     "C#N));!$(N(C(=O)));!$(N(C(=S)));!$(N(C(=N)));!$(N(#C));!$(N-c)]",
     0, ""},
    {"gte_4_nitro",
     "[$([N+](=O)[O-]),$(N(=O)=O)].[$([N+](=O)[O-]),$(N(=O)=O)].[$([N+](=O)[O-]"
     "),$(N(=O)=O)].[$([N+](=O)[O-]),$(N(=O)=O)]",
     0, ""},
    {"gte_5_phenolic_OH", "a[O;D1].a[O;D1].a[O;D1].a[O;D1].a[O;D1]", 0, ""},
    {"gte_7_aliphatic_OH",
     "C[O;D1].C[O;D1].C[O;D1].C[O;D1].C[O;D1].C[O;D1].C[O;D1]", 0, ""},
    {"gte_7_total_hal",
     "[Cl,Br,I].[Cl,Br,I].[Cl,Br,I].[Cl,Br,I].[Cl,Br,I].[Cl,Br,I].[Cl,Br,I]", 0,
     ""},
    {"gte_8_CF2_or_CH2",
     "[CH2,$(C(F)(F));R0][CH2,$(C(F)(F));R0][CH2,$(C(F)(F));R0][CH2,$(C(F)(F));"
     "R0][CH2,$(C(F)(F));R0][CH2,$(C(F)(F));R0][CH2,$(C(F)(F));R0][CH2,$(C(F)("
     "F));R0]",
     0, "Edited"},
    {"halo_5heterocycle_bis_EWG",
     "[#7,#8,#16]1[#6]([$(S(=O)(=O)),$([F,Cl]),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)"
     "),$([N+](=O)[O-]),$(C(=O))])[#6]([$(S(=O)(=O)),$([F,Cl]),$(C(F)(F)(F)),$("
     "C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O))])[#7][#6]1([Cl,Br,I])",
     0, ""},
    {"halo_acrylate",
     "[$([C;H2]),$([C&H1;$(C-F)]),$([C&H1;$(C-Cl)]),$([C&H1;$(C-Br)]),$([C&H1;$"
     "(C-I)]),$(C(F)F),$(C(Cl)Cl),$(C(Br)Br),$(C(I)I),$(C(F)Cl),$(C(F)Br),$(C("
     "F)I),$(C(Cl)Br),$(C(Br)I)](=[$([C&H1;$(C(-C(=O)))]),$(C(F)(C(=O))),$(C("
     "Cl)(C(=O))),$(C(Br)(C(=O))),$(C(I)(C(=O))),$(C(C)(C(=O))),$(C(c)(C(=O)))]"
     ")",
     0, ""},
    {"halo_imino", "C(=[#7])([Cl,Br,I,$(O(S(=O)(=O)))])", 0, ""},
    {"halo_olefin_bis_EWG",
     "C([Cl,Br,I,$(O(S(=O)(=O)))])=C([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)"
     "(=O)),$([N+](=O)[O-]),$(C=O)])([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)"
     "(=O)),$([N+](=O)[O-]),$(C=O)])",
     0, ""},
    {"halo_phenolic_carbonyl",
     "C(=O)Oc1c([Cl,F])[cH1,$(c[F,Cl])]c([F,Cl])[cH1,$(c[F,Cl])]c1([F,Cl])", 0,
     ""},
    {"halo_phenolic_sulfonyl",
     "S(=O)Oc1c([Cl,F])[cH1,$(c[F,Cl])]c([F,Cl])[cH1,$(c[F,Cl])]c1([F,Cl])", 0,
     ""},
    {"halogen_heteroatom", "[!C;!c;!H][F,Cl,Br,I]", 0, ""},
    {"hetero_silyl", "[Si]~[!#6]", 0, ""},
    {"hydrazine",
     "[N;X3;!$(N-S(=O)(=O));!$(N-C(F)(F)(F));!$(N-C#N);!$(N-C(=O));!$(N-C(=S));"
     "!$(N-C(=N))]-[N;X3;!$(N-S(=O)(=O));!$(N-C(F)(F)(F));!$(N-C#N);!$(N-C(=O))"
     ";!$(N-C(=S));!$(N-C(=N))]",
     0, ""},
    {"hydrazothiourea", "[N;!R]=NC(=S)N", 0, ""},
    {"hydroxamate_warhead", "C([N;H1]([O;D1]))=O", 0, ""},
    {"hyperval_sulfur", "[$([#16&D3]),$([#16&D4])]=,:[#6]", 0, ""},
    {"isonitrile", "[N+]#[C-]", 0, ""},
    {"keto_def_heterocycle",
     "[$(c([C;!R;!$(C-[N,O,S]);!$(C-[H])](=O))1naaaa1),$(c([C;!R;!$(C-[N,O,S]);"
     "!$(C-[H])](=O))1naa[n,s,o]1)]",
     0, ""},
    {"linear_polycyclic_aromatic_I",
     "[$(a12aaaaa1aa3a(aa(aaaa4)a4a3)a2),$(a12aaaaa1aa3a(aaa4a3aaaa4)a2),$("
     "a12aaaaa1a(aa5)a3a(aaa4a3a5aaa4)a2)]",
     0, ""},
    {"linear_polycyclic_aromatic_II",
     "[$(a12aaaa4a1a3a(aaaa3aa4)aa2),$(a12aaaaa1a3a(aaa4a3aaaa4)aa2),$(a1(a("
     "aaaa4)a4a3a2aaaa3)a2aaaa1)]",
     0, ""},
    {"maleimide_etc",
     "[$([C;H1]),$(C(-[F,Cl,Br,I]))]1=[$([C;H1]),$(C(-[F,Cl,Br,I]))]C(=O)[N,O,"
     "S]C(=O)1",
     0, ""},
    {"meldrums_acid_deriv", "O=C1OC(C)(C)OC(C1)=O", 0, ""},
    {"monofluoroacetate", "[C;H2](F)C(=O)[O,N,S]", 0, ""},
    {"nitrone", "[C;!R]=[N+][O;D1]", 0, ""},
    {"nitrosamine", "N-[N;X2](=O)", 0, ""},
    {"non_ring_CH2O_acetal", "[O,N,S;!$(S~O)]!@[CH2]!@[O,S,N;!$(S~O)]", 0, ""},
    {"non_ring_acetal", "[O,N,S;!$(S~O)]!@[C;H1;X4]!@[O,N,S;!$(S~O)]", 0, ""},
    {"non_ring_ketal", "[O,N,S;!$(S~O)]!@[C;H0;X4](!@[O,N,S;!$(S~O)])(C)", 0,
     ""},
    {"ortho_hydroiminoquinone", "c1c([N;D1])c([N;D1])c[cH1][cH1]1", 0, ""},
    {"ortho_hydroquinone", "a1c([O,S;D1])c([O,S;D1])a[cH1][cH1]1", 0, ""},
    {"ortho_nitrophenyl_carbonyl",
     "[#6]1(-O-[C;!R](=[O,N;!R]))[#6]([$(N(=O)(=O)),$([N+](=O)[O-])])[#6][#6][#"
     "6][#6]1",
     0, ""},
    {"ortho_quinone",
     "[CH1,$(C(-[Cl,Br,I]))]1=CC(=[O,N,S;!R])C(=[O,N,S])C=[CH1,$(C(-[Cl,Br,I]))"
     "]1",
     0, ""},
    {"oxaziridine", "C1~[O,S]~N1", 0, ""},
    {"oxime", "[$(C=N[O;D1]);!$(C=[N+])][#6][#6]", 0, ""},
    {"oxonium", "[o+,O+]", 0, ""},
    {"para_hydroiminoquinone", "a1[cH1]c([N;D1])[cH1]ac([N;D1])1", 0, ""},
    {"para_hydroquinone", "a1[cH1]c([O,S;D1])[cH1]ac([O,S;D1])1", 0, ""},
    {"para_nitrophenyl_ester",
     "[#6]1(-O(-[C;!R](-[!N])(=[O,N;!R])))[#6][#6][#6]([$(N(=O)(=O)),$([N+](=O)"
     "[O-])])[#6][#6]1",
     0, ""},
    {"para_quinone",
     "[CH1,$(C(-[Cl,Br,I]))]1=[CH1,$(C(-[Cl,Br,I]))]C(=[O,N,S])[CH1,$(C(-[Cl,"
     "Br,I]))]=[CH1,$(C(-[Cl,Br,I]))]C1(=[O,N,S])",
     0, ""},
    {"paraquat_like",
     "[#6]1[#6][#6]([#6]2[#6][#6][#7;+][#6][#6]2)[#6][#6][#7;+]1", 0, ""},
    {"pentafluorophenylester", "C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)", 0, ""},
    {"perchloro_cp", "C1(Cl)(Cl)C(Cl)C(Cl)=C(Cl)C1(Cl)", 0, ""},
    {"perhalo_dicarbonyl_phenyl",
     "c1(C=O)c([Br,Cl,I])c([Br,Cl,I])c([Br,Cl,I])c([Br,Cl,I])c1(C=O)", 0, ""},
    {"perhalo_phenyl",
     "c1c([Br,Cl,I])c([Br,Cl,I])c([Br,Cl,I])c([Br,Cl,I])c1([Br,Cl,I])", 0, ""},
    {"peroxide", "[#8]~[#8]", 0, ""},
    {"phenolate_bis_EWG",
     "O=[C,S]Oc1aaa([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)["
     "O-]),$(C(=O)O),$(C(=O)N)])aa([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(="
     "O)),$([N+](=O)[O-]),$(C(=O)O),$(C(=O)N)])1",
     0, ""},
    {"phos_serine_warhead", "NC(COP(O)(O)=O)C(O)=O", 0, ""},
    {"phos_threonine_warhead", "NC(C(C)OP(O)(O)=O)C(O)=O", 0, ""},
    {"phos_tyrosine_warhead", "NC(Cc1ccc(OP(O)(O)=O)cc1)C(O)=O", 0, ""},
    {"phosphite", "[c,C]-[P;v3]", 0, ""},
    {"phosphonium", "[#15;+]~[!O]", 0, ""},
    {"phosphorane", "C=P", 0, ""},
    {"phosphorous_nitrogen_bond", "[#15]~[N,n]", 0, ""},
    {"phosphorus_phosphorus_bond", "P~P", 0, ""},
    {"phosphorus_sulfur_bond", "P~S", 0, ""},
    {"polyene", "C=[C;!R][C;!R]=[C;!R][C;!R]=[C;!R]", 0, ""},
    {"polyhalo_phenol_a",
     "c1c([O;D1])c(-[Cl,Br,I])c(-[Cl,Br,I])cc1.c1c([O;D1])c(-[Cl,Br,I])c(-[Cl,"
     "Br,I])cc1",
     0, ""},
    {"polyhalo_phenol_b",
     "c1c([O;D1])c(-[Cl,Br,I])cc(-[Cl,Br,I])c1.c1c([O;D1])c(-[Cl,Br,I])cc(-[Cl,"
     "Br,I])c1",
     0, ""},
    {"polyhalo_phenol_c",
     "c1c([O;D1])ccc(-[Cl,Br,I])c(-[Cl,Br,I])1.c1c([O;D1])ccc(-[Cl,Br,I])c(-["
     "Cl,Br,I])1",
     0, ""},
    {"polyhalo_phenol_d",
     "c(-[Cl,Br,I])1c([O;D1])c(-[Cl,Br,I])ccc1.c(-[Cl,Br,I])1c([O;D1])c(-[Cl,"
     "Br,I])ccc1",
     0, ""},
    {"polyhalo_phenol_e",
     "c1c([O;D1])ccc(-[Cl,Br,I])c(-[Cl,Br,I])1.c1c([O;D1])ccc(-[Cl,Br,I])c(-["
     "Cl,Br,I])1",
     0, ""},
    {"polysulfide", "[S;D2]-[S;D2]-[S;D2]", 0, ""},
    {"porphyrin", "[#6;r16,r17,r18]~[#6]1~[#6]~[#6]~[#6](~[#6])~[#7]1", 0, ""},
    {"primary_halide_sulfate",
     "[CH2][Cl,Br,I,$(O(S(=O)(=O)[!$(N);!$([O&D1])]))]", 0, ""},
    {"quat_N_N", "[N,n;R;+]!@[N,n]", 0, ""},
    {"quat_N_acyl", "[N,n;+]!@C(=O)", 0, ""},
    {"quinone_methide",
     "[#6;!$([#6](-[N,O,S]))]1=[#6;!$([#6](-[N,O,S]))][#6](=[#6])[#6;!$([#6](-["
     "N,O,S]))]=[#6;!$([#6](-[N,O,S]))][#6]1(=[O,N,S])",
     0, ""},
    {"rhodanine", "C(=C)1SC(=S)NC(=O)1", 0, ""},
    {"secondary_halide_sulfate",
     "[CH;!$(C=C)][Cl,Br,I,$(O(S(=O)(=O)[!$(N);!$([O&D1])]))]", 0, ""},
    {"sulf_D2_nitrogen",
     "[S;D2](-[N;!$(N(=C));!$(N(-S(=O)(=O)));!$(N(-C(=O)))])", 0, ""},
    {"sulf_D2_oxygen_D2", "[S;D2][O;D2]", 0, ""},
    {"sulf_D3_nitrogen", "[S;D3](-N)(-[c,C])(-[c,C])", 0, ""},
    {"sulfite_sulfate_ester", "[C,c]OS(=O)O[C,c]", 0, ""},
    {"sulfonium", "[S+;X3;$(S-C);!$(S-[O;D1])]", 0, ""},
    {"sulfonyl_anhydride", "[$(C(=O)),$(S(=O)(=O))][O,S](S(=O)(=O))", 0, ""},
    {"sulfonyl_halide", "S(=O)(=O)[F,Cl,Br,I]", 0, ""},
    {"sulfonyl_heteroatom", "[!#6;!#1;!#11;!#19]O(S(=O)(=O)(-[C,c]))", 0, ""},
    {"sulphonyl_cyanide", "S(=O)(=O)C#N", 0, ""},
    {"tertiary_halide_sulfate",
     "[C;X4](-[Cl,Br,I,$(O(S(=O)(=O)[!$(N);!$([O&D1])]))])(-[c,C])(-[c,C])(-[c,"
     "C])",
     0, ""},
    {"thio_hydroxamate", "[S;D2]([$(N(=C)),$(N(-S(=O)(=O))),$(N(-C(=O)))])", 0,
     ""},
    {"thio_xanthate", "[S;!R]-[C;!R](=[S;!R])(-[S;!R])", 0, ""},
    {"thiocarbonate", "SC(=O)[O,S]", 0, ""},
    {"thioester", "[S;!R;H0]C(=[S,O;!R])([!O;!S;!N])", 0, ""},
    {"thiol_warhead", "NC(C[S;D1])C(O)=O", 0, ""},
    {"thiopyrylium", "c1[S,s;+]cccc1", 0, ""},
    {"thiosulfoxide", "[C,c][S;X3](~O)-S", 0, ""},
    {"triamide",
     "[$(N(-C(=O))(-C(=O))(-C(=O))),$(n([#6](=O))([#6](=O))([#6](=O)))]", 0,
     ""},
    {"triaryl_phosphine_oxide", "P(=O)(a)(a)(a)", 0, ""},
    {"trichloromethyl_ketone",
     "[$(C(=O));!$(C-N);!$(C-O);!$(C-S)]C(Cl)(Cl)(Cl)", 0, ""},
    {"triflate", "OS(=O)(=O)(C(F)(F)(F))", 0, ""},
    {"trifluoroacetate_ester", "C(F)(F)(F)C(=O)O", 0, ""},
    {"trifluoroacetate_thioester", "C(F)(F)(F)C(=O)S", 0, ""},
    {"trifluoromethyl_ketone", "[$(C(=O));!$(C-N);!$(C-O);!$(C-S)]C(F)(F)(F)",
     0, ""},
    {"trihalovinyl_heteroatom",
     "C(-[Cl,Br,I])(-[Cl,Br,I])=C(-[Cl,Br,I])(-[N,O,S])", 0, ""},
    {"trinitro_aromatic",
     "[$(a1aaa([$(N(=O)(=O)),$([N+](=O)[O-])])a([$(N(=O)(=O)),$([N+](=O)[O-])])"
     "a1([$(N(=O)(=O)),$([N+](=O)[O-])])),$(a1aa([$(N(=O)(=O)),$([N+](=O)[O-])]"
     ")a([$(N(=O)(=O)),$([N+](=O)[O-])])aa1([$(N(=O)(=O)),$([N+](=O)[O-])])),$("
     "a1a([$(N(=O)(=O)),$([N+](=O)[O-])])aa([$(N(=O)(=O)),$([N+](=O)[O-])])aa1("
     "[$(N(=O)(=O)),$([N+](=O)[O-])]))]",
     0, ""},
    {"trinitromethane_derivative",
     "C([$([N+](=O)[O-]),$(N(=O)=O)])([$([N+](=O)[O-]),$(N(=O)=O)])([$([N+](=O)"
     "[O-]),$(N(=O)=O)])",
     0, ""},
    {"tris_activated_aryl_ester",
     "[$(O=[C,S]Oc1a([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](="
     "O)[O-]),$(C(=O)O),$(C(=O)N)])a([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(="
     "O)(=O)),$([N+](=O)[O-]),$(C(=O)O),$(C(=O)N)])a([$(S(=O)(=O)),F,$(C(F)(F)("
     "F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O)O),$(C(=O)N)])aa1),$(O=[C,"
     "S]Oc1a([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),"
     "$(C(=O)O),$(C(=O)N)])a([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),"
     "$([N+](=O)[O-]),$(C(=O)O),$(C(=O)N)])aaa([$(S(=O)(=O)),F,$(C(F)(F)(F)),$("
     "C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O)O),$(C(=O)N)])1),$(O=[C,S]Oc1a(["
     "$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O)"
     "O),$(C(=O)N)])aa([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+]"
     "(=O)[O-]),$(C(=O)O),$(C(=O)N)])a([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$("
     "N(=O)(=O)),$([N+](=O)[O-]),$(C(=O)O),$(C(=O)N)])a1),$(O=[C,S]Oc1a([$(S(="
     "O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O)O),$("
     "C(=O)N)])aa([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)["
     "O-]),$(C(=O)O),$(C(=O)N)])aa([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),$(N(=O)"
     "(=O)),$([N+](=O)[O-]),$(C(=O)O),$(C(=O)N)])1)]",
     0, ""},
    {"trisub_bis_act_olefin",
     "[CH;!R;!$(C-N)]=C([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+]("
     "=O)[O-]),$(C(=O))])([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+"
     "](=O)[O-]),$(C(=O))])",
     0, "Edited"},
    {"vinyl_carbonyl_EWG",
     "[C;!R]([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$("
     "C=O)])([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$("
     "C=O)])=[C;!R]([C;!R](=O))([!$([#8]);!$([#7])])",
     0, ""}};
const unsigned int NUM_NIH =
    static_cast<unsigned int>(sizeof(NIH) / sizeof(FilterData_t));

const FilterProperty_t NIH_PROPS[] = {
    {"FilterSet", "NIH"},
    {"Scope", "annotate compounds with problematic functional groups"},
    {"Reference",
     "Doveston R, et al. A Unified Lead-oriented Synthesis of over Fifty "
     "Molecular Scaffolds. Org Biomol Chem 13 (2014) 859D65. "
     "doi:10.1039/C4OB02287D."},
    {"Reference",
     "Jadhav A, et al. Quantitative Analyses of Aggregation, Autofluorescence, "
     "and Reactivity Artifacts in a Screen for Inhibitors of a Thiol Protease. "
     "J Med Chem 53 (2009) 37D51. doi:10.1021/jm901070c."}};
const unsigned int NUM_NIH_PROPS =
    static_cast<unsigned int>(sizeof(NIH_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// PAINS_A data
// # Reference: Baell JB, Holloway GA. New Substructure Filters for Removal of
// Pan Assay Interference Compounds (PAINS) from Screening Libraries and for
// Their Exclusion in Bioassays. J Med Chem 53 (2010) 2719D40.
// doi:10.1021/jm901137j.
// # Scope: PAINS filters (family A)
//
#include "pains_a.in"

const unsigned int NUM_PAINS_A =
    static_cast<unsigned int>(sizeof(PAINS_A) / sizeof(FilterData_t));

const FilterProperty_t PAINS_A_PROPS[] = {
    {"FilterSet", "PAINS_A"},
    {"Reference",
     "Baell JB, Holloway GA. New Substructure Filters for Removal of Pan Assay "
     "Interference Compounds (PAINS) from Screening Libraries and for Their "
     "Exclusion in Bioassays. J Med Chem 53 (2010) 2719D40. "
     "doi:10.1021/jm901137j."},
    {"Scope", "PAINS filters (family A)"}};
const unsigned int NUM_PAINS_A_PROPS =
    static_cast<unsigned int>(sizeof(PAINS_A_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// PAINS_B data
// # Reference: Baell JB, Holloway GA. New Substructure Filters for Removal of
// Pan Assay Interference Compounds (PAINS) from Screening Libraries and for
// Their Exclusion in Bioassays. J Med Chem 53 (2010) 2719D40.
// doi:10.1021/jm901137j.
// # Scope: PAINS filters (family B)
// # sulfonamide_B(41) c:1:c:c(:c:c:c:1-[#8]-[#1])-[#7](-[#1])-[#16](=[#8])=[#8]
// 0
// # sulfonamide_B(41) [N;H1](c1ccc([O;H1])cc1)S(=O)=O 0
// # imidazole_A(19)
// n:1:c(:n(:c(:c:1-c:2:c:c:c:c:c:2)-c:3:c:c:c:c:c:3)-[#1])-[#6]:,=[!#1] 0
//
#include "pains_b.in"

const unsigned int NUM_PAINS_B =
    static_cast<unsigned int>(sizeof(PAINS_B) / sizeof(FilterData_t));

const FilterProperty_t PAINS_B_PROPS[] = {
    {"FilterSet", "PAINS_B"},
    {"Reference",
     "Baell JB, Holloway GA. New Substructure Filters for Removal of Pan Assay "
     "Interference Compounds (PAINS) from Screening Libraries and for Their "
     "Exclusion in Bioassays. J Med Chem 53 (2010) 2719D40. "
     "doi:10.1021/jm901137j."},
    {"Scope", "PAINS filters (family B)"},
};
const unsigned int NUM_PAINS_B_PROPS =
    static_cast<unsigned int>(sizeof(PAINS_B_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// PAINS_C data
// # Reference: Baell JB, Holloway GA. New Substructure Filters for Removal of
// Pan Assay Interference Compounds (PAINS) from Screening Libraries and for
// Their Exclusion in Bioassays. J Med Chem 53 (2010) 2719D40.
// doi:10.1021/jm901137j.
// # Scope: PAINS filters (family C)
//
#include "pains_c.in"

const unsigned int NUM_PAINS_C =
    static_cast<unsigned int>(sizeof(PAINS_C) / sizeof(FilterData_t));

const FilterProperty_t PAINS_C_PROPS[] = {
    {"FilterSet", "PAINS_C"},
    {"Reference",
     "Baell JB, Holloway GA. New Substructure Filters for Removal of Pan Assay "
     "Interference Compounds (PAINS) from Screening Libraries and for Their "
     "Exclusion in Bioassays. J Med Chem 53 (2010) 2719D40. "
     "doi:10.1021/jm901137j."},
    {"Scope", "PAINS filters (family C)"}};
const unsigned int NUM_PAINS_C_PROPS =
    static_cast<unsigned int>(sizeof(PAINS_C_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// ZINC data
// # Reference: http://blaster.docking.org/filtering/
// # Scope: drug-likeness and unwanted functional group filters
//

const FilterData_t ZINC[] = {
    {"Non-Hydrogen_atoms", "[a,A]", 40, ""},
    {"carbons", "[#6]", 40, ""},
    {"N,O,S", "[#7,#8,#16]", 20, ""},
    {"Sulfonyl_halides", "S(=O)(=O)[Cl,Br]", 1, ""},
    {"Acid_halides", "[S,C](=[O,S])[F,Br,Cl,I]", 1, ""},
    {"Alkyl_halides", "[Br,Cl,I][CX4;CH,CH2]", 1, ""},
    {"Phosphenes", "cPc", 0, ""},
    {"Heptanes", "[CD1][CD2][CD2][CD2][CD2][CD2][CD2]", 0, ""},
    {"Perchlorates", "OCl(O)(O)(O)", 0, ""},
    {"Fluorines", "F", 7, ""},
    {"Cl,Br,I", "[Cl,Br,I]", 6, ""},
    {"Carbazides", "O=CN=[N+]=[N-]", 0, ""},
    {"Acid_anhydrides", "C(=O)OC(=O)", 0, ""},
    {"Peroxides", "OO", 0, ""},
    {"Iso(thio)cyanates", "N=C=[S,O]", 1, ""},
    {"Thiocyanates", "SC#N", 1, ""},
    {"Phosphoranes", "C=P", 0, ""},
    {"P/S_halides", "[P,S][Cl,Br,F,I]", 0, ""},
    {"Cyanohydrines", "N#CC[OH]", 0, ""},
    {"Carbazides", "O=CN=[N+]=[N-]", 0, ""},
    {"Sulfate_esters", "COS(=O)O[C,c]", 1, ""},
    {"Sulfonates", "COS(=O)(=O)[C,c]", 1, ""},
    {"Pentafluorophenyl_esters", "C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)", 0, ""},
    {"Paranitrophenyl_esters", "C(=O)Oc1ccc(N(=O)=O)cc1", 0, ""},
    {"HOBt_esters", "C(=O)Onnn", 0, ""},
    {"Triflates", "OS(=O)(=O)C(F)(F)F", 0, ""},
    {"Lawesson's_reagents", "P(=S)(S)S", 0, ""},
    {"Phosphoramides", "NP(=O)(N)N", 0, ""},
    {"Aromatic_azides", "cN=[N+]=[N-]", 0, ""},
    {"Quaternary_C,Cl,I,P,S", "[C+,Cl+,I+,P+,S+]", 2, ""},
    {"Beta_carbonyl_quaternary_N", "C(=O)C[N+,n+]", 2, ""},
    {"Acylhydrazides", "[N;R0][N;R0]C(=O)", 2, ""},
    {"Chloramidines", "[Cl]C([C&R0])=N", 0, ""},
    {"Isonitriles", "[N+]#[C-]", 0, ""},
    {"Triacyloximes", "C(=O)N(C(=O))OC(=O)", 0, ""},
    {"Acyl_cyanides", "N#CC(=O)", 0, ""},
    {"Sulfonyl_cyanides", "S(=O)(=O)C#N", 0, ""},
    {"Cyanophosphonates", "P(OCC)(OCC)(=O)C#N", 0, ""},
    {"Azocyanamides", "[N;R0]=[N;R0]C#N", 0, ""},
    {"Azoalkanals", "[N;R0]=[N;R0]CC=O", 0, ""},
    {"(Thio)epoxides,aziridines", "C1[O,S,N]C1", 2, ""},
    {"Benzylic_quaternary_N", "cC[N+]", 2, ""},
    {"Thioesters", "C[O,S;R0][C;R0](=S)", 2, ""},
    {"Diand_Triphosphates", "P(=O)([OH])OP(=O)[OH]", 3, ""},
    {"Aminooxy(oxo)", "[#7]O[#6,#16]=O", 2, ""},
    {"nitros", "N(~[OD1])~[OD1]", 2, ""},
    {"Imines", "C=[N;R0]*", 2, ""},
    {"Acrylonitriles", "N#CC=C", 2, ""},
    {"Propenals", "C=CC(=O)[!#7;!#8]", 2, ""},
    {"Quaternary_N", "[ND4+]", 1, ""}};
const unsigned int NUM_ZINC =
    static_cast<unsigned int>(sizeof(ZINC) / sizeof(FilterData_t));

const FilterProperty_t ZINC_PROPS[] = {
    {"FilterSet", "ZINC"},
    {"Reference", "http://blaster.docking.org/filtering/"},
    {"Scope", "drug-likeness and unwanted functional group filters"}};
const unsigned int NUM_ZINC_PROPS =
    static_cast<unsigned int>(sizeof(ZINC_PROPS) / sizeof(FilterProperty_t));

////////////////////////////////////////////////////////////////////////
// API
unsigned int GetNumEntries(FilterCatalogParams::FilterCatalogs catalog) {
  switch (catalog) {
    case FilterCatalogParams::BRENK:
      return NUM_BRENK;
    case FilterCatalogParams::NIH:
      return NUM_NIH;
    case FilterCatalogParams::PAINS_A:
      return NUM_PAINS_A;
    case FilterCatalogParams::PAINS_B:
      return NUM_PAINS_B;
    case FilterCatalogParams::PAINS_C:
      return NUM_PAINS_C;
    case FilterCatalogParams::ZINC:
      return NUM_ZINC;
    default:
      return 0;
  }
}

const FilterData_t* GetFilterData(FilterCatalogParams::FilterCatalogs catalog) {
  switch (catalog) {
    case FilterCatalogParams::BRENK:
      return BRENK;
    case FilterCatalogParams::NIH:
      return NIH;
    case FilterCatalogParams::PAINS_A:
      return PAINS_A;
    case FilterCatalogParams::PAINS_B:
      return PAINS_B;
    case FilterCatalogParams::PAINS_C:
      return PAINS_C;
    case FilterCatalogParams::ZINC:
      return ZINC;
    default:
      return nullptr;
  }
}

unsigned GetNumPropertyEntries(FilterCatalogParams::FilterCatalogs catalog) {
  switch (catalog) {
    case FilterCatalogParams::BRENK:
      return NUM_BRENK_PROPS;
    case FilterCatalogParams::NIH:
      return NUM_NIH_PROPS;
    case FilterCatalogParams::PAINS_A:
      return NUM_PAINS_A_PROPS;
    case FilterCatalogParams::PAINS_B:
      return NUM_PAINS_B_PROPS;
    case FilterCatalogParams::PAINS_C:
      return NUM_PAINS_C_PROPS;
    case FilterCatalogParams::ZINC:
      return NUM_ZINC_PROPS;
    default:
      return 0;
  }
}

const FilterProperty_t* GetFilterProperties(
    FilterCatalogParams::FilterCatalogs catalog) {
  switch (catalog) {
    case FilterCatalogParams::BRENK:
      return BRENK_PROPS;
    case FilterCatalogParams::NIH:
      return NIH_PROPS;
    case FilterCatalogParams::PAINS_A:
      return PAINS_A_PROPS;
    case FilterCatalogParams::PAINS_B:
      return PAINS_B_PROPS;
    case FilterCatalogParams::PAINS_C:
      return PAINS_C_PROPS;
    case FilterCatalogParams::ZINC:
      return ZINC_PROPS;
    default:
      return nullptr;
  }
}
}  // namespace RDKit
