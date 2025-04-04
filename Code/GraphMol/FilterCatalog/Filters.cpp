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
#include "nih.in"

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

/////////////////////////////////////////////////////////////////////////////////////////
// # Chembl 23 structural alert, rule set CHEMBL_Glaxo
//
#include "chembl_glaxo.in"

const unsigned int NUM_CHEMBL_Glaxo =
    static_cast<unsigned int>(sizeof(CHEMBL_Glaxo) / sizeof(FilterData_t));

const FilterProperty_t CHEMBL_Glaxo_PROPS[] = {
    {"FilterSet", "ChEMBL23_Glaxo"},
    {"Reference", "https://github.com/PatWalters/rd_filters"},
    {"Scope", "ChEMBL 23 structural alerts rule-set: CHEMBL_Glaxo"}};

const unsigned int NUM_CHEMBL_Glaxo_PROPS = static_cast<unsigned int>(
    sizeof(CHEMBL_Glaxo_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// # Chembl 23 structural alert, rule set CHEMBL_Dundee
//
#include "chembl_dundee.in"

const unsigned int NUM_CHEMBL_Dundee =
    static_cast<unsigned int>(sizeof(CHEMBL_Dundee) / sizeof(FilterData_t));

const FilterProperty_t CHEMBL_Dundee_PROPS[] = {
    {"FilterSet", "ChEMBL23_Dundee"},
    {"Reference", "https://github.com/PatWalters/rd_filters"},
    {"Scope", "ChEMBL 23 structural alerts rule-set: CHEMBL_Dundee"}};

const unsigned int NUM_CHEMBL_Dundee_PROPS = static_cast<unsigned int>(
    sizeof(CHEMBL_Dundee_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// # Chembl 23 structural alert, rule set CHEMBL_BMS
//
#include "chembl_bms.in"

const unsigned int NUM_CHEMBL_BMS =
    static_cast<unsigned int>(sizeof(CHEMBL_BMS) / sizeof(FilterData_t));

const FilterProperty_t CHEMBL_BMS_PROPS[] = {
    {"FilterSet", "ChEMBL23_BMS"},
    {"Reference", "https://github.com/PatWalters/rd_filters"},
    {"Scope", "ChEMBL 23 structural alerts rule-set: CHEMBL_BMS"}};

const unsigned int NUM_CHEMBL_BMS_PROPS = static_cast<unsigned int>(
    sizeof(CHEMBL_BMS_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// # Chembl 23 structural alert, rule set CHEMBL_SureChEMBL
//
#include "chembl_surechembl.in"

const unsigned int NUM_CHEMBL_SureChEMBL =
    static_cast<unsigned int>(sizeof(CHEMBL_SureChEMBL) / sizeof(FilterData_t));

const FilterProperty_t CHEMBL_SureChEMBL_PROPS[] = {
    {"FilterSet", "ChEMBL23_SureChEMBL"},
    {"Reference", "https://github.com/PatWalters/rd_filters"},
    {"Scope", "ChEMBL 23 structural alerts rule-set: CHEMBL_SureChEMBL"}};

const unsigned int NUM_CHEMBL_SureChEMBL_PROPS = static_cast<unsigned int>(
    sizeof(CHEMBL_SureChEMBL_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// # Chembl 23 structural alert, rule set CHEMBL_MLSMR
//
#include "chembl_mlsmr.in"

const unsigned int NUM_CHEMBL_MLSMR =
    static_cast<unsigned int>(sizeof(CHEMBL_MLSMR) / sizeof(FilterData_t));

const FilterProperty_t CHEMBL_MLSMR_PROPS[] = {
    {"FilterSet", "ChEMBL23_MLSMR"},
    {"Reference", "https://github.com/PatWalters/rd_filters"},
    {"Scope", "ChEMBL 23 structural alerts rule-set: CHEMBL_MLSMR"}};

const unsigned int NUM_CHEMBL_MLSMR_PROPS = static_cast<unsigned int>(
    sizeof(CHEMBL_MLSMR_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// # Chembl 23 structural alert, rule set CHEMBL_Inpharmatica
//
#include "chembl_inpharmatica.in"

const unsigned int NUM_CHEMBL_Inpharmatica = static_cast<unsigned int>(
    sizeof(CHEMBL_Inpharmatica) / sizeof(FilterData_t));

const FilterProperty_t CHEMBL_Inpharmatica_PROPS[] = {
    {"FilterSet", "ChEMBL23_Inpharmatica"},
    {"Reference", "https://github.com/PatWalters/rd_filters"},
    {"Scope", "ChEMBL 23 structural alerts rule-set: CHEMBL_Inpharmatica"}};

const unsigned int NUM_CHEMBL_Inpharmatica_PROPS = static_cast<unsigned int>(
    sizeof(CHEMBL_Inpharmatica_PROPS) / sizeof(FilterProperty_t));

/////////////////////////////////////////////////////////////////////////////////////////
// # Chembl 23 structural alert, rule set CHEMBL_LINT
//
#include "chembl_lint.in"

const unsigned int NUM_CHEMBL_LINT =
    static_cast<unsigned int>(sizeof(CHEMBL_LINT) / sizeof(FilterData_t));

const FilterProperty_t CHEMBL_LINT_PROPS[] = {
    {"FilterSet", "ChEMBL23_LINT"},
    {"Reference", "https://github.com/PatWalters/rd_filters"},
    {"Scope", "ChEMBL 23 structural alerts rule-set: CHEMBL_LINT"}};

const unsigned int NUM_CHEMBL_LINT_PROPS = static_cast<unsigned int>(
    sizeof(CHEMBL_LINT_PROPS) / sizeof(FilterProperty_t));

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
    case FilterCatalogParams::CHEMBL_Glaxo:
      return NUM_CHEMBL_Glaxo;
    case FilterCatalogParams::CHEMBL_Dundee:
      return NUM_CHEMBL_Dundee;
    case FilterCatalogParams::CHEMBL_BMS:
      return NUM_CHEMBL_BMS;
    case FilterCatalogParams::CHEMBL_SureChEMBL:
      return NUM_CHEMBL_SureChEMBL;
    case FilterCatalogParams::CHEMBL_MLSMR:
      return NUM_CHEMBL_MLSMR;
    case FilterCatalogParams::CHEMBL_Inpharmatica:
      return NUM_CHEMBL_Inpharmatica;
    case FilterCatalogParams::CHEMBL_LINT:
      return NUM_CHEMBL_LINT;
    default:
      return 0;
  }
}

const FilterData_t *GetFilterData(FilterCatalogParams::FilterCatalogs catalog) {
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
    case FilterCatalogParams::CHEMBL_Glaxo:
      return CHEMBL_Glaxo;
    case FilterCatalogParams::CHEMBL_Dundee:
      return CHEMBL_Dundee;
    case FilterCatalogParams::CHEMBL_BMS:
      return CHEMBL_BMS;
    case FilterCatalogParams::CHEMBL_SureChEMBL:
      return CHEMBL_SureChEMBL;
    case FilterCatalogParams::CHEMBL_MLSMR:
      return CHEMBL_MLSMR;
    case FilterCatalogParams::CHEMBL_Inpharmatica:
      return CHEMBL_Inpharmatica;
    case FilterCatalogParams::CHEMBL_LINT:
      return CHEMBL_LINT;
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
    case FilterCatalogParams::CHEMBL_Glaxo:
      return NUM_CHEMBL_Glaxo_PROPS;
    case FilterCatalogParams::CHEMBL_Dundee:
      return NUM_CHEMBL_Dundee_PROPS;
    case FilterCatalogParams::CHEMBL_BMS:
      return NUM_CHEMBL_BMS_PROPS;
    case FilterCatalogParams::CHEMBL_SureChEMBL:
      return NUM_CHEMBL_SureChEMBL_PROPS;
    case FilterCatalogParams::CHEMBL_MLSMR:
      return NUM_CHEMBL_MLSMR_PROPS;
    case FilterCatalogParams::CHEMBL_Inpharmatica:
      return NUM_CHEMBL_Inpharmatica_PROPS;
    case FilterCatalogParams::CHEMBL_LINT:
      return NUM_CHEMBL_LINT_PROPS;
    default:
      return 0;
  }
}

const FilterProperty_t *GetFilterProperties(
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
    case FilterCatalogParams::CHEMBL_Glaxo:
      return CHEMBL_Glaxo_PROPS;
    case FilterCatalogParams::CHEMBL_Dundee:
      return CHEMBL_Dundee_PROPS;
    case FilterCatalogParams::CHEMBL_BMS:
      return CHEMBL_BMS_PROPS;
    case FilterCatalogParams::CHEMBL_SureChEMBL:
      return CHEMBL_SureChEMBL_PROPS;
    case FilterCatalogParams::CHEMBL_MLSMR:
      return CHEMBL_MLSMR_PROPS;
    case FilterCatalogParams::CHEMBL_Inpharmatica:
      return CHEMBL_Inpharmatica_PROPS;
    case FilterCatalogParams::CHEMBL_LINT:
      return CHEMBL_LINT_PROPS;
    default:
      return nullptr;
  }
}
}  // namespace RDKit
