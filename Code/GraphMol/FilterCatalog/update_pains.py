# must be run from this directory

import csv
import os
import sys

py3 = sys.version_info[0] == 3

pains_a = [
  'ene_six_het_A(483)', 'hzone_phenol_A(479)', 'anil_di_alk_A(478)', 'indol_3yl_alk(461)',
  'quinone_A(370)', 'azo_A(324)', 'imine_one_A(321)', 'mannich_A(296)', 'anil_di_alk_B(251)',
  'anil_di_alk_C(246)', 'ene_rhod_A(235)', 'hzone_phenol_B(215)', 'ene_five_het_A(201)',
  'anil_di_alk_D(198)', 'imine_one_isatin(189)', 'anil_di_alk_E(186)'
]

pains_b = [
  'thiaz_ene_A(128)', 'pyrrole_A(118)', 'catechol_A(92)', 'ene_five_het_B(90)',
  'imine_one_fives(89)', 'ene_five_het_C(85)', 'hzone_pipzn(79)', 'keto_keto_beta_A(68)',
  'hzone_pyrrol(64)', 'ene_one_ene_A(57)', 'cyano_ene_amine_A(56)', 'ene_five_one_A(55)',
  'cyano_pyridone_A(54)', 'anil_alk_ene(51)', 'amino_acridine_A(46)', 'ene_five_het_D(46)',
  'thiophene_amino_Aa(45)', 'ene_five_het_E(44)', 'sulfonamide_A(43)', 'thio_ketone(43)',
  'sulfonamide_B(41)', 'anil_no_alk(40)', 'thiophene_amino_Ab(40)', 'het_pyridiniums_A(39)',
  'anthranil_one_A(38)', 'cyano_imine_A(37)', 'diazox_sulfon_A(36)', 'hzone_anil_di_alk(35)',
  'rhod_sat_A(33)', 'hzone_enamin(30)', 'pyrrole_B(29)', 'thiophene_hydroxy(28)',
  'cyano_pyridone_B(27)', 'imine_one_sixes(27)', 'dyes5A(27)', 'naphth_amino_A(25)',
  'naphth_amino_B(25)', 'ene_one_ester(24)', 'thio_dibenzo(23)', 'cyano_cyano_A(23)',
  'hzone_acyl_naphthol(22)', 'het_65_A(21)', 'imidazole_A(19)', 'ene_cyano_A(19)',
  'anthranil_acid_A(19)', 'dyes3A(19)', 'dhp_bis_amino_CN(19)', 'het_6_tetrazine(18)',
  'ene_one_hal(17)', 'cyano_imine_B(17)', 'thiaz_ene_B(17)', 'ene_rhod_B(16)',
  'thio_carbonate_A(15)', 'anil_di_alk_furan_A(15)', 'ene_five_het_F(15)'
]

pains_c = [
  'anil_di_alk_F(14)', 'hzone_anil(14)', 'het_5_pyrazole_OH(14)', 'het_thio_666_A(13)',
  'styrene_A(13)', 'ene_rhod_C(13)', 'dhp_amino_CN_A(13)', 'cyano_imine_C(12)', 'thio_urea_A(12)',
  'thiophene_amino_B(12)', 'keto_keto_beta_B(12)', 'keto_phenone_A(11)', 'cyano_pyridone_C(11)',
  'thiaz_ene_C(11)', 'hzone_thiophene_A(11)', 'ene_quin_methide(10)', 'het_thio_676_A(10)',
  'ene_five_het_G(10)', 'acyl_het_A(9)', 'anil_di_alk_G(9)', 'dhp_keto_A(9)', 'thio_urea_B(9)',
  'anil_alk_bim(9)', 'imine_imine_A(9)', 'thio_urea_C(9)', 'imine_one_fives_B(9)',
  'dhp_amino_CN_B(9)', 'anil_OC_no_alk_A(8)', 'het_thio_66_one(8)', 'styrene_B(8)',
  'het_thio_5_A(8)', 'anil_di_alk_ene_A(8)', 'ene_rhod_D(8)', 'ene_rhod_E(8)', 'anil_OH_alk_A(8)',
  'pyrrole_C(8)', 'thio_urea_D(8)', 'thiaz_ene_D(8)', 'ene_rhod_F(8)', 'thiaz_ene_E(8)',
  'het_65_B(7)', 'keto_keto_beta_C(7)', 'het_66_A(7)', 'thio_urea_E(7)', 'thiophene_amino_C(7)',
  'hzone_phenone(7)', 'ene_rhod_G(7)', 'ene_cyano_B(7)', 'dhp_amino_CN_C(7)', 'het_5_A(7)',
  'ene_five_het_H(6)', 'thio_amide_A(6)', 'ene_cyano_C(6)', 'hzone_furan_A(6)', 'anil_di_alk_H(6)',
  'het_65_C(6)', 'thio_urea_F(6)', 'ene_five_het_I(6)', 'keto_keto_gamma(5)', 'quinone_B(5)',
  'het_6_pyridone_OH(5)', 'hzone_naphth_A(5)', 'thio_ester_A(5)', 'ene_misc_A(5)',
  'cyano_pyridone_D(5)', 'het_65_Db(5)', 'het_666_A(5)', 'diazox_sulfon_B(5)', 'anil_NH_alk_A(5)',
  'sulfonamide_C(5)', 'het_thio_N_55(5)', 'keto_keto_beta_D(5)', 'ene_rhod_H(5)', 'imine_ene_A(5)',
  'het_thio_656a(5)', 'pyrrole_D(5)', 'pyrrole_E(5)', 'thio_urea_G(5)', 'anisol_A(5)',
  'pyrrole_F(5)', 'dhp_amino_CN_D(5)', 'thiazole_amine_A(4)', 'het_6_imidate_A(4)',
  'anil_OC_no_alk_B(4)', 'styrene_C(4)', 'azulene(4)', 'furan_acid_A(4)', 'cyano_pyridone_E(4)',
  'anil_alk_thio(4)', 'anil_di_alk_I(4)', 'het_thio_6_furan(4)', 'anil_di_alk_ene_B(4)',
  'imine_one_B(4)', 'anil_OC_alk_A(4)', 'ene_five_het_J(4)', 'pyrrole_G(4)', 'ene_five_het_K(4)',
  'cyano_ene_amine_B(4)', 'thio_ester_B(4)', 'ene_five_het_L(4)', 'hzone_thiophene_B(4)',
  'dhp_amino_CN_E(4)', 'het_5_B(4)', 'imine_imine_B(3)', 'thiazole_amine_B(3)',
  'imine_ene_one_A(3)', 'diazox_A(3)', 'ene_one_A(3)', 'anil_OC_no_alk_C(3)', 'thiazol_SC_A(3)',
  'het_666_B(3)', 'furan_A(3)', 'colchicine_A(3)', 'thiophene_C(3)', 'anil_OC_alk_B(3)',
  'het_thio_66_A(3)', 'rhod_sat_B(3)', 'ene_rhod_I(3)', 'keto_thiophene(3)', 'imine_imine_C(3)',
  'het_65_pyridone_A(3)', 'thiazole_amine_C(3)', 'het_thio_pyr_A(3)', 'melamine_A(3)',
  'anil_NH_alk_B(3)', 'rhod_sat_C(3)', 'thiophene_amino_D(3)', 'anil_OC_alk_C(3)',
  'het_thio_65_A(3)', 'het_thio_656b(3)', 'thiazole_amine_D(3)', 'thio_urea_H(3)',
  'cyano_pyridone_F(3)', 'rhod_sat_D(3)', 'ene_rhod_J(3)', 'imine_phenol_A(3)',
  'thio_carbonate_B(3)', 'het_thio_N_5A(3)', 'het_thio_N_65A(3)', 'anil_di_alk_J(3)',
  'pyrrole_H(3)', 'ene_cyano_D(3)', 'cyano_cyano_B(3)', 'ene_five_het_M(3)', 'cyano_ene_amine_C(3)',
  'thio_urea_I(3)', 'dhp_amino_CN_F(3)', 'anthranil_acid_B(3)', 'diazox_B(3)', 'thio_aldehyd_A(3)',
  'thio_amide_B(2)', 'imidazole_B(2)', 'thiazole_amine_E(2)', 'thiazole_amine_F(2)',
  'thio_ester_C(2)', 'ene_one_B(2)', 'quinone_C(2)', 'keto_naphthol_A(2)', 'thio_amide_C(2)',
  'phthalimide_misc(2)', 'sulfonamide_D(2)', 'anil_NH_alk_C(2)', 'het_65_E(2)', 'hzide_naphth(2)',
  'anisol_B(2)', 'thio_carbam_ene(2)', 'thio_amide_D(2)', 'het_65_Da(2)', 'thiophene_D(2)',
  'het_thio_6_ene(2)', 'cyano_keto_A(2)', 'anthranil_acid_C(2)', 'naphth_amino_C(2)',
  'naphth_amino_D(2)', 'thiazole_amine_G(2)', 'het_66_B(2)', 'coumarin_A(2)', 'anthranil_acid_D(2)',
  'het_66_C(2)', 'thiophene_amino_E(2)', 'het_6666_A(2)', 'sulfonamide_E(2)', 'anil_di_alk_K(2)',
  'het_5_C(2)', 'ene_six_het_B(2)', 'steroid_A(2)', 'het_565_A(2)', 'thio_imine_ium(2)',
  'anthranil_acid_E(2)', 'hzone_furan_B(2)', 'thiophene_E(2)', 'ene_misc_B(2)', 'het_thio_5_B(2)',
  'thiophene_amino_F(2)', 'anil_OC_alk_D(2)', 'tert_butyl_A(2)', 'thio_urea_J(2)',
  'het_thio_65_B(2)', 'coumarin_B(2)', 'thio_urea_K(2)', 'thiophene_amino_G(2)', 'anil_NH_alk_D(2)',
  'het_thio_5_C(2)', 'thio_keto_het(2)', 'het_thio_N_5B(2)', 'quinone_D(2)',
  'anil_di_alk_furan_B(2)', 'ene_six_het_C(2)', 'het_55_A(2)', 'het_thio_65_C(2)', 'hydroquin_A(2)',
  'anthranil_acid_F(2)', 'pyrrole_I(2)', 'thiophene_amino_H(2)', 'imine_one_fives_C(2)',
  'keto_phenone_zone_A(2)', 'dyes7A(2)', 'het_pyridiniums_B(2)', 'het_5_D(2)',
  'thiazole_amine_H(1)', 'thiazole_amine_I(1)', 'het_thio_N_5C(1)', 'sulfonamide_F(1)',
  'thiazole_amine_J(1)', 'het_65_F(1)', 'keto_keto_beta_E(1)', 'ene_five_one_B(1)',
  'keto_keto_beta_zone(1)', 'thio_urea_L(1)', 'het_thio_urea_ene(1)', 'cyano_amino_het_A(1)',
  'tetrazole_hzide(1)', 'imine_naphthol_A(1)', 'misc_anisole_A(1)', 'het_thio_665(1)',
  'anil_di_alk_L(1)', 'colchicine_B(1)', 'misc_aminoacid_A(1)', 'imidazole_amino_A(1)',
  'phenol_sulfite_A(1)', 'het_66_D(1)', 'misc_anisole_B(1)', 'tetrazole_A(1)', 'het_65_G(1)',
  'misc_trityl_A(1)', 'misc_pyridine_OC(1)', 'het_6_hydropyridone(1)', 'misc_stilbene(1)',
  'misc_imidazole(1)', 'anil_NH_no_alk_A(1)', 'het_6_imidate_B(1)', 'anil_alk_B(1)',
  'styrene_anil_A(1)', 'misc_aminal_acid(1)', 'anil_no_alk_D(1)', 'anil_alk_C(1)',
  'misc_anisole_C(1)', 'het_465_misc(1)', 'anthranil_acid_G(1)', 'anil_di_alk_M(1)',
  'anthranil_acid_H(1)', 'thio_urea_M(1)', 'thiazole_amine_K(1)', 'het_thio_5_imine_A(1)',
  'thio_amide_E(1)', 'het_thio_676_B(1)', 'sulfonamide_G(1)', 'thio_thiomorph_Z(1)',
  'naphth_ene_one_A(1)', 'naphth_ene_one_B(1)', 'amino_acridine_A(1)', 'keto_phenone_B(1)',
  'hzone_acid_A(1)', 'sulfonamide_H(1)', 'het_565_indole(1)', 'pyrrole_J(1)', 'pyrazole_amino_B(1)',
  'pyrrole_K(1)', 'anthranil_acid_I(1)', 'thio_amide_F(1)', 'ene_one_C(1)', 'het_65_H(1)',
  'cyano_imine_D(1)', 'cyano_misc_A(1)', 'ene_misc_C(1)', 'het_66_E(1)', 'keto_keto_beta_F(1)',
  'misc_naphthimidazole(1)', 'naphth_ene_one_C(1)', 'keto_phenone_C(1)', 'coumarin_C(1)',
  'thio_est_cyano_A(1)', 'het_65_imidazole(1)', 'anthranil_acid_J(1)', 'colchicine_het(1)',
  'ene_misc_D(1)', 'indole_3yl_alk_B(1)', 'anil_OH_no_alk_A(1)', 'thiazole_amine_L(1)',
  'pyrazole_amino_A(1)', 'het_thio_N_5D(1)', 'anil_alk_indane(1)', 'anil_di_alk_N(1)',
  'het_666_C(1)', 'ene_one_D(1)', 'anil_di_alk_indol(1)', 'anil_no_alk_indol_A(1)',
  'dhp_amino_CN_G(1)', 'anil_di_alk_dhp(1)', 'anthranil_amide_A(1)', 'hzone_anthran_Z(1)',
  'ene_one_amide_A(1)', 'het_76_A(1)', 'thio_urea_N(1)', 'anil_di_alk_coum(1)',
  'ene_one_amide_B(1)', 'het_thio_656c(1)', 'het_5_ene(1)', 'thio_imide_A(1)', 'dhp_amidine_A(1)',
  'thio_urea_O(1)', 'anil_di_alk_O(1)', 'thio_urea_P(1)', 'het_pyraz_misc(1)', 'diazox_C(1)',
  'diazox_D(1)', 'misc_cyclopropane(1)', 'imine_ene_one_B(1)', 'coumarin_D(1)', 'misc_furan_A(1)',
  'rhod_sat_E(1)', 'rhod_sat_imine_A(1)', 'rhod_sat_F(1)', 'het_thio_5_imine_B(1)',
  'het_thio_5_imine_C(1)', 'ene_five_het_N(1)', 'thio_carbam_A(1)', 'misc_anilide_A(1)',
  'misc_anilide_B(1)', 'mannich_B(1)', 'mannich_catechol_A(1)', 'anil_alk_D(1)', 'het_65_I(1)',
  'misc_urea_A(1)', 'imidazole_C(1)', 'styrene_imidazole_A(1)', 'thiazole_amine_M(1)',
  'misc_pyrrole_thiaz(1)', 'pyrrole_L(1)', 'het_thio_65_D(1)', 'ene_misc_E(1)', 'thio_cyano_A(1)',
  'cyano_amino_het_B(1)', 'cyano_pyridone_G(1)', 'het_65_J(1)', 'ene_one_yne_A(1)',
  'anil_OH_no_alk_B(1)', 'hzone_acyl_misc_A(1)', 'thiophene_F(1)', 'anil_OC_alk_E(1)',
  'anil_OC_alk_F(1)', 'het_65_K(1)', 'het_65_L(1)', 'coumarin_E(1)', 'coumarin_F(1)',
  'coumarin_G(1)', 'coumarin_H(1)', 'het_thio_67_A(1)', 'sulfonamide_I(1)', 'het_65_mannich(1)',
  'anil_alk_A(1)', 'het_5_inium(1)', 'anil_di_alk_P(1)', 'thio_urea_Q(1)', 'thio_pyridine_A(1)',
  'melamine_B(1)', 'misc_phthal_thio_N(1)', 'hzone_acyl_misc_B(1)', 'tert_butyl_B(1)',
  'diazox_E(1)', 'anil_NH_no_alk_B(1)', 'anil_no_alk_A(1)', 'anil_no_alk_B(1)',
  'thio_ene_amine_A(1)', 'het_55_B(1)', 'cyanamide_A(1)', 'ene_one_one_A(1)', 'ene_six_het_D(1)',
  'ene_cyano_E(1)', 'ene_cyano_F(1)', 'hzone_furan_C(1)', 'anil_no_alk_C(1)', 'hzone_acid_D(1)',
  'hzone_furan_E(1)', 'het_6_pyridone_NH2(1)', 'imine_one_fives_D(1)', 'pyrrole_M(1)',
  'pyrrole_N(1)', 'pyrrole_O(1)', 'ene_cyano_G(1)', 'sulfonamide_J(1)', 'misc_pyrrole_benz(1)',
  'thio_urea_R(1)', 'ene_one_one_B(1)', 'dhp_amino_CN_H(1)', 'het_66_anisole(1)',
  'thiazole_amine_N(1)', 'het_pyridiniums_C(1)', 'het_5_E(1)'
]

DIRNAME = os.path.join(os.path.dirname(sys.argv[0]))
PAINS_CSV = os.path.join(DIRNAME, "..", "..", "..", "Data", "Pains", "wehi_pains.csv")
PAINS_A_FILENAME = os.path.join(DIRNAME, "pains_a.in")
PAINS_B_FILENAME = os.path.join(DIRNAME, "pains_b.in")
PAINS_C_FILENAME = os.path.join(DIRNAME, "pains_c.in")

for fn in (PAINS_CSV, PAINS_A_FILENAME, PAINS_B_FILENAME, PAINS_B_FILENAME):
  if not os.path.exists(fn):
    raise IOError("Could not find necessary file: %s", fn)

sa = set(pains_a)
sb = set(pains_b)
sc = set(pains_c)

PAINS = {}

with open(PAINS_CSV) as fh:
  for smiles, name in csv.reader(fh):
    name = name.replace("<regId=", "").replace(">", "")
    PAINS[name] = smiles

PAINS_A = []
for n in pains_a:
  PAINS_A.append((n, PAINS[n]))
  del PAINS[n]
PAINS_B = []
for n in pains_b:
  PAINS_B.append((n, PAINS[n]))
  del PAINS[n]
PAINS_C = []
for n in pains_c:
  PAINS_C.append((n, PAINS[n]))
  del PAINS[n]
assert not PAINS

PAINS_A = ['{"%s","%s",0,""}' % (name, smiles) for name, smiles in PAINS_A]
PAINS_A = "const FilterData_t PAINS_A[] = {\n%s\n};" % ",\n".join(PAINS_A)

PAINS_B = ['{"%s","%s",0,""}' % (name, smiles) for name, smiles in PAINS_B]
PAINS_B = "const FilterData_t PAINS_B[] = {\n%s\n};" % ",\n".join(PAINS_B)

PAINS_C = ['{"%s","%s",0,""}' % (name, smiles) for name, smiles in PAINS_C]
PAINS_C = "const FilterData_t PAINS_C[] = {\n%s\n};" % ",\n".join(PAINS_C)


def write_pains(filename, data):
  with open(filename) as fh:
    t = fh.read()
  if t != data:
    if py3:
      import io

      # newline = don't convert to windows style
      with io.open(filename, 'w', newline='') as f:
        f.write(data)
    else:
      # wb = means don't convert newline to windows style
      with open(filename, 'wb') as f:
        f.write(data)


write_pains(PAINS_A_FILENAME, PAINS_A)
write_pains(PAINS_B_FILENAME, PAINS_B)
write_pains(PAINS_C_FILENAME, PAINS_C)

print("== Done updating pains files")
