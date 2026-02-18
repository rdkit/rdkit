"""This code allows the user to calculate CNS-MPO for a single (or multiple) molecule with SMILES and pKa values deposited in the lists.
   There is no need to create any file with SMILES and pKa values, so it is ideal for quick analysis.
   To use this code, user must have installed pandas and RDKit in virtual environment.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from math import log10


class CNS_MPO_single_molecule:
    def __init__(self, smiles_list, pKa_list):
        if len(smiles_list) != len(pKa_list):
            raise ValueError(
                "Length of smiles_list must be equal to length of pKa_list"
            )
        self.smiles_list = smiles_list
        self.pKa_list = pKa_list
        self._df = None
        self.calculate()

    def clogD(self, logP, pKa, pH=7.4):
        return logP - log10(1 + 10 ** (pH - pKa))

    def csv_file_preparation(self):
        dictionary = {"MW": [], "LogP": [], "HBD": [], "TPSA": []}

        for cpd in self.smiles_list:
            molecule = Chem.MolFromSmiles(cpd)
            if molecule is None:
                # Handle invalid SMILES
                dictionary["MW"].append(None)
                dictionary["LogP"].append(None)
                dictionary["HBD"].append(None)
                dictionary["TPSA"].append(None)
                continue
            mol_mw = Descriptors.MolWt(molecule)
            mol_logp = Crippen.MolLogP(molecule)
            mol_hbd = rdMolDescriptors.CalcNumHBD(molecule)
            mol_tpsa = Descriptors.TPSA(molecule)

            dictionary["MW"].append(mol_mw)
            dictionary["LogP"].append(mol_logp)
            dictionary["HBD"].append(mol_hbd)
            dictionary["TPSA"].append(mol_tpsa)

        df_descriptors = pd.DataFrame(dictionary)
        df_descriptors["pKa"] = self.pKa_list
        df_descriptors["LogD"] = df_descriptors.apply(
            lambda x: self.clogD(x["LogP"], x["pKa"]), axis=1
        )
        return df_descriptors

    def mw_score_func(self, mw):
        if mw is None:
            return 0
        if mw <= 360:
            return 1
        elif 360 < mw <= 500:
            return -0.005 * mw + 2.5
        else:
            return 0

    def logp_score_func(self, logp):
        if logp is None:
            return 0
        if logp <= 3:
            return 1
        elif 3 < logp <= 5:
            return -0.5 * logp + 2.5
        else:
            return 0

    def logd_score_func(self, logd):
        if logd is None:
            return 0
        if logd <= 2:
            return 1
        elif 2 < logd <= 4:
            return -0.5 * logd + 2
        else:
            return 0

    def pka_score_func(self, pka):
        if pka is None:
            return 0
        if pka <= 8:
            return 1
        elif 8 < pka <= 10:
            return -0.5 * pka + 5
        else:
            return 0

    def tpsa_score_func(self, tpsa):
        if tpsa is None:
            return 0
        if 40 <= tpsa <= 90:
            return 1
        elif 90 < tpsa <= 120:
            return -0.0333 * tpsa + 4
        elif 20 <= tpsa < 40:
            return 0.05 * tpsa - 1
        else:
            return 0

    def hbd_score_func(self, hbd):
        if hbd is None:
            return 0
        if hbd == 0:
            return 1
        elif hbd == 1:
            return 0.75
        elif hbd == 2:
            return 0.5
        elif hbd == 3:
            return 0.25
        else:
            return 0

    def calcCNS_MPO(self):
        df_descriptors = self.csv_file_preparation()
        df_descriptors["MW_score"] = df_descriptors["MW"].apply(
            self.mw_score_func
        )
        df_descriptors["LogP_score"] = df_descriptors["LogP"].apply(
            self.logp_score_func
        )
        df_descriptors["LogD_score"] = df_descriptors["LogD"].apply(
            self.logd_score_func
        )
        df_descriptors["pKa_score"] = df_descriptors["pKa"].apply(
            self.pka_score_func
        )
        df_descriptors["TPSA_score"] = df_descriptors["TPSA"].apply(
            self.tpsa_score_func
        )
        df_descriptors["HBD_score"] = df_descriptors["HBD"].apply(
            self.hbd_score_func
        )

        df_descriptors["CNS_MPO"] = (
            df_descriptors["MW_score"]
            + df_descriptors["LogP_score"]
            + df_descriptors["LogD_score"]
            + df_descriptors["pKa_score"]
            + df_descriptors["TPSA_score"]
            + df_descriptors["HBD_score"]
        )

        return df_descriptors[
            ["MW", "LogP", "HBD", "TPSA", "pKa", "LogD", "CNS_MPO"]
        ]

    def calculate(self):
        self._df = self.calcCNS_MPO()

    def __getitem__(self, key):
        if isinstance(key, (int, slice)):
            return self._df.iloc[key]
        elif isinstance(key, str):
            return self._df[key]
        else:
            raise KeyError(f"Unsupported key type: {type(key)}")

    def __repr__(self):
        return repr(self._df)

    def __str__(self):
        return str(self._df)



"""Example of use"""
# from cns_mpo_single_molecule import CNS_MPO_single_molecule
#
# x = CNS_MPO_single_molecule(
#     smiles_list=[
#         "C(C=1CCN(C5)CCC(C5)c(n4)c(c3o4)ccc(c3)F)(=O)N(C2)C(CCC2)=NC1C",
#         "C1CNCC(C1C2=CC=C(C=C2)F)COC3=CC4=C(C=C3)OCO4",
#     ],
#     pKa_list=[8.765, 9.365],
# )
# print(x)

## It is worth mentioning that, SMILES and pKa values have to be deposited in lists,
## regardless of whether we administer one or more molecules.

### Authors: Adam Mazur, Rafał Kurczab