"""This code enable user to calculate CNS-MPO based on the CSV file, that contains SMILES and pKa values of the molecules.
   The advantage of this code is that the output can be easily used as input to the DataFrame, giving you access to free record manipulation.
   To use this code, user must have installed pandas and RDKit in virtual environment.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from math import log10


class CNS_MPO_csv_to_df:
    def __init__(self, csv_file):
        self.csv_file = csv_file
        self._df = None
        self.calculate()


    def read_csv(self):
        """This method loads CSV file with columns: 'Id', 'Smiles' and 'pKa', respectively.
           Method requires separating the contents of the CSV file with a semicolon (';') rather than a comma (',')"""
        df = pd.read_csv(self.csv_file, sep=";")
        df["pKa"] = df["pKa"].str.replace(",", ".").astype(float)
        return df

    def clogD(self, logP, pKa, pH=7.4):
        return logP - log10(1 + 10 ** (pH - pKa))

    def csv_file_preparation(self):
        dictionary = {"MW": [], "LogP": [], "HBD": [], "TPSA": []}

        for cpd in self._df["Smiles"]:
            molecule = Chem.MolFromSmiles(cpd)
            mol_mw = Descriptors.MolWt(molecule)
            mol_logp = Crippen.MolLogP(molecule)
            mol_hbd = rdMolDescriptors.CalcNumHBD(molecule)
            mol_tpsa = Descriptors.TPSA(molecule)

            dictionary["MW"].append(mol_mw)
            dictionary["LogP"].append(mol_logp)
            dictionary["HBD"].append(mol_hbd)
            dictionary["TPSA"].append(mol_tpsa)

        df_descriptors = pd.DataFrame(dictionary)
        df_descriptors["pKa"] = self._df["pKa"]
        df_descriptors["LogD"] = df_descriptors.apply(
            lambda x: self.clogD(x["LogP"], x["pKa"]), axis=1
        )
        return df_descriptors

    def mw_score_func(self, mw):
        if mw <= 360:
            return 1
        elif 360 < mw <= 500:
            return -0.005 * mw + 2.5
        else:
            return 0

    def logp_score_func(self, logp):
        if logp <= 3:
            return 1
        elif 3 < logp <= 5:
            return -0.5 * logp + 2.5
        else:
            return 0

    def logd_score_func(self, logd):
        if logd <= 2:
            return 1
        elif 2 < logd <= 4:
            return -0.5 * logd + 2
        else:
            return 0

    def pka_score_func(self, pka):
        if pka <= 8:
            return 1
        elif 8 < pka <= 10:
            return -0.5 * pka + 5
        else:
            return 0

    def tpsa_score_func(self, tpsa):
        if 40 <= tpsa <= 90:
            return 1
        elif 90 < tpsa <= 120:
            return -0.0333 * tpsa + 4
        elif 20 <= tpsa < 40:
            return 0.05 * tpsa - 1
        else:
            return 0

    def hbd_score_func(self, hbd):
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
        df_descriptors["Id"] = self._df["Id"]

        return df_descriptors[
            ["Id", "MW", "LogP", "LogD", "pKa", "TPSA", "HBD", "CNS_MPO"]
        ]

    def calculate(self):
        self._df = self.read_csv()
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

    def __iter__(self):
        return iter(self._df.to_dict(orient="records"))

    def to_dict(self):
        return self._df.to_dict(orient="records")

    def to_dataframe(self):
        return self._df



"""Example of use"""
# import pandas as pd
# from cns_mpo_csv_to_df import CNS_MPO_csv_to_df
#
# x = CNS_MPO_csv_to_df(csv_file_path)
# df = pd.DataFrame(x)
# print(df)

### Authors: Adam Mazur, Rafał Kurczab