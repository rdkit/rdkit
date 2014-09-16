# $Id$
#
# Copyright (C) 2001-2010 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
from rdkit.Chem import rdchem
import shlex



def _from_json_table(f):
    columns = []
    while True:
        line = f.readline().strip()
        if ":::" in line: break
        if "#" in line: continue
        columns.append(line)
        
    types = []
    for k in columns: 
        if k[0] == "i": types.append(int)
        elif k[0] == "r": types.append(float)
        elif k[0] == "b": types.append(bool)
        elif k[0] == "s": types.append(str)        
        
    rows = []
    while True:
        line = f.readline().strip()
        if ":::" in line: break
        if "}" in line: break
        rows.append(shlex.split(line))

    results = []

    for r in rows:
        row = {}
        for c, callback, v in zip(columns, types, r[1:]):
            # "<>" appears to mark a null value of some sort
            if v == "<>": row[c] = None
            else: row[c] = callback(v)
        results.append(row)
    return results
  

def _from_json_dict(f):
    keys = []
    while True:
        line = f.readline()
        if ":::" in line: break
        if "#" in line: continue
        keys.append(line.strip())

    types = []
    for k in keys: 
        if k[0] == "i": types.append(int)
        elif k[0] == "r": types.append(float)
        elif k[0] == "b": types.append(bool)
        elif k[0] == "s": types.append(str)
        
    values = []
    
    m_blocks = {}
    
    while True:
        
        line = f.readline().strip()
        if line.startswith("m_"): 
            name, rest = line.split("[")
            elements = int(rest.split("]")[0])
            d = _from_json_table(f)
            assert(len(d) == elements)
            m_blocks[name] = d
            f.readline().strip()
        if "}" in line: break
        if not line: continue
        # if sm_title is empty it is stored as ""
        values.append(line.replace('"','')) 
    
    assert(len(types) == len(keys))
    
    values = [callback(v) for callback, v in zip(types, values)]
    
    ret = dict(zip(keys, values))
    ret.update(m_blocks)
    return ret


def _from_json(f):
    line = f.readline()
    blocks = []
    while line:
        line = line.strip()
        # some_name {
        if "{" in line and not ("[" in line or "]" in line):             
            blocks.append(_from_json_dict(f))
        # some_name[27] {
        elif "{" in line and ("[" in line and "]" in line):
            blocks.append(_from_json_table(f))
        line = f.readline()
    return blocks
  
class ForwardMAEParser(object):
    def __init__(self, f):
        if isinstance(f, str):
            self.file = open(f) 
        elif isinstance(f, file):
            self.file = f
        elif isinstance(f, gzip.GzipFile):
            self.file = f
        else: raise TypeError("f has to be str, file or gzip.GzipFile")    
        
    def __iter__(self):
        return self
    
    def next(self):
        f = self.file
        line = f.readline()
        
        while line:
            line = line.strip()
            # some_name {
            if "{" in line and not ("[" in line or "]" in line):             
                return _from_json_dict(f)
            # some_name[27] {
            elif "{" in line and ("[" in line and "]" in line):
                return _from_json_table(f)
            line = f.readline()
        raise StopIteration


from rdkit.Geometry import Point3D
class ForwardMAESupplier(ForwardMAEParser):
    def next(self):
        # read maestro block as a dictionary
        block = super(ForwardMAESupplier, self).next()
        
        # create editable molecule
        mol = rdchem.Mol()
        mol = rdchem.EditableMol(mol)
        
        if not "m_atom" in block: return None
        
        # Add atoms
        for row in block["m_atom"]: 
            a = rdchem.Atom(row["i_m_atomic_number"])
            [a.SetProp(k, str(v)) for k, v in row.items() if not "coord" in k]
            mol.AddAtom(a)
        
        # Add bonds
        bond_orders = {1: rdchem.BondType.SINGLE, 2: rdchem.BondType.DOUBLE, 3: rdchem.BondType.TRIPLE}
        for row in block["m_bond"]:
            mol.AddBond(row["i_m_from"]-1, row["i_m_to"]-1, bond_orders[row["i_m_order"]] )
        
        # Add conformer coordinates
        mol = mol.GetMol()
        c = rdchem.Conformer(len(block["m_atom"]))
        for i, row in enumerate(block["m_atom"]): c.SetAtomPosition(i, Point3D(row["r_m_x_coord"], row["r_m_y_coord"], row["r_m_z_coord"]))
        
        # Set molecule properties
        [mol.SetProp(k, str(v)) for k, v in block.items() if not k.startswith("m_")] 
        
        
        return mol
        