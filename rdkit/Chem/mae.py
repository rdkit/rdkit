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
from rdkit.Chem.rdmolops import SanitizeMol
from rdkit.Geometry.rdGeometry import Point3D
import gzip
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
        elif k[0] == "b": types.append(int)
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
        values.append(line) 
    
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
        
        self.counter = 0
        
    def __iter__(self):
        return self
    
    def next(self):
        f = self.file
        line = f.readline()
        
        while line:
            line = line.strip()
            # some_name {
            if "{" in line and not ("[" in line or "]" in line):  
                self.counter += 1 
                if not self.counter % 1000: print self.counter           
                return _from_json_dict(f)
            # some_name[27] {
            elif "{" in line and ("[" in line and "]" in line):
                self.counter += 1
                if not self.counter % 1000: print self.counter
                return _from_json_table(f)
            line = f.readline()
        self.file.close()
        raise StopIteration



class ForwardMAEMolSupplier(ForwardMAEParser):
    def next(self):
        # read maestro block as a dictionary
        block = super(ForwardMAEMolSupplier, self).next()
        
        # create editable molecule
        mol = rdchem.Mol()
        mol = rdchem.EditableMol(mol)
        
        if not "m_atom" in block: return self.next()
        
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
        mol.AddConformer(c)
        
        # Set molecule properties
        [mol.SetProp(k, str(v)) for k, v in block.items() if not k.startswith("m_")] 
        
        ret = SanitizeMol(mol, catchErrors=True)
        
        if ret: return None
        
        return mol
        
def MolToDict(mol, confId=-1):
    types = {"i" : int, "r": float, "b": int, "s": str}
    properties = [(e, types[e[0]], m.GetProp(e)) for e in m.GetPropNames()]    
    properties = dict([(name, callback(value)) for name, callback, value in properties if callback])
    
    bond_types = {rdchem.BondType.SINGLE:1, rdchem.BondType.DOUBLE:2, rdchem.BondType.TRIPLE:3}
    bonds = [{'i_m_from': b.GetBeginAtomIdx()+1, 
               'i_m_to': b.GetEndAtomIdx()+1, 
               'i_m_order':bond_types[ b.GetBondType()]} for b in m.GetBonds()]
    
    atoms = [[ (p, types[p[0]], a.GetProp(p)) for p in a.GetPropNames() ] for a in m.GetAtoms()]
    atoms = [dict([(name, callback(value))for name, callback, value in a]) for a in atoms]
    
    c = mol.GetConformer(confId)
    positions = [c.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    for a, p in zip(atoms, positions): 
        a["r_m_x_coord"] = p.x
        a["r_m_y_coord"] = p.y
        a["r_m_z_coord"] = p.z
        
    properties["m_atom"] = atoms
    properties["m_bond"] = bonds    
    
    return properties

def MolToMAEBlock(mol, confId=-1, convert=True):
    if convert:
        mol = MolToDict(mol)
    keys = sorted([k for k in mol.keys() if not k.startswith("m_")])
    lines = ["f_m_ct { "] + ["  {}".format(k)for k in keys] + ["  :::"] + ["  {}".format(str(mol[k])) for k in keys]
    
    
    if "m_depend" in mol:
        depend_columns = ["i_m_depend_dependency","s_m_depend_property"]
        lines = lines + ["  m_depend[{}] {{".format(len(mol["m_depend"]))]
        lines = lines + ["    {}".format(c) for c in depend_columns] + ["    :::"]
        lines = lines + ["    {} {} {}".format(*tuple([i+1]+[a[c] for c in depend_columns])) for i, a in enumerate(mol["m_depend"])]
        lines = lines + ["    :::\n  }"]
    
    atom_columns = ["i_m_mmod_type","r_m_x_coord","r_m_y_coord","r_m_z_coord","i_m_residue_number","i_m_color","r_m_charge1","r_m_charge2","i_m_atomic_number","i_i_constraint"]
    lines = lines + ["  m_atom[{}] {{".format(len(mol["m_atom"]))]
    lines = lines + ["    {}".format(c) for c in atom_columns] + ["    :::"]
    lines = lines + ["    {} {} {} {} {} {} {} {} {} {} {}".format(*tuple([i+1]+[a[c] for c in atom_columns])) for i, a in enumerate(mol["m_atom"])]
    lines = lines + ["    :::\n  }"]
    
    bond_columns = ["i_m_from", "i_m_to", "i_m_order"]
    lines = lines + ["  m_bond[{}] {{".format(len(mol["m_bond"]))]
    lines = lines + ["    {}".format(c) for c in bond_columns] + ["    :::"]
    lines = lines + ["    {} {} {} {}".format(*tuple([i+1]+[a[c] for c in bond_columns])) for i, a in enumerate(mol["m_bond"])]
    lines = lines + ["    :::\n  }"]
    
    lines = lines + ["}\n\n\n\n"]
    return "\n".join(lines)
  
class DictWriter():
    def __init__(self, filename):
        f = open(filename, "w")
        header = """{ 
  s_m_m2io_version
  :::
  2.0.0 
}\n\n"""
        f.write(header)
        self.f = f 
        
    def write(self, mol, confId=-1):
        self.f.write(MolToMAEBlock(mol, confId=confId, convert=False))
  
class MAEWriter(DictWriter):
    def write(self, mol, confId=-1):
        self.f.write(MolToMAEBlock(mol=mol, confId=confId, convert=True))  
  