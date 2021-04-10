import ctypes
rdk = ctypes.cdll.LoadLibrary('./lib/librdkitcffi.so')
rdk.get_smiles.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_char_p]
rdk.get_smiles.restype = ctypes.c_char_p
rdk.get_mol.restype = ctypes.c_void_p
rdk.free_ptr.argtypes = [ctypes.c_void_p]

sz = ctypes.c_size_t()
mol = rdk.get_mol(b'CCC(O)C(C(=O)O)C', ctypes.byref(sz), b'')
smi = rdk.get_smiles(mol, sz, b'')
assert (smi == b'CCC(O)C(C)C(=O)O')
print(smi)
rdk.free_ptr(mol)
