import sys
from rdkit import Chem
import threading
import multiprocessing

# this just tests some threading stuff to ensure it doesn't crash with python
#  releasing the GIL smarts are recursive...
ref_sdf = '\n     RDKit          3D\n\n 22 23  0  0  1  0  0  0  0  0999 V2000\n   -6.1917   -1.9517    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.0664   -1.3009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.9401   -1.9499    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8148   -1.2991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6885   -1.9483    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5632   -1.2973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5642    0.0027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6905    0.6517    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6916    1.9517    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8158    0.0009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.9422    0.6501    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.0685    1.2991    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5632   -1.9465    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6885   -1.2955    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6874    0.0046    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8148   -1.9447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9401   -1.2936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9391    0.0064    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0644    0.6572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.1907    0.0082    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.1917   -1.2918    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0664   -1.9429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4 10  2  0\n  4  5  1  0\n  5  6  2  0\n  6  7  1  0\n  6 13  1  0\n  7  8  2  0\n  8  9  1  0\n  8 10  1  0\n 10 11  1  0\n 11 12  3  0\n 13 14  1  0\n 14 15  2  0\n 14 16  1  0\n 16 17  1  0\n 17 22  1  0\n 17 18  2  0\n 18 19  1  0\n 19 20  2  0\n 20 21  1  0\n 21 22  2  0\nM  END'
ref_mol = Chem.MolFromMolBlock(ref_sdf)

core_smarts = '[#6]-!@[#6]-!@[#8]-!@[#6]:1:[#6](-!@[#6]#!@[#7]):[#6](-!@[#7]):[#6]:[#6](-!@[#7]-!@[#6](-!@[#6]-!@[#6]:2:[#6]:[#6]:[#6]:[#6]:[#6]:2)=!@[#8]):[#7]:1'
if ref_mol is None:
    raise ValueError('Bad ref structure')
core_mol = Chem.MolFromSmarts(core_smarts)
if core_mol is None:
    raise ValueError('Bad core structure')

expected = {}


def runner(func, args):
    if args:
        res = getattr(ref_mol, func)(args)
    else:
        res = getattr(ref_mol, func)()
    if func in expected:
        assert res == expected[func], "Got %r expected %r" % (ers, expected[func])
    return res


funcs = ["GetSubstructMatch", "GetSubstructMatches", "HasSubstructMatch"]

# get the expected results from the non-thread version
for func in funcs:
    expected[func] = runner(func, core_mol)

nthreads = int(multiprocessing.cpu_count() * 100 / 4)  # 100 threads per cpu
threads = []
for i in range(0, nthreads):
    for func in funcs:
        t = threading.Thread(target=runner, args=(func, core_mol))
        t.start()
        threads.append(t)
    t = threading.Thread(target=runner, args=("ToBinary", None))
    t.start()
    threads.append(t)
for t in threads:
    t.join()
