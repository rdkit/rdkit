from rdkit import Chem
from rdkit.Chem import AllChem
from collections import defaultdict
import math
from itertools import combinations

def compute_wl_signatures(mol, max_iter=5):
    """Weisfeiler-Lehman算法生成原子环境特征"""
    atoms = mol.GetAtoms()
    num_atoms = len(atoms)
    # 预计算每个原子的邻居索引
    neighbors = [sorted(n.GetIdx() for n in atom.GetNeighbors()) for atom in atoms]
    # 使用列表存储颜色，索引对应原子ID
    colors = [atom.GetSymbol() for atom in atoms]
    for _ in range(max_iter):
        color_map = {}
        new_colors = [None] * num_atoms
        for a_idx in range(num_atoms):
            neighbor_list = neighbors[a_idx]
            neighbor_colors = sorted(colors[n] for n in neighbor_list)
            new_color = (colors[a_idx], tuple(neighbor_colors))
            if new_color not in color_map:
                color_map[new_color] = len(color_map)
            new_colors[a_idx] = color_map[new_color]
        if new_colors == colors:
            break
        colors = new_colors.copy()
    # 转换为字典
    return dict(enumerate(colors))

def smiles_to_n(smiles):
    """转换SMILES并计算分子复杂度"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    adj_dict = defaultdict(float)
    for bond in mol.GetBonds():
        start, end = sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
        bond_type = bond.GetBondType()
        cnt = 1.5 if bond_type == Chem.BondType.AROMATIC else \
              1.0 if bond_type == Chem.BondType.SINGLE else \
              2.0 if bond_type == Chem.BondType.DOUBLE else \
              3.0 if bond_type == Chem.BondType.TRIPLE else 0.0
        adj_dict[(start, end)] += cnt

    # 计算eta
    atom_neighbors = defaultdict(dict)
    for (a, b), cnt in adj_dict.items():
        atom_neighbors[a][b] = cnt
        atom_neighbors[b][a] = cnt

    eta_case1 = 0.0
    for b in atom_neighbors:
        counts = list(atom_neighbors[b].values())
        sum_counts = sum(counts)
        sum_sq = sum(c * c for c in counts)
        eta_case1 += (sum_counts ** 2 - sum_sq) / 2

    eta_case2 = sum(k * (k - 1) / 2 for k in adj_dict.values() if k > 1)
    eta = eta_case1 + eta_case2

    # 计算原子环境特征
    wl_colors = compute_wl_signatures(mol)

    # 优化eta_i计算
    eta_i = defaultdict(float)
    # case 1：A-B-C chain，使用combinations生成对
    for b_atom in atom_neighbors:
        b_color = wl_colors[b_atom]
        neighbors = list(atom_neighbors[b_atom].items())
        for (a1, cnt1), (a2, cnt2) in combinations(neighbors, 2):
            key = ('chain', b_color, tuple(sorted((wl_colors[a1], wl_colors[a2]))))
            eta_i[key] += cnt1 * cnt2

    # case 2：A-B loop
    for (a, b), cnt in adj_dict.items():
        if cnt > 1:
            key = ('loop', tuple(sorted((wl_colors[a], wl_colors[b]))))
            eta_i[key] += cnt * (cnt - 1) / 2

    # 计算修正项 C(eta)
    entropy_term = sum(v * math.log2(v) if v > 0 else 0 for v in eta_i.values())
    c_eta = 2 * eta * math.log2(eta) - entropy_term if eta != 0 else 0

    # 计算E和E_i
    atom_types = defaultdict(int)
    for atom in mol.GetAtoms():
        atom_types[atom.GetSymbol()] += 1
    E = sum(atom_types.values())
    entropy_term_E = sum(v * math.log2(v) for v in atom_types.values() if v > 0)
    c_E = E * math.log2(E) - entropy_term_E if E != 0 else 0

    C_T = c_eta + c_E
    return eta, dict(eta_i), c_eta, E, atom_types, c_E, C_T

# 使用示例
if __name__ == "__main__":
    smiles = "O=C1[C@@]23C([C@]([C@]4([H])C1)([H])CC[C@@]([C@]4(C)C5)([H])CC6=C5N=C(C[C@@](CC[C@]7([H])[C@]8([H])C[C@@H](O)[C@@]9(C)C7=C[C@H]%10[C@]9(O)[C@H](C)[C@]%11([C@H](O)C[C@](CO)(C)O%11)O%10)([H])[C@]8(C)C%12)C%12=N6)=CC[C@]2([H])[C@H](C)[C@]%13(OC(C)(C)C[C@H]%13O)OC3"
    eta, eta_i, c_eta, E, E_i, c_E, C_T = smiles_to_n(smiles)
    print(f"Total eta: {eta}")
    print("Eta_i categories:")
    for k, v in eta_i.items():
        print(f"{k}: {v}")
    print(f"Correction term C(eta): {c_eta}")
    print(f"Total number of atoms (E): {E}")
    print("E_i categories:")
    for k, v in E_i.items():
        print(f"{k}: {v}")
    print(f"Correction term C(E): {c_E}")
    print(f"Total Complexity (C_T): {C_T}")