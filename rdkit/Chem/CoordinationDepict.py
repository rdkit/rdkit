"""2D coordinate generation helpers for coordination complexes.

The standard RDKit coordinate generators are optimized for organic molecular
graphs.  In a coordination complex that often places independently drawable
ligands on top of each other around the metal atom.  This module keeps ligand
depiction delegated to CoordGen and only handles the missing coordination-
sphere layout step.

The initial implementation is intentionally limited to mononuclear complexes.
Molecules without exactly one bonded metal centre fall back to CoordGen.
"""

import math

import numpy as np

from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Geometry import Point3D


_METAL_ATOMIC_NUMBERS = (
    set(range(21, 31))
    | set(range(39, 49))
    | set(range(57, 81))
    | set(range(89, 104))
    | {13, 31, 49, 50, 51, 81, 82, 83, 84}
)


def _metal_indices(mol):
  return [
    atom.GetIdx() for atom in mol.GetAtoms()
    if atom.GetAtomicNum() in _METAL_ATOMIC_NUMBERS
  ]


def _rotation(angle):
  c, s = math.cos(angle), math.sin(angle)
  return np.array(((c, -s), (s, c)))


def _coordinates(mol):
  conf = mol.GetConformer()
  return np.array([
    (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
    for i in range(mol.GetNumAtoms())
  ])


def _median_bond_length(mol, coords):
  lengths = [
    np.linalg.norm(coords[bond.GetBeginAtomIdx()] -
                   coords[bond.GetEndAtomIdx()])
    for bond in mol.GetBonds()
  ]
  return float(np.median(lengths)) if lengths else 1.0


def _donor_groups(mol, donors):
  """Return connected components of donor atoms.

  Three or more contiguous donors represent one haptic interaction.  Keeping
  them in one group prevents an eta-bound ring from consuming five independent
  slots around the coordination sphere.
  """
  donors = set(donors)
  groups = []
  while donors:
    seed = donors.pop()
    group = {seed}
    pending = [seed]
    while pending:
      atom_idx = pending.pop()
      for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx in donors:
          donors.remove(neighbor_idx)
          group.add(neighbor_idx)
          pending.append(neighbor_idx)
    if len(group) >= 3:
      groups.append(sorted(group))
    else:
      groups.extend([idx] for idx in sorted(group))
  return sorted(groups, key=lambda group: group[0])


def _anchors(coords, groups):
  return np.array([
    coords[group].mean(axis=0) if len(group) >= 3 else coords[group[0]]
    for group in groups
  ])


def _component_without_bond(mol, start, blocked):
  component = {start}
  pending = [start]
  while pending:
    atom_idx = pending.pop()
    for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
      edge = frozenset((atom_idx, neighbor.GetIdx()))
      if edge == blocked or neighbor.GetIdx() in component:
        continue
      component.add(neighbor.GetIdx())
      pending.append(neighbor.GetIdx())
  return component


def _reflect(coords, begin, end, atom_indices):
  axis = coords[end] - coords[begin]
  length = np.linalg.norm(axis)
  if length < 1.0e-8:
    return coords
  axis /= length
  reflected = coords.copy()
  relative = reflected[atom_indices] - reflected[begin]
  reflected[atom_indices] = (
    reflected[begin] + 2.0 * np.outer(relative @ axis, axis) - relative)
  return reflected


def _nonbonded_clashes(mol, coords):
  distances = np.linalg.norm(coords[:, None] - coords[None, :], axis=-1)
  mask = np.triu(np.ones(distances.shape, dtype=bool), 1)
  for bond in mol.GetBonds():
    mask[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()] = False
    mask[bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()] = False
  return int(np.count_nonzero(distances[mask] < 0.55))


def _compact_chelate(mol, coords, donors):
  """Choose a compact 2D rotamer for a chelating ligand."""
  if len(donors) < 2:
    return coords

  candidates = []
  for first, second in zip(donors[:-1], donors[1:]):
    path = Chem.GetShortestPath(mol, int(first), int(second))
    for begin, end in zip(path[:-1], path[1:]):
      bond = mol.GetBondBetweenAtoms(begin, end)
      if bond.IsInRing() or bond.GetBondType() != Chem.BondType.SINGLE:
        continue
      blocked = frozenset((begin, end))
      side = _component_without_bond(mol, end, blocked)
      if begin in side:
        continue
      candidate = (begin, end, sorted(side))
      if candidate not in candidates:
        candidates.append(candidate)
  if not candidates:
    return coords

  # Enumerating the few rotatable bonds between donors is deterministic and
  # avoids a local search selecting an outward-facing chelate.
  candidates = candidates[:9]
  best = coords
  best_score = math.inf
  for mask in range(1 << len(candidates)):
    trial = coords
    for bit, (begin, end, side) in enumerate(candidates):
      if mask & (1 << bit):
        trial = _reflect(trial, begin, end, side)
    donor_coords = trial[donors]
    separation = max(
      np.linalg.norm(donor_coords[i] - donor_coords[j])
      for i in range(len(donors)) for j in range(i))
    score = abs(float(separation) - 1.8) + 10.0 * _nonbonded_clashes(
      mol, trial)
    if score < best_score:
      best = trial
      best_score = score
  return best


def _fit_ligand(coords, anchors, metal_bond_length):
  """Place a ligand in a canonical sector centred on the positive x axis."""
  if len(anchors) == 1:
    donor = anchors[0]
    ligand_centre = coords.mean(axis=0)
    outward = ligand_centre - donor
    if np.linalg.norm(outward) < 1.0e-8:
      outward = np.array((1.0, 0.0))
    angle = -math.atan2(outward[1], outward[0])
    return (coords - donor) @ _rotation(angle).T + np.array(
      (metal_bond_length, 0.0))

  source_centre = anchors.mean(axis=0)
  centred = anchors - source_centre
  if len(anchors) == 2:
    donor_distance = float(np.linalg.norm(anchors[1] - anchors[0]))
    half_angle = math.asin(min(donor_distance / (2.0 * metal_bond_length),
                               0.95))
    target_angles = np.array((-half_angle, half_angle))
  else:
    span = min(math.pi, math.radians(42.0 * (len(anchors) - 1)))
    target_angles = np.linspace(-span / 2.0, span / 2.0, len(anchors))
  target = metal_bond_length * np.column_stack(
    (np.cos(target_angles), np.sin(target_angles)))

  target_centre = target.mean(axis=0)
  candidates = []
  for target_points in (target, target[::-1]):
    covariance = centred.T @ (target_points - target_centre)
    left, _, right_t = np.linalg.svd(covariance)
    rotation = right_t.T @ left.T
    if np.linalg.det(rotation) < 0:
      right_t[-1] *= -1
      rotation = right_t.T @ left.T
    candidates.append(
      (coords - source_centre) @ rotation.T + target_centre)

  def placement_score(candidate):
    radii = np.linalg.norm(candidate, axis=1)
    return (np.count_nonzero(radii < 0.8 * metal_bond_length),
            -float(candidate[:, 0].mean()))

  return min(candidates, key=placement_score)


def _angular_width(coords):
  angles = np.arctan2(coords[:, 1], coords[:, 0])
  return max(float(np.max(np.abs(angles))), math.radians(8.0))


def _sector_centres(widths):
  padding = math.radians(12.0)
  weights = np.array([max(2.0 * width + padding, padding) for width in widths])
  weights *= 2.0 * math.pi / weights.sum()
  centres = []
  cursor = -math.pi
  for weight in weights:
    centres.append(cursor + weight / 2.0)
    cursor += weight
  return centres


def _separate_ligands(placements, centres, minimum_distance=0.65):
  if len(placements) < 2:
    return placements
  directions = [np.array((math.cos(angle), math.sin(angle)))
                for angle in centres]
  for offset in np.linspace(0.0, 3.0, 31):
    shifted = [coords + offset * direction
               for coords, direction in zip(placements, directions)]
    clearance = min(
      float(np.min(np.linalg.norm(first[:, None] - second[None, :], axis=-1)))
      for i, first in enumerate(shifted)
      for second in shifted[i + 1:])
    if clearance >= minimum_distance:
      return shifted
  return shifted


def _coordgen(mol):
  mol.RemoveAllConformers()
  rdCoordGen.AddCoords(mol)
  return mol.GetConformer().GetId()


def Compute2DCoordinationCoords(mol, metalBondLength=1.5):
  """Generate readable 2D coordinates for a mononuclear coordination complex.

  The molecule is modified in place and the new conformer id is returned,
  matching the coordinate-generator APIs elsewhere in RDKit.  Metal-donor
  bonds must be present in the molecular graph.  Molecules with zero or more
  than one bonded metal centre use regular CoordGen.

  Args:
    mol: molecule to modify.
    metalBondLength: initial distance between the metal and donor atoms.  The
      layout may increase all metal-ligand distances to prevent clashes.
  """
  if mol is None:
    raise ValueError("mol must not be None")
  if metalBondLength <= 0:
    raise ValueError("metalBondLength must be positive")

  metals = [
    idx for idx in _metal_indices(mol)
    if mol.GetAtomWithIdx(idx).GetDegree()
  ]
  if len(metals) != 1:
    return _coordgen(mol)

  metal_idx = metals[0]
  metal_bonds = [
    bond for bond in mol.GetAtomWithIdx(metal_idx).GetBonds()
    if bond.GetOtherAtomIdx(metal_idx) not in metals
  ]
  if not metal_bonds:
    return _coordgen(mol)

  donors = [bond.GetOtherAtomIdx(metal_idx) for bond in metal_bonds]
  fragmented = Chem.FragmentOnBonds(
    mol, [bond.GetIdx() for bond in metal_bonds], addDummies=False)
  mappings = []
  pieces = Chem.GetMolFrags(
    fragmented, asMols=True, sanitizeFrags=False,
    fragsMolAtomMapping=mappings)

  ligands = []
  spectators = []
  for piece, mapping in zip(pieces, mappings):
    mapping = list(mapping)
    if mapping == [metal_idx]:
      continue
    local_donors = [mapping.index(idx) for idx in donors if idx in mapping]
    if local_donors:
      ligands.append((piece, mapping, local_donors))
    else:
      spectators.append((piece, mapping))
  if not ligands:
    return _coordgen(mol)

  arranged = []
  for piece, mapping, local_donors in ligands:
    piece.RemoveAllConformers()
    rdCoordGen.AddCoords(piece)
    coords = _coordinates(piece)
    bond_length = _median_bond_length(piece, coords)
    if bond_length > 1.0e-8:
      coords /= bond_length
    groups = _donor_groups(piece, local_donors)
    coords = _compact_chelate(piece, coords, [group[0] for group in groups])
    fitted = _fit_ligand(coords, _anchors(coords, groups), metalBondLength)
    arranged.append((mapping, fitted, _angular_width(fitted)))

  centres = _sector_centres([entry[2] for entry in arranged])
  placements = [coords @ _rotation(angle).T
                for (_, coords, _), angle in zip(arranged, centres)]
  placements = _separate_ligands(placements, centres)
  output = np.zeros((mol.GetNumAtoms(), 2), dtype=float)
  assigned = np.zeros(mol.GetNumAtoms(), dtype=bool)
  output[metal_idx] = (0.0, 0.0)
  assigned[metal_idx] = True
  for (mapping, _, _), rotated in zip(arranged, placements):
    for local_idx, original_idx in enumerate(mapping):
      output[original_idx] = rotated[local_idx]
      assigned[original_idx] = True

  right_edge = float(output[assigned, 0].max()) if assigned.any() else 0.0
  for piece, mapping in spectators:
    piece.RemoveAllConformers()
    rdCoordGen.AddCoords(piece)
    coords = _coordinates(piece)
    coords -= coords.mean(axis=0)
    coords[:, 0] += right_edge + 3.0
    for local_idx, original_idx in enumerate(mapping):
      output[original_idx] = coords[local_idx]
      assigned[original_idx] = True
    right_edge = float(coords[:, 0].max())

  conformer = Chem.Conformer(mol.GetNumAtoms())
  conformer.Set3D(False)
  for idx, (x, y) in enumerate(output):
    conformer.SetAtomPosition(idx, Point3D(float(x), float(y), 0.0))
  mol.RemoveAllConformers()
  return mol.AddConformer(conformer, assignId=True)


def PrepareCoordinationMolForDrawing(mol, addCoords=True,
                                     metalBondLength=1.5):
  """Return a drawable copy with contiguous eta bonds collapsed to centroids.

  RDKit's existing ``DativeBondsToHaptic`` conversion is used, so endpoint
  metadata and dummy-atom conventions remain identical to other RDKit APIs.
  Plain single bonds are not reinterpreted as haptic interactions.
  """
  if mol is None:
    raise ValueError("mol must not be None")
  result = Chem.Mol(mol)
  if addCoords or not result.GetNumConformers():
    Compute2DCoordinationCoords(result, metalBondLength=metalBondLength)
  return Chem.DativeBondsToHaptic(result)


__all__ = [
  "Compute2DCoordinationCoords",
  "PrepareCoordinationMolForDrawing",
]
