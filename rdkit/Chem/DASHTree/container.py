#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""The .dash container in python: its constants, a reader and the writer.

The tools write containers through this module, so the layout lives here
once. The C++ side of it is Code/GraphMol/DASHTree/DASHTreeImpl.h. All data is
little-endian, and every array is 64-byte aligned.
"""

import struct

import numpy as np
from rdkit.Chem import rdDASHTree

MAGIC = b"DASHTREE"
FORMAT_VERSION = 3
ENDIAN_ID = 0xDEADBEEF

HEADER_SIZE = 128
# the fields of DASHHeader, up to the feature table. The rest is zero padding.
HEADER_FORMAT = "<8sIIIIIIdQQQQQQhhhHQ"
PROP_ENTRY_FORMAT = "<32sB7xQ"
PROP_ENTRY_SIZE = 48
PROP_NAME_LEN = 32
ALIGN = 64

RECORD_DTYPE = np.dtype([("key", "<u2"), ("nChildren", "u1"), ("flags", "u1"),
                         ("firstChild", "<u4")])
FEATURE_DTYPE = np.dtype([("atomicNum", "u1"), ("degree", "u1"), ("formalCharge", "i1"),
                          ("conjugated", "u1"), ("numHs", "u1"), ("reserved", "3u1")])

DTYPE_F16, DTYPE_F32, DTYPE_I32, DTYPE_F64 = 1, 2, 3, 4
# what an int32 column stores for "no value here". Float columns use NaN.
MISSING_INT32 = np.iinfo(np.int32).min
NP_DTYPE = {
  DTYPE_F16: np.dtype("<f2"),
  DTYPE_F32: np.dtype("<f4"),
  DTYPE_I32: np.dtype("<i4"),
  DTYPE_F64: np.dtype("<f8"),
}

# key packing, as in MatchKey.h
KEY_ATOM_BITS, KEY_CON_ATOM_BITS, KEY_CON_TYPE_BITS = 8, 4, 3
MAX_ATOM_TYPE = (1 << KEY_ATOM_BITS) - 1  # 255
MAX_CON_ATOM = (1 << KEY_CON_ATOM_BITS) - 2  # 14, stored as conAtom + 1
MAX_CON_TYPE = (1 << KEY_CON_TYPE_BITS) - 2  # 6,  stored as conType + 1

MAX_CHILDREN = 255  # numChildren is one byte
MAX_NODES = 0xFFFFFFFF  # firstChild is four bytes

NODE_FLAG_STOP = 0x01

# the cumulative-attention threshold of the precomputed stop flag
DEFAULT_ATTENTION_THRESHOLD = 10.0


def packKey(atomType, conAtom, conType):
  conTypeShift = KEY_ATOM_BITS + KEY_CON_ATOM_BITS
  return atomType | ((conAtom + 1) << KEY_ATOM_BITS) | ((conType + 1) << conTypeShift)


def unpackKey(key):
  """(atomType, conAtom, conType) of a packed key. It also works on arrays."""
  key = np.asarray(key, dtype=np.int64)
  atom = key & MAX_ATOM_TYPE
  conAtom = ((key >> KEY_ATOM_BITS) & ((1 << KEY_CON_ATOM_BITS) - 1)) - 1
  conType = ((key >> (KEY_ATOM_BITS + KEY_CON_ATOM_BITS)) & ((1 << KEY_CON_TYPE_BITS) - 1)) - 1
  return atom, conAtom, conType


def alignUp(n):
  return (n + ALIGN - 1) // ALIGN * ALIGN


class Container:
  """A mapped .dash file that exposes its arrays as read-only numpy views.

  Attributes: nBranches, nNodes, threshold, branchRoot (uint32), features
  (FEATURE_DTYPE, one per branch), records (RECORD_DTYPE), attn (float32) and
  props, an ordered dict name -> (dtype tag, array).
  """

  def __init__(self, path):
    rdDASHTree.DASHTree(str(path))  # the C++ reader checks the layout
    self.path = path
    self._mm = np.memmap(path, dtype=np.uint8, mode="r")
    (_, _, _, self.nBranches, nProps, self.nNodes, _, self.threshold, offRoots, offRec, offAttn, _,
     offDir, _, _, _, _, nFeatures, offFeatures) = struct.unpack_from(HEADER_FORMAT, self._mm, 0)
    self.branchRoot = self._view(offRoots, np.dtype("<u4"), self.nBranches)
    self.features = self._view(offFeatures, FEATURE_DTYPE, nFeatures)
    self.records = self._view(offRec, RECORD_DTYPE, self.nNodes)
    self.attn = self._view(offAttn, np.dtype("<f4"), self.nNodes)
    self.props = {}
    for i in range(nProps):
      name, tag, off = struct.unpack_from(PROP_ENTRY_FORMAT, self._mm, offDir + i * PROP_ENTRY_SIZE)
      name = name.rstrip(b"\0").decode("utf-8")
      self.props[name] = (tag, self._view(off, NP_DTYPE[tag], self.nNodes))

  def _view(self, offset, dtype, count):
    return np.frombuffer(self._mm, dtype=dtype, count=count, offset=offset)

  def parentAndLevel(self):
    """parent[n] (-1 at a root) and level[n] (0 at a root) for every node."""
    rec = self.records
    parent = np.full(self.nNodes, -1, dtype=np.int64)
    level = np.full(self.nNodes, -1, dtype=np.int16)
    frontier = self.branchRoot.astype(np.int64)
    level[frontier] = 0
    depth = 0
    reached = frontier.size
    while frontier.size:
      depth += 1
      has = rec["nChildren"][frontier] > 0
      parents = frontier[has]
      counts = rec["nChildren"][parents].astype(np.int64)
      starts = rec["firstChild"][parents].astype(np.int64)
      # with the check after the loop, every node is reached exactly once
      reached += counts.sum()
      if reached > self.nNodes or (starts + counts > self.nNodes).any():
        raise ValueError(f"{self.path}: the child links do not form a tree")
      children = (np.repeat(starts, counts) + np.arange(counts.sum()) -
                  np.repeat(np.cumsum(counts) - counts, counts))
      parent[children] = np.repeat(parents, counts)
      level[children] = depth
      frontier = children
    if (level < 0).any():
      raise ValueError(f"{self.path}: {int((level < 0).sum())} nodes are reachable from no root")
    return parent, level


def writeContainer(path, branchRoot, features, records, attn, props, threshold):
  """Writes a container.

  features has one FEATURE_DTYPE row per branch. props is a list of (name,
  dtype tag, array).
  """
  branchRoot = np.ascontiguousarray(branchRoot, dtype="<u4")
  features = np.ascontiguousarray(features, dtype=FEATURE_DTYPE)
  if features.size != branchRoot.size:
    raise ValueError(f"{features.size} atom-feature classes for {branchRoot.size} branches")
  records = np.ascontiguousarray(records, dtype=RECORD_DTYPE)
  attn = np.ascontiguousarray(attn, dtype="<f4")
  nNodes = records.size
  if branchRoot.size == 0 or nNodes == 0:
    raise ValueError("a container needs at least one branch and one node")
  if attn.size != nNodes:
    raise ValueError("every per-node array must have one entry per node")

  blocks = [("branchRoot", branchRoot), ("featureTable", features), ("nodeRecord", records),
            ("nodeAttn", attn)]
  names = [name for name, _, _ in props]
  if len(set(names)) != len(names):
    raise ValueError(f"property names repeat: {names}")
  for name, tag, arr in props:
    if len(name.encode("utf-8")) >= PROP_NAME_LEN:
      raise ValueError(f"property name is too long: {name}")
    arr = np.ascontiguousarray(arr, dtype=NP_DTYPE[tag])
    if arr.size != nNodes:
      raise ValueError(f"property {name} has {arr.size} values for {nNodes} nodes")
    blocks.append((f"prop:{name}", arr))

  propDirOff = HEADER_SIZE
  off = alignUp(HEADER_SIZE + PROP_ENTRY_SIZE * len(props))
  offsets = {}
  for name, arr in blocks:
    offsets[name] = off
    off = alignUp(off + arr.nbytes)
  fileSize = off

  with open(path, "wb") as f:
    f.write(
      struct.pack(HEADER_FORMAT, MAGIC, FORMAT_VERSION, ENDIAN_ID, branchRoot.size, len(props),
                  nNodes, int(records["nChildren"].max()), float(threshold), offsets["branchRoot"],
                  offsets["nodeRecord"], offsets["nodeAttn"], 0, propDirOff, fileSize,
                  *keyRangesOf(records, branchRoot), features.size, offsets["featureTable"]))
    f.write(b"\0" * (HEADER_SIZE - f.tell()))
    for name, tag, _ in props:
      f.write(struct.pack(PROP_ENTRY_FORMAT, name.encode("utf-8"), tag, offsets[f"prop:{name}"]))
    for name, arr in blocks:
      f.write(b"\0" * (offsets[name] - f.tell()))
      f.write(arr.tobytes())
    f.write(b"\0" * (fileSize - f.tell()))


def keyRangesOf(records, branchRoot):
  """(maxAtomType, maxConAtom, maxConType) over every node that is not a root,
  as the header records them."""
  isChild = np.ones(records.size, dtype=bool)
  isChild[np.asarray(branchRoot, dtype=np.int64)] = False
  if not isChild.any():
    return -1, -1, -1
  atom, conAtom, conType = unpackKey(records["key"][isChild])
  return int(atom.max()), int(conAtom.max()), int(conType.max())
