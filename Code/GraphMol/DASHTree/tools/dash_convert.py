#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Convert a DASH tree into RDKit's .dash container.

The DASH-tree python package stores a tree as a pair of files per atom-feature
branch:

  <b>.gz   gzipped pickle of a list of node tuples
           (id, atom_type, con_atom, con_type, attention, [child ids])
  <b>.h5   pandas HDF5 table (Fixed format), one row per node

Most of that is redundant, which this script verifies rather than assumes: the
h5 columns ``atom_type``, ``con_atom`` and ``con_type`` are identical to node
fields 1/2/3, ``level`` is the node's depth, ``max_attention`` is node field 4
rounded to float16, and the node tuple's leading id is always its own index.
Only the child lists and the value columns carry information.

The container it writes is described in ../DESIGN.md. Two decisions shape it:

  * Nodes are renumbered in breadth-first order across the whole forest. That
    makes each node's children a contiguous run of ids, so no child-index array
    is needed at all, and it puts the shallow levels every descent passes
    through into one small region at the front of the file -- which is what
    keeps the working set to a couple of pages per atom instead of ten.
  * A node's match key, child range and stop flag live in one 8-byte record, so
    a descent step is a single dependent load.

Usage:
    python dash_convert.py <tree_folder> <out.dash> [--props result,std,...]
"""

import argparse
import gzip
import os
import pickle
import struct
import sys

import numpy as np

MAGIC = b"DASHTREE"
FORMAT_VERSION = 2
ENDIAN_ID = 0xDEADBEEF

HEADER_SIZE = 128
PROP_ENTRY_SIZE = 48
PROP_NAME_LEN = 32
ALIGN = 64
NODE_RECORD_SIZE = 8

DTYPE_F16, DTYPE_F32, DTYPE_I32, DTYPE_F64 = 1, 2, 3, 4
_NP_DTYPE = {DTYPE_F16: np.float16, DTYPE_F32: np.float32,
             DTYPE_I32: np.int32, DTYPE_F64: np.float64}

# key packing, mirroring AtomFeatures.h
KEY_ATOM_BITS, KEY_CON_ATOM_BITS, KEY_CON_TYPE_BITS = 8, 4, 3
MAX_ATOM_TYPE = (1 << KEY_ATOM_BITS) - 1     # 255
MAX_CON_ATOM = (1 << KEY_CON_ATOM_BITS) - 2  # 14, stored as conAtom + 1
MAX_CON_TYPE = (1 << KEY_CON_TYPE_BITS) - 2  # 6,  stored as conType + 1

MAX_CHILDREN = 255          # numChildren is one byte
MAX_NODES = 0xFFFFFFFF      # firstChild is four bytes

NODE_FLAG_STOP = 0x01

# the default cumulative-attention threshold the stop flag is precomputed for
DEFAULT_ATTENTION_THRESHOLD = 10.0

# h5 columns that merely duplicate the topology and are therefore never stored
REDUNDANT_COLUMNS = ("level", "atom_type", "con_atom", "con_type", "max_attention")

# atomic number of the feature class whose atoms are matched through their heavy
# neighbour; the descent spends its first level getting there and does not
# accumulate that level's attention, which the stop flag has to account for
HYDROGEN = 1


def pack_key(atom_type, con_atom, con_type):
    return (atom_type | ((con_atom + 1) << KEY_ATOM_BITS)
            | ((con_type + 1) << (KEY_ATOM_BITS + KEY_CON_ATOM_BITS)))


def align_up(n):
    return (n + ALIGN - 1) // ALIGN * ALIGN


def load_branch(folder, b):
    with gzip.open(os.path.join(folder, f"{b}.gz"), "rb") as f:
        return pickle.load(f)


def hydrogen_branch(n_branches):
    """Branch index of the hydrogen feature class, or -1."""
    try:
        sys.path.insert(0, os.environ.get("DASH_TREE_REPO", ""))
        from serenityff.charge.tree.atom_features import AtomFeatures
        for i, f in enumerate(AtomFeatures.feature_list[:n_branches]):
            if f[0] == HYDROGEN:
                return i
    except ImportError:
        pass
    return -1


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tree_folder")
    ap.add_argument("out")
    ap.add_argument("--props", default=None,
                    help="comma separated property columns to keep "
                         "(default: every non-redundant column in the h5 files)")
    ap.add_argument("--float32", action="store_true",
                    help="narrow every float column to float32. Halves the "
                         "float64 columns at the cost of reproducing them exactly.")
    ap.add_argument("--no-source-ids", action="store_true",
                    help="omit the map back to the python package's node numbering. "
                         "It is never read unless a caller asks for it, and the "
                         "bit-exactness tests do.")
    ap.add_argument("--hydrogen-branch", type=int, default=None,
                    help="branch index of the hydrogen feature class "
                         "(default: read it from the serenityff package)")
    args = ap.parse_args(argv)

    import pandas as pd

    folder = args.tree_folder
    n_branches = len([f for f in os.listdir(folder) if f.endswith(".gz")])
    print(f"branches: {n_branches}")

    h_branch = (args.hydrogen_branch if args.hydrogen_branch is not None
                else hydrogen_branch(n_branches))
    if h_branch < 0:
        sys.exit("could not determine the hydrogen branch index; pass "
                 "--hydrogen-branch (it decides where the attention "
                 "accumulation starts, so it must be right)")
    print(f"hydrogen branch: {h_branch}")

    # ---- pass 1: sizes and the property column list ------------------------
    branch_n = np.zeros(n_branches, dtype=np.int64)
    for b in range(n_branches):
        branch_n[b] = len(load_branch(folder, b))
    n_nodes = int(branch_n.sum())
    if n_nodes > MAX_NODES:
        sys.exit(f"{n_nodes} nodes exceeds the container's 32-bit node ids")
    branch_root = np.zeros(n_branches, dtype=np.uint32)
    acc = 0
    for b in range(n_branches):
        branch_root[b] = acc
        acc += int(branch_n[b])
    print(f"nodes: {n_nodes:,}")

    df0 = pd.read_hdf(os.path.join(folder, "0.h5"), key="df", mode="r")
    if args.props:
        props = [p.strip() for p in args.props.split(",") if p.strip()]
        missing = [p for p in props if p not in df0.columns]
        if missing:
            sys.exit(f"requested columns not in the data: {missing}\n"
                     f"available: {list(df0.columns)}")
        redundant = [p for p in props if p in REDUNDANT_COLUMNS]
        if redundant:
            sys.exit(f"these columns duplicate the topology and are not stored: "
                     f"{redundant}")
    else:
        props = [c for c in df0.columns if c not in REDUNDANT_COLUMNS]
    print(f"properties: {props}")
    print(f"  (skipped as redundant: "
          f"{[c for c in df0.columns if c in REDUNDANT_COLUMNS]})")

    def out_dtype(col):
        """Keep each column at the width its source actually needs.

        The MBIS charges are float16 and round-trip through it exactly; the
        DASH-properties columns (AM1BCC, DFTD4, ...) are genuine float64 and do
        not survive narrowing, so storing them as float32 would quietly lose
        bits. --float32 narrows them anyway, for when size matters more than
        reproducing the source exactly.
        """
        d = df0[col].dtype
        if np.issubdtype(d, np.integer):
            return DTYPE_I32
        if args.float32:
            return DTYPE_F32
        if d == np.float16:
            return DTYPE_F16
        if d == np.float32:
            return DTYPE_F32
        return DTYPE_F64

    prop_dtypes = [out_dtype(c) for c in props]
    for c, dt in zip(props, prop_dtypes):
        print(f"    {c:<28} {df0[c].dtype} -> {_NP_DTYPE[dt].__name__}")

    # ---- pass 2: renumber breadth-first, one branch at a time --------------
    # Within a branch the roots come first and then every node's children in
    # tree order, which is exactly breadth-first: processing nodes in id order
    # and appending their children makes each node's children a contiguous run,
    # so no child-index array is needed, and it gathers the shallow levels every
    # descent passes through into a small region at the start of the branch.
    node_key = np.zeros(n_nodes, dtype=np.uint16)
    node_children = np.zeros(n_nodes, dtype=np.uint8)
    node_first = np.zeros(n_nodes, dtype=np.uint32)
    node_flags = np.zeros(n_nodes, dtype=np.uint8)
    node_attn = np.zeros(n_nodes, dtype=np.float32)
    node_source = np.zeros(n_nodes, dtype=np.uint32)
    prop_arrays = [np.empty(n_nodes, dtype=_NP_DTYPE[dt]) for dt in prop_dtypes]

    max_atom_type = max_con_atom = max_con_type = -(10 ** 9)
    min_con_atom = min_con_type = 10 ** 9
    max_children = 0
    n_stop = 0

    for b in range(n_branches):
        tree = load_branch(folder, b)
        base = int(branch_root[b])
        count = len(tree)
        cum = np.zeros(count, dtype=np.float64)

        root = tree[0]
        node_key[base] = pack_key(int(root[1]), int(root[2]), int(root[3]))
        node_source[base] = 0
        node_attn[base] = np.float32(root[4])
        cum[0] = 0.0  # a root's own attention is never accumulated

        next_local = 1
        for local in range(count):
            src = int(node_source[base + local])
            # The root of every branch lists itself as its own first child. That
            # edge can never be matched: its key carries conAtom = -1, while
            # every candidate the descent generates carries a real position, and
            # the one candidate that does use -1 -- a hydrogen's heavy neighbour
            # -- cannot carry the hydrogen feature class. So it is dropped
            # instead of encoded.
            children = [c for c in tree[src][5] if c != src]
            if src == 0 and len(children) != len(tree[0][5]) - 1:
                sys.exit(f"branch {b}: root does not list itself exactly once")
            if src != 0 and len(children) != len(tree[src][5]):
                sys.exit(f"branch {b}: node {src} references itself")
            if len(children) > MAX_CHILDREN:
                sys.exit(f"node {src} of branch {b} has {len(children)} children, "
                         f"more than the container's {MAX_CHILDREN}")
            max_children = max(max_children, len(children))
            node_children[base + local] = len(children)
            node_first[base + local] = base + next_local

            # the first step out of a hydrogen root only walks to the heavy
            # neighbour, and the descent does not accumulate its attention
            accumulate = not (src == 0 and b == h_branch)
            for c in children:
                cn = tree[c]
                at, ca, ct = int(cn[1]), int(cn[2]), int(cn[3])
                max_atom_type = max(max_atom_type, at)
                max_con_atom = max(max_con_atom, ca)
                min_con_atom = min(min_con_atom, ca)
                max_con_type = max(max_con_type, ct)
                min_con_type = min(min_con_type, ct)
                if not 0 <= at <= MAX_ATOM_TYPE:
                    sys.exit(f"atom_type {at} does not fit {KEY_ATOM_BITS} bits")
                if not -1 <= ca <= MAX_CON_ATOM:
                    sys.exit(f"con_atom {ca} does not fit {KEY_CON_ATOM_BITS} bits")
                if not -1 <= ct <= MAX_CON_TYPE:
                    sys.exit(f"con_type {ct} does not fit {KEY_CON_TYPE_BITS} bits")

                attention = np.float32(cn[4])
                node_key[base + next_local] = pack_key(at, ca, ct)
                node_source[base + next_local] = c
                node_attn[base + next_local] = attention
                cum[next_local] = cum[local] + (float(attention) if accumulate
                                                else 0.0)
                next_local += 1

        if next_local != count:
            sys.exit(f"branch {b}: renumbering covered {next_local} of {count} "
                     f"nodes -- the tree is not a single rooted tree")

        stop = cum > DEFAULT_ATTENTION_THRESHOLD
        node_flags[base:base + count] = np.where(stop, NODE_FLAG_STOP,
                                                 0).astype(np.uint8)
        n_stop += int(stop.sum())
        del tree, cum

        # the property columns, indexed by the new numbering
        d = pd.read_hdf(os.path.join(folder, f"{b}.h5"), key="df", mode="r")
        if len(d) != count:
            sys.exit(f"branch {b}: {count} nodes but {len(d)} data rows")
        rows = node_source[base:base + count]
        for arr, col, dt in zip(prop_arrays, props, prop_dtypes):
            arr[base:base + count] = d[col].to_numpy()[rows].astype(
                _NP_DTYPE[dt], copy=False)
        del d

    print(f"  nodes past the default attention threshold: {n_stop:,}")
    print(f"key ranges: atom_type <= {max_atom_type}   "
          f"con_atom [{min_con_atom},{max_con_atom}]   "
          f"con_type [{min_con_type},{max_con_type}]   "
          f"max children {max_children}")

    # ---- write -------------------------------------------------------------
    records = np.zeros(n_nodes, dtype=[("key", "<u2"), ("nChildren", "u1"),
                                      ("flags", "u1"), ("firstChild", "<u4")])
    records["key"] = node_key[:n_nodes]
    records["nChildren"] = node_children[:n_nodes]
    records["flags"] = node_flags[:n_nodes]
    records["firstChild"] = node_first[:n_nodes]
    assert records.itemsize == NODE_RECORD_SIZE, records.itemsize

    blocks = [("branchRoot", branch_root),
              ("nodeRecord", records),
              ("nodeAttn", node_attn[:n_nodes])]
    if not args.no_source_ids:
        blocks.append(("nodeSourceId", node_source[:n_nodes]))
    blocks += [(f"prop:{c}", a) for c, a in zip(props, prop_arrays)]

    prop_dir_off = HEADER_SIZE
    off = align_up(HEADER_SIZE + PROP_ENTRY_SIZE * len(props))
    offsets = {}
    for name, arr in blocks:
        offsets[name] = off
        off = align_up(off + arr.nbytes)
    file_size = off
    have_source = 0 if args.no_source_ids else offsets["nodeSourceId"]

    with open(args.out, "wb") as f:
        f.write(struct.pack(
            "<8sIIIIIIdQQQQQQ",
            MAGIC, FORMAT_VERSION, ENDIAN_ID, n_branches, len(props),
            n_nodes, max_children, DEFAULT_ATTENTION_THRESHOLD,
            offsets["branchRoot"], offsets["nodeRecord"], offsets["nodeAttn"],
            have_source, prop_dir_off, file_size))
        f.write(struct.pack("<hhh", max_atom_type, max_con_atom, max_con_type))
        f.write(b"\0" * (HEADER_SIZE - f.tell()))

        for col, dt in zip(props, prop_dtypes):
            nm = col.encode("utf-8")
            if len(nm) >= PROP_NAME_LEN:
                sys.exit(f"property name too long: {col}")
            f.write(nm + b"\0" * (PROP_NAME_LEN - len(nm)))
            f.write(struct.pack("<B7xQ", dt, offsets[f"prop:{col}"]))

        for name, arr in blocks:
            f.write(b"\0" * (offsets[name] - f.tell()))
            f.write(arr.tobytes())
        f.write(b"\0" * (file_size - f.tell()))

    print(f"\nwrote {args.out}  {os.path.getsize(args.out)/1e6:.1f} MB")
    for name, arr in blocks:
        print(f"    {name:<28} n={arr.size:>10,}  {arr.nbytes/1e6:>7.1f} MB "
              f"@ {offsets[name]}")


if __name__ == "__main__":
    main()
