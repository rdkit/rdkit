#!/usr/bin/env python
# Converts the pickled fragment-contribution table used by
# $RDBASE/Contrib/SA_Score/sascorer.py into the flat binary form read by
# RDKit::Descriptors::SAScoreFragmentTable.  Re-run after regenerating
# fpscores.pkl.gz with $RDBASE/Contrib/SA_Score/build_sascore_db.py:
#
#   python build_fpscores_bin.py ../../Contrib/SA_Score/fpscores.pkl.gz fpscores.bin
#
# Layout (little-endian):
#     magic     8s                b"RDSASCR1"
#     nScores   uint32
#     nKeys     uint32
#     scores    float32 * nScores   distinct contributions, ascending
#     keys      uint32  * nKeys     Morgan bit ids, ascending
#     scoreIdx  uint16  * nKeys     index into scores[]
#
# The scores are stored once and referenced by index because the table holds
# far more fragments than distinct contributions (705k vs 3.5k).

import argparse
import gzip
import pickle
import struct
from pathlib import Path

MAGIC = b"RDSASCR1"


def encode(rows):
  """ packs the pickle's [[score, bitId, ...], ...] rows into the binary table

  """
  contribs = {}
  for row in rows:
    for bit in row[1:]:
      contribs[bit] = float(row[0])

  scores = sorted(set(contribs.values()))
  if len(scores) > 0xFFFF:
    raise ValueError(f"{len(scores)} distinct scores will not fit a uint16 index")
  rank = {s: i for i, s in enumerate(scores)}

  keys = sorted(contribs)
  idx = [rank[contribs[k]] for k in keys]

  return b"".join([
    MAGIC,
    struct.pack("<II", len(scores), len(keys)),
    struct.pack(f"<{len(scores)}f", *scores),
    struct.pack(f"<{len(keys)}I", *keys),
    struct.pack(f"<{len(keys)}H", *idx),
  ])


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("pickle_file", type=Path, help="input fpscores.pkl.gz")
  parser.add_argument("output_file", type=Path, help="output fpscores.bin")
  args = parser.parse_args()

  rows = pickle.load(gzip.open(args.pickle_file))
  blob = encode(rows)
  args.output_file.write_bytes(blob)

  n_scores, n_keys = struct.unpack_from("<II", blob, 8)
  print(f"wrote {args.output_file}: {len(blob):,} bytes, "
        f"{n_keys:,} fragments, {n_scores:,} distinct scores")


if __name__ == "__main__":
  main()
