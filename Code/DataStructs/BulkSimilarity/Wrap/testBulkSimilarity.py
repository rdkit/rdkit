#
#  Copyright (C) 2026 Clay Moore and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import random
import unittest

import numpy as np

from rdkit import DataStructs


def _random_fp(num_bits, density, rng):
  fp = DataStructs.ExplicitBitVect(num_bits)
  for i in range(num_bits):
    if rng.random() < density:
      fp.SetBit(i)
  return fp


def _random_fps(count, num_bits, rng, density=0.1):
  return [_random_fp(num_bits, density, rng) for _ in range(count)]


def _reference_matrix(probes, targets):
  out = np.empty((len(probes), len(targets)), dtype=np.float64)
  for i, p in enumerate(probes):
    for j, t in enumerate(targets):
      out[i, j] = DataStructs.TanimotoSimilarity(p, t)
  return out


class BulkSimilarityTests(unittest.TestCase):

  def setUp(self):
    self.rng = random.Random(0xC0FFEE)

  def test_returns_2d_float64_array_of_expected_shape(self):
    probes = _random_fps(3, 256, self.rng)
    targets = _random_fps(5, 256, self.rng)
    result = DataStructs.BulkTanimotoMatrix(probes, targets)
    self.assertIsInstance(result, np.ndarray)
    self.assertEqual(result.dtype, np.float64)
    self.assertEqual(result.shape, (3, 5))

  def test_matches_per_pair_TanimotoSimilarity_reference(self):
    probes = _random_fps(7, 512, self.rng)
    targets = _random_fps(11, 512, self.rng)
    actual = DataStructs.BulkTanimotoMatrix(probes, targets)
    expected = _reference_matrix(probes, targets)
    np.testing.assert_allclose(actual, expected, atol=1e-12)

  def test_handles_non_eight_word_fp_sizes(self):
    # 320 bits = 5 uint64_t words; exercises the AVX-512 tail mask path
    # when that kernel is active.
    probes = _random_fps(4, 320, self.rng)
    targets = _random_fps(6, 320, self.rng)
    actual = DataStructs.BulkTanimotoMatrix(probes, targets)
    expected = _reference_matrix(probes, targets)
    np.testing.assert_allclose(actual, expected, atol=1e-12)

  def test_identical_inputs_score_one(self):
    fp = _random_fp(128, 0.2, self.rng)
    result = DataStructs.BulkTanimotoMatrix([fp], [fp])
    self.assertEqual(result.shape, (1, 1))
    self.assertAlmostEqual(result[0, 0], 1.0)

  def test_zero_vectors_score_zero(self):
    a = DataStructs.ExplicitBitVect(128)
    b = DataStructs.ExplicitBitVect(128)
    result = DataStructs.BulkTanimotoMatrix([a], [b])
    self.assertEqual(result[0, 0], 0.0)

  def test_empty_inputs_return_empty_matrix(self):
    fp = _random_fp(128, 0.1, self.rng)
    self.assertEqual(DataStructs.BulkTanimotoMatrix([], [fp]).shape, (0, 1))
    self.assertEqual(DataStructs.BulkTanimotoMatrix([fp], []).shape, (1, 0))
    self.assertEqual(DataStructs.BulkTanimotoMatrix([], []).shape, (0, 0))

  def test_mismatched_bit_sizes_raise(self):
    a = DataStructs.ExplicitBitVect(128)
    b = DataStructs.ExplicitBitVect(256)
    with self.assertRaises(Exception):
      DataStructs.BulkTanimotoMatrix([a], [b])

  def test_non_explicit_bv_inputs_raise(self):
    a = DataStructs.ExplicitBitVect(128)
    with self.assertRaises(Exception):
      DataStructs.BulkTanimotoMatrix([a], ["not a fingerprint"])

  def test_active_kernel_is_known(self):
    self.assertIn(DataStructs.BulkSimilarityActiveKernel(),
                  ("scalar", "avx512vpopcntdq"))


if __name__ == "__main__":
  unittest.main()
