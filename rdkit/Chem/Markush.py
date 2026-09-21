#  Copyright (C) 2026 RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
"""Utilities for working with finite Markush spaces.

The generic-group query support in :mod:`rdkit.Chem` supplies the Markush
matching semantics.  This module provides a small, explicit representation of
a *finite disjunction* of such queries and operations on a supplied candidate
library.  A library is required for enumeration because generic groups (for
example ``ARY``) describe an unbounded chemical space.
"""

from rdkit import Chem


class MarkushFormula:
  """A disjunction of RDKit query molecules.

  Args:
    queries: Query molecules defining the alternatives in the formula.

  The queries are retained by reference.  Use :func:`MakeMarkushFormula` to
  build a formula which explicitly covers a collection of molecules.
  """

  def __init__(self, queries):
    self.queries = tuple(queries)
    if not self.queries:
      raise ValueError("a Markush formula needs at least one query")
    if any(query is None for query in self.queries):
      raise ValueError("Markush queries cannot be None")


class _ExactMarkushFormula(MarkushFormula):

  def __init__(self, queries, identities):
    super().__init__(queries)
    self._identities = frozenset(identities)


def _generic_match_parameters():
  params = Chem.SubstructMatchParameters()
  params.useGenericMatchers = True
  params.useChirality = True
  return params


def _as_formula(formula):
  if not isinstance(formula, MarkushFormula):
    formula = MarkushFormula((formula, ))
  return formula


def IsInMarkushScope(formula, molecule, params=None):
  """Return whether ``molecule`` matches at least one formula alternative.

  Generic Markush matchers are enabled by default. Pass an explicit
  :class:`~rdkit.Chem.rdchem.SubstructMatchParameters` instance to change the
  matching policy.
  """
  if molecule is None:
    raise ValueError("molecule cannot be None")
  formula = _as_formula(formula)
  if isinstance(formula, _ExactMarkushFormula):
    return Chem.MolToSmiles(molecule) in formula._identities
  if params is None:
    params = _generic_match_parameters()
  return any(molecule.HasSubstructMatch(query, params) for query in formula.queries)


def EnumerateMarkush(formula, candidates, params=None):
  """Return unique members of ``candidates`` that belong to ``formula``.

  Candidates are processed in input order. Chemical identity is determined by
  canonical isomeric SMILES, so duplicate graph representations occur once in
  the result.
  """
  result = []
  seen = set()
  formula = _as_formula(formula)
  for molecule in candidates:
    if molecule is None:
      raise ValueError("candidate molecules cannot be None")
    if IsInMarkushScope(formula, molecule, params):
      key = Chem.MolToSmiles(molecule)
      if key not in seen:
        seen.add(key)
        result.append(molecule)
  return tuple(result)


def MakeMarkushFormula(molecules):
  """Build an explicit Markush formula that covers ``molecules``.

  Each distinct input molecule becomes one SMARTS alternative. This is a
  conservative construction: it always covers all supplied molecules and
  avoids inventing unverified R-group generalizations.
  """
  queries = []
  identities = []
  seen = set()
  for molecule in molecules:
    if molecule is None:
      raise ValueError("molecules cannot contain None")
    key = Chem.MolToSmiles(molecule)
    if key not in seen:
      seen.add(key)
      identities.append(key)
      queries.append(Chem.MolFromSmarts(Chem.MolToSmarts(molecule)))
  return _ExactMarkushFormula(queries, identities)
