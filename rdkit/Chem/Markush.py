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
    self.queries = _as_tuple(queries, "Markush queries")
    if not self.queries:
      raise ValueError("a Markush formula needs at least one query")
    if any(query is None for query in self.queries):
      raise ValueError("Markush queries cannot be None")
    if any(not isinstance(query, Chem.Mol) for query in self.queries):
      raise TypeError("Markush queries must be RDKit molecules")


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


def _as_tuple(values, name):
  try:
    return tuple(values)
  except TypeError:
    raise TypeError(f"{name} must be iterable") from None


def _require_iterable(values, name):
  try:
    iter(values)
  except TypeError:
    raise TypeError(f"{name} must be iterable") from None
  return values


def _require_molecule(molecule, none_message, type_message):
  if molecule is None:
    raise ValueError(none_message)
  if not isinstance(molecule, Chem.Mol):
    raise TypeError(type_message)


def _match_parameters(params):
  if params is None:
    return _generic_match_parameters()
  if not isinstance(params, Chem.SubstructMatchParameters):
    raise TypeError("params must be SubstructMatchParameters or None")
  return params


def IsInMarkushScope(formula, molecule, params=None):
  """Return whether ``molecule`` matches at least one formula alternative.

  Generic Markush matchers are enabled by default. Pass an explicit
  :class:`~rdkit.Chem.rdchem.SubstructMatchParameters` instance to change the
  matching policy.
  """
  _require_molecule(molecule, "molecule cannot be None",
                    "molecule must be an RDKit molecule")
  formula = _as_formula(formula)
  params = _match_parameters(params)
  if isinstance(formula, _ExactMarkushFormula):
    return Chem.MolToSmiles(molecule) in formula._identities
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
  params = _match_parameters(params)
  for molecule in _require_iterable(candidates, "candidates"):
    _require_molecule(molecule, "candidate molecules cannot be None",
                      "candidate molecules must be RDKit molecules")
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
  for molecule in _require_iterable(molecules, "molecules"):
    _require_molecule(molecule, "molecules cannot contain None",
                      "molecules must be RDKit molecules")
    key = Chem.MolToSmiles(molecule)
    if key not in seen:
      seen.add(key)
      identities.append(key)
      queries.append(Chem.MolFromSmarts(Chem.MolToSmarts(molecule)))
  return _ExactMarkushFormula(queries, identities)
