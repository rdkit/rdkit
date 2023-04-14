# $Id$
#
#  Copyright (C) 2007-2009 Greg Landrum
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os

from sqlalchemy import (Column, Float, Integer, LargeBinary, String, Text,
                        create_engine)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski

decBase = declarative_base()


class Compound(decBase):
  __tablename__ = 'molecules'
  guid = Column(Integer, primary_key=True)
  molpkl = Column(LargeBinary)


def RegisterSchema(dbUrl, echo=False):
  engine = create_engine(dbUrl, echo=echo)
  decBase.metadata.create_all(engine)
  return sessionmaker(bind=engine)


ConnectToSchema = RegisterSchema


def _ConnectToSchema(dbUrl, echo=False):
  engine = create_engine(dbUrl, echo=echo)
  decBase.metadata.create_all(engine)
  return sessionmaker(bind=engine)


#set up the logger:

import rdkit.RDLogger as logging

logger = logging.logger()
logger.setLevel(logging.INFO)


def ProcessMol(session, mol, globalProps, nDone, nameProp='_Name', nameCol='compound_id',
               redraw=False, keepHs=False, skipProps=False, addComputedProps=False,
               skipSmiles=False):
  if not mol:
    raise ValueError('no molecule')
  if keepHs:
    Chem.SanitizeMol(mol)
  try:
    nm = mol.GetProp(nameProp)
  except KeyError:
    nm = None
  if not nm:
    nm = f'Mol_{nDone}'

  cmpd = Compound()
  session.add(cmpd)

  if redraw:
    AllChem.Compute2DCoords(m)

  if not skipSmiles:
    cmpd.smiles = Chem.MolToSmiles(mol, True)
  cmpd.molpkl = mol.ToBinary()
  setattr(cmpd, nameCol, nm)

  if not skipProps:
    if addComputedProps:
      cmpd.DonorCount = Lipinski.NumHDonors(mol)
      cmpd.AcceptorCount = Lipinski.NumHAcceptors(mol)
      cmpd.RotatableBondCount = Lipinski.NumRotatableBonds(mol)
      cmpd.AMW = Descriptors.MolWt(mol)
      cmpd.MolLogP = Crippen.MolLogP(mol)
    pns = list(mol.GetPropNames())
    for pn in pns:
      if pn.lower() == nameCol.lower():
        continue
      pv = mol.GetProp(pn).strip()
      if pn in globalProps:
        setattr(cmpd, pn.lower(), pv)
  return cmpd


def LoadDb(suppl, dbName, nameProp='_Name', nameCol='compound_id', silent=False, redraw=False,
           errorsTo=None, keepHs=False, defaultVal='N/A', skipProps=False, regName='molecules',
           skipSmiles=False, maxRowsCached=-1, uniqNames=False, addComputedProps=False,
           lazySupplier=False, numForPropScan=10, startAnew=True):
  if not lazySupplier:
    nMols = len(suppl)
  else:
    nMols = -1
  if not silent:
    logger.info(f"Generating molecular database in file {dbName}")
    if not lazySupplier:
      logger.info(f"  Processing {nMols} molecules")

  globalProps = {}
  if startAnew:
    if os.path.exists(dbName):
      for _ in range(5):
        try:
          os.unlink(dbName)
          break
        except Exception:
          import time
          time.sleep(2)
    if os.path.exists(dbName):
      raise IOError(f'could not delete old database {dbName}')
  sIter = iter(suppl)
  setattr(Compound, nameCol.lower(),
          Column(nameCol.lower(), String, default=defaultVal, unique=uniqNames))
  if not skipSmiles:
    Compound.smiles = Column(Text, unique=True)
  if not skipProps:
    while numForPropScan > 0:
      try:
        m = next(sIter)
      except StopIteration:
        numForPropScan = 0
        break
      if not m:
        continue
      for pn in m.GetPropNames():
        if pn.lower() == nameCol.lower():
          continue
        if pn not in globalProps:
          globalProps[pn] = 1
          setattr(Compound, pn.lower(), Column(pn.lower(), String, default=defaultVal))
      numForPropScan -= 1
    if addComputedProps:
      Compound.DonorCount = Column(Integer)
      Compound.AcceptorCount = Column(Integer)
      Compound.RotatableBondCount = Column(Integer)
      Compound.AMW = Column(Float)
      Compound.MolLogP = Column(Float)
  session = RegisterSchema(f'sqlite:///{dbName}')()

  nDone = 0
  cache = []
  for m in suppl:
    nDone += 1
    if not m:
      if errorsTo:
        if hasattr(suppl, 'GetItemText'):
          d = suppl.GetItemText(nDone - 1)
          errorsTo.write(d)
        else:
          logger.warning('full error file support not complete')
      continue

    cmpd = ProcessMol(session, m, globalProps, nDone, nameProp=nameProp, nameCol=nameCol,
                      redraw=redraw, keepHs=keepHs, skipProps=skipProps,
                      addComputedProps=addComputedProps, skipSmiles=skipSmiles)
    if cmpd is not None:
      cache.append(cmpd)

    if not silent and not nDone % 100:
      logger.info(f'  done {nDone}')
      try:
        session.commit()
      except Exception:
        session.rollback()
        for cmpd in cache:
          try:
            session.add(cmpd)
            session.commit()
          except Exception:
            session.rollback()
          except BaseException:
            # Rollback even with KeyboardInterrupt
            session.rollback()
            raise
      cache = []

  try:
    session.commit()
  except BaseException as exc:
    import traceback
    traceback.print_exc()
    session.rollback()
    for cmpd in cache:
      try:
        session.add(cmpd)
        session.commit()
      except Exception:
        session.rollback()
      except BaseException:
        session.rollback()
        raise
    if not isinstance(exc, Exception):
      # Re-raise on KeyboardInterrupt, SystemExit, etc.
      raise exc


if __name__ == '__main__':
  import sys
  sdf = Chem.SDMolSupplier(sys.argv[1])
  db = sys.argv[2]
  LoadDb(sdf, db, addComputedProps=False)
  session = RegisterSchema(f'sqlite:///{db}')()
  print('>>>>', len(session.query(Compound).all()))
