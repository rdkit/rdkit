#
#  Copyright (C) 2000  greg Landrum
#
""" utility functions with "chemical know-how"

"""

import os
import re

from rdkit import RDConfig

if not RDConfig.usePgSQL:
  _atomDbName = os.path.join(RDConfig.RDDataDir, 'atomdb.gdb')
else:
  _atomDbName = "::RDData"


def GetAtomicData(atomDict, descriptorsDesired, dBase=_atomDbName, table='atomic_data', where='',
                  user='sysdba', password='masterkey', includeElCounts=0):
  """ pulls atomic data from a database

      **Arguments**

        - atomDict: the dictionary to populate

        - descriptorsDesired: the descriptors to pull for each atom

        - dBase: the DB to use

        - table: the DB table to use

        - where: the SQL where clause

        - user: the user name to use with the DB

        - password: the password to use with the DB

        - includeElCounts: if nonzero, valence electron count fields are added to
           the _atomDict_

    """
  extraFields = ['NVAL', 'NVAL_NO_FULL_F', 'NVAL_NO_FULL_D', 'NVAL_NO_FULL']
  from rdkit.Dbase import DbModule
  cn = DbModule.connect(dBase, user, password)
  c = cn.cursor()
  descriptorsDesired = [s.upper() for s in descriptorsDesired]
  if 'NAME' not in descriptorsDesired:
    descriptorsDesired.append('NAME')
  if includeElCounts and 'CONFIG' not in descriptorsDesired:
    descriptorsDesired.append('CONFIG')
  for field in extraFields:
    if field in descriptorsDesired:
      descriptorsDesired.remove(field)
  toPull = ','.join(descriptorsDesired)
  command = 'select %s from atomic_data %s' % (toPull, where)
  try:
    c.execute(command)
  except Exception:
    print('Problems executing command:', command)
    return
  res = c.fetchall()
  for atom in res:
    tDict = {}
    for i in range(len(descriptorsDesired)):
      desc = descriptorsDesired[i]
      val = atom[i]
      tDict[desc] = val
    name = tDict['NAME']
    atomDict[name] = tDict
    if includeElCounts:
      config = atomDict[name]['CONFIG']
      atomDict[name]['NVAL'] = ConfigToNumElectrons(config)
      atomDict[name]['NVAL_NO_FULL_F'] = ConfigToNumElectrons(config, ignoreFullF=1)
      atomDict[name]['NVAL_NO_FULL_D'] = ConfigToNumElectrons(config, ignoreFullD=1)
      atomDict[name]['NVAL_NO_FULL'] = ConfigToNumElectrons(config, ignoreFullF=1, ignoreFullD=1)


def SplitComposition(compStr):
  """ Takes a simple chemical composition and turns into a list of element,# pairs.

        i.e. 'Fe3Al' -> [('Fe',3),('Al',1)]

        **Arguments**

         - compStr: the composition string to be processed

        **Returns**

         - the *composVect* corresponding to _compStr_

        **Note**

          -this isn't smart enough by half to deal with anything even
              remotely subtle, so be gentle.

    """
  target = r'([A-Z][a-z]?)([0-9\.]*)'

  theExpr = re.compile(target)

  matches = theExpr.findall(compStr)
  res = []
  for match in matches:
    if len(match[1]) > 0:
      res.append((match[0], float(match[1])))
    else:
      res.append((match[0], 1))

  return res


def ConfigToNumElectrons(config, ignoreFullD=0, ignoreFullF=0):
  """ counts the number of electrons appearing in a configuration string

      **Arguments**

        - config: the configuration string (e.g. '2s^2 2p^4')

        - ignoreFullD: toggles not counting full d shells

        - ignoreFullF: toggles not counting full f shells

      **Returns**

        the number of valence electrons

    """
  arr = config.split(' ')

  nEl = 0
  for i in range(1, len(arr)):
    l = arr[i].split('^')
    incr = int(l[1])
    if ignoreFullF and incr == 14 and l[0].find('f') != -1 and len(arr) > 2:
      incr = 0
    if ignoreFullD and incr == 10 and l[0].find('d') != -1 and len(arr) > 2:
      incr = 0
    nEl = nEl + incr
  return nEl


if __name__ == '__main__':  # pragma: nocover

  print(SplitComposition('Fe'))
  print(SplitComposition('Fe3Al'))
  print(SplitComposition('Fe99PdAl'))
  print(SplitComposition('TiNiSiSO12P'))
  temp = [
    '[Xe] 4f^12 6s^2', '[Xe] 4f^14 5d^6 6s^2', '[Xe] 4f^14 5d^10 6s^2',
    '[Xe] 4f^14 5d^10 6s^2 6p^1', '[Xe] 5d^10'
  ]
  print('ignore all')
  for entry in temp:
    print(entry, '\t\t\t\t', ConfigToNumElectrons(entry, ignoreFullD=1, ignoreFullF=1))
  print('ignore d')
  for entry in temp:
    print(entry, '\t\t\t\t', ConfigToNumElectrons(entry, ignoreFullD=1, ignoreFullF=0))
  print('ignore f')
  for entry in temp:
    print(entry, '\t\t\t\t', ConfigToNumElectrons(entry, ignoreFullD=0, ignoreFullF=1))
  print('ignore None')
  for entry in temp:
    print(entry, '\t\t\t\t', ConfigToNumElectrons(entry, ignoreFullD=0, ignoreFullF=0))
