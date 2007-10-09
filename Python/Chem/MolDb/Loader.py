# $Id$
#
#  Copyright (C) 2007 Greg Landrum
#   @@ All Rights Reserved @@
#
import Chem
from Chem import AllChem
from Dbase.DbConnection import DbConnect
from Dbase import DbModule

#set up the logger:
import RDLogger as logging
logger = logging.logger()
logger.setLevel(logging.INFO)

def LoadDb(suppl,dbName,nameProp='_Name',nameCol='compound_id',silent=False,
           redraw=False,errorsTo=None,keepHs=False,defaultVal='N/A',
           regName='molecules'):
  nMols = len(suppl)
  if not silent:
    logger.info("Generating molecular database in file %s"%dbName)
    logger.info("  Processing %d molecules"%nMols)
  rows = []
  globalProps = {}
  
  nDone = 0
  typeConversions={0:('varchar',str),1:('float',float),2:('int',int)}
  for m in suppl:
    nDone +=1
    if not m:
      if errorsTo:
        if hasattr(suppl,'GetItemText'):
          d = suppl.GetItemText(nDone-1)
          errorsTo.write(d)
        else:
          logger.warning('full error file support not complete')
      continue
    if keepHs:
      Chem.SanitizeMol(m)
    pns = list(m.GetPropNames())
    try:
      nm = m.GetProp(nameProp)
    except KeyError:
      nm = None
    if not nm:
      nm = 'Mol_%d'%nDone
    row = [nm]
    pD={}
    for pi,pn in enumerate(pns):
      if pn==nameCol: continue
      pv = m.GetProp(pn).strip()
      
      if pv.find('>')<0 and pv.find('<')<0:
        colTyp = globalProps.get(pn,2)
        while colTyp>0:
          try:
            tpi = typeConversions[colTyp][1](pv)
          except:
            colTyp-=1
          else:
            break
        globalProps[pn]=colTyp
        pD[pn]=typeConversions[colTyp][1](pv)
      else:
        pD[pn]=pv
    if redraw:
      AllChem.Compute2DCoords(m)
    row.append(Chem.MolToSmiles(m,True))
    row.append(DbModule.binaryHolder(m.ToBinary()))
    row.append(pD)
    rows.append(row)
    if not silent and not nDone%100:
      logger.info('  done %d'%nDone)

  typs = ['%s varchar not null'%nameCol]
  pns = []
  for pn,v in globalProps.iteritems():
    addNm = pn
    addNm = addNm.replace('.','_')
    addNm = addNm.replace('-','_')
    addNm = addNm.replace(' ','_')
    addNm = addNm.replace(':','_')
    typs.append('%s %s'%(addNm,typeConversions[v][0]))
  typs.append('smiles varchar')
  typs.append('molpkl %s'%(DbModule.binaryTypeName))

  for i,row in enumerate(rows):
    newRow = []
    newRow.append(row[0])
    pD=row[3]
    for pn in globalProps:
      pv = pD.get(pn,defaultVal)
      newRow.append(pv)
    newRow.append(row[1])
    newRow.append(row[2])
    rows[i] = newRow

  conn = DbConnect(dbName)
  curs = conn.GetCursor()
  try:
    curs.execute('drop table %s'%regName)
  except:
    pass
  curs.execute('create table %s (%s)'%(regName,','.join(typs)))
  qs = ','.join([DbModule.placeHolder for x in typs])
  curs.executemany('insert into %s values (%s)'%(regName,qs),rows)
  
  conn.Commit()

