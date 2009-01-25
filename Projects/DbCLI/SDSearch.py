# $Id$
#
#  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Greg Landrum Nov 2006
#
_version = "0.5.2"

_usage="""
  SDSearch [optional arguments] <sdfilename>


  NOTES:

    - the property names for the database are the union of those for
      all molecules

    - missing property values will be set to 'N/A'
    
    - The property names may be altered on loading the database.  Any
      non-alphanumeric character in a property name will be replaced
      with '_'. e.g. "Gold.Goldscore.Constraint.Score" becomes
      "Gold_Goldscore_Constraint_Score".  This is important to know
      when querying.

    - Property names are not case sensitive in the database; this may
      cause some problems if they are case sensitive in the sd file
      
"""
from pyRDKit import Chem
from pyRDKit.Chem import AllChem
from pyRDKit.Dbase.DbConnection import DbConnect
from pyRDKit.Dbase import DbModule
from pyRDKit.Chem.MolDb import Loader
import sys,os,time
from pyRDKit import RDConfig

# set up the logger:
from pyRDKit import RDLogger as logging
logger = logging.logger()
logger.setLevel(logging.INFO)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --- --- ---  --- --- ---  --- --- --- 
# set up the option parser
from optparse import OptionParser
parser=OptionParser(_usage,version='%prog '+_version)
parser.add_option('--smarts','--sma',default='',
                  help='provide the SMILES or SMARTS to be used as a substructure query')
parser.add_option('-q','--query',default='',
                  help='provide the property query (see the NOTE about property names')
parser.add_option('--sdfOut','--sdOut',default='',
                  help="""provide the name of an SD file to be created
                        with matching molecules. The default is not to
                        produce an SD file Providing the filename "-"
                        will result in the output going to stdout (the
                        console).""")
parser.add_option('-n','--negate',default=False,action='store_true',
                  help="negate the smarts query (e.g. only print molecules that do *not* match the query)")
parser.add_option('--forceLoad',default=False,action='store_true',
                  help="force the database to be regenerated")
parser.add_option('--noNames',default=False,action='store_true',
                  help="do not print out the names of matching molecules")
parser.add_option('--nameOut',default='',
                  help="provide the name of a text file to store the names of matching molecules. The default is to print the names to the console.")
parser.add_option('--nameCol',default='_Name',
                  help="provide an alternative property to be used as the molecule name. The default is to use the mol-block name")
parser.add_option('--dbName',default='',
                  help="""provide the name of the database to be created and/or
                        used.  The default is to remove '.sdf' from
                        the name of the input SD file and appending
                        '.sqlt' in its place.  e.g. mols.sdf ->
                        mols.sqlt
                        If the dbname is provide and the database is 
                        already populated, the sd filename is optional""") 
parser.add_option('--redraw',default=False,action='store_true',
                  help='generate new 2D coordinates for the molecules')
parser.add_option('--errorFile','--errs',default='',
                  help="""provide the name of an sd file to store molecules which cannot be parsed. The default is to just ignore those molecules.""")
parser.add_option('--keepHs',default=False,action='store_true',
                  help='do not remove Hs from the molecules before storing them in the db (NOTE: this is a bit fragile)')
parser.add_option('--silent',default=False,action='store_true',
                  help='switch off the printing of status messages to the console.')


if __name__=='__main__':
  options,args = parser.parse_args()
  if len(args)!=1 and not options.dbName:
    parser.error('please provide either an sd file or a db name')

  if len(args):
    sdName = args[0]
    if not options.dbName:
      options.dbName = sdName.split('.sdf')[0]+'.sqlt'
    
  if options.sdfOut=='-':
    # disable info messages, which normally go to stdout:
    logger.setLevel(logging.WARNING)
    options.silent=True
  
  if options.forceLoad or not os.path.exists(options.dbName):
    if options.errorFile: 
      options.errorFile=open(options.errorFile,'w+')
    suppl = Chem.SDMolSupplier(sdName,sanitize=not options.keepHs)
    Loader.LoadDb(suppl,options.dbName,nameProp=options.nameCol,redraw=options.redraw,
                  errorsTo=options.errorFile,keepHs=options.keepHs)
    options.redraw=False
    

  patt=None
  if options.smarts:
    try:
      patt = Chem.MolFromSmarts(options.smarts)
    except:
      logger.error('Could not parse smarts',exc_info=True)
      sys.exit(2)

  if options.sdfOut and options.sdfOut!='-':
    sdW = file(options.sdfOut,'w+')
  elif options.sdfOut=='-':
    sdW = sys.stdout
  else:
    sdW = None

  if not options.noNames and options.nameOut:
    outF = file(options.nameOut,'w+')
  else:
    outF = sys.stdout
  if not options.silent: logger.info("Searching database")
  t1 = time.time()
  conn = DbConnect(options.dbName)
  cns = conn.GetColumnNames('molecules')
  curs = conn.GetCursor()
  if options.query:
    curs.execute('select * from molecules where %s'%options.query)
  else:
    curs.execute('select * from molecules')
  nHits = 0
  hits = []
  sz=0
  row=curs.fetchone()
  while row:
    sz += 1
    if not options.silent and not sz%1000: logger.info('searched %d molecules, found %d hits'%(sz,nHits))
    if patt or sdW:
      m = Chem.Mol(str(row[-1]))
      m.SetProp('_Name',str(row[0]))
    else:
      m = None

    if patt:
      match = m.HasSubstructMatch(patt)
      if options.negate:
        match=not match
      if not match:
        row = curs.fetchone()
        continue
      if options.redraw:
        AllChem.Compute2DCoords(m)
    if sdW:
      print >>sdW,Chem.MolToMolBlock(m)
    if sdW:
      for i in range(1,len(cns)-2):
        pn = cns[i]
        pv = str(row[i])
        if sdW:
          print >>sdW,'> <%s>\n%s\n'%(pn,pv)
      if sdW:
        print >>sdW,'$$$$'  
    if not options.noNames:
      print >>outF,row[0]
    row = curs.fetchone()
    hits.append(m)
    nHits+=1
  t2 = time.time()
  if not options.silent: 
    logger.info("%d matches found in %d molecules."%(nHits,sz))
    logger.info("\tThe search took %.2f seconds."%(t2-t1))

      
  
