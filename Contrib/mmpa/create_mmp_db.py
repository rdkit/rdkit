# Copyright (c) 2013, GlaxoSmithKline Research & Development Ltd.
# All rights reserved.
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
#     * Neither the name of GlaxoSmithKline Research & Development Ltd.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
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
# Created by Jameed Hussain, July 2013

import os
import re
import sqlite3
import subprocess
import sys
from optparse import OptionParser

from rdkit import Chem


def heavy_atom_count(smi):
  m = Chem.MolFromSmiles(smi)
  return m.GetNumAtoms()


def add_to_db(id, core, context, context_size):

  cursor.execute("INSERT OR IGNORE INTO context_table(context_smi) VALUES(?)", (context, ))

  the_id_of_the_row = None
  cursor.execute("SELECT context_id FROM context_table WHERE context_smi = ?", (context, ))
  the_id_of_the_row = cursor.fetchone()[0]

  core_ni = re.sub(r'\[\*\:1\]', '[*]', core)
  core_ni = re.sub(r'\[\*\:2\]', '[*]', core_ni)
  core_ni = re.sub(r'\[\*\:3\]', '[*]', core_ni)

  cursor.execute(
    "INSERT INTO core_table(context_id, cmpd_id, core_smi, core_smi_ni) VALUES(?,?,?,?);",
    (the_id_of_the_row, id, core, core_ni))
  #print "the_id_of_the_row:%s, id:%s, core:%s" % (the_id_of_the_row,id,core)

  #check if heavy atom count for context_id has been calculated
  cursor.execute("select context_size from context_table where context_id = ?",
                 (the_id_of_the_row, ))
  heavy_calculated = cursor.fetchone()
  #add heavy atom count if needed
  if (heavy_calculated[0] is None):
    cursor.execute("update context_table set context_size = ? where context_id = ? ",
                   (context_size, the_id_of_the_row))


def index_hydrogen_change():
  #Algorithm details
  #1) Loop through context smiles in the db
  #2) If context smiles is the result of a single cut (so contains only 1 *) replace the * with H, and cansmi
  #3) If full smiles matches this new context smi, add *-H to that core_smi index.

  cursor.execute("SELECT context_smi FROM context_table")
  res = cursor.fetchall()

  for row in res:
    key = row[0]
    attachments = key.count('*')
    #print attachments

    if (attachments == 1):
      smi = str(key)
      #simple method
      smi = re.sub(r'\[\*\:1\]', '[H]', smi)

      #now cansmi it
      temp = Chem.MolFromSmiles(smi)

      if temp is None:
        sys.stderr.write('Error with key: %s, Added H: %s\n' % (key, smi))
      else:
        c_smi = Chem.MolToSmiles(temp, isomericSmiles=True)

        cursor.execute("SELECT cmpd_id FROM cmpd_smisp WHERE smiles = ?", (c_smi, ))
        cmpd_id = cursor.fetchone()
        if (cmpd_id):
          core = "[*:1][H]"
          cmpd_id = cmpd_id[0]
          key_size = temp.GetNumAtoms() - 1
          add_to_db(cmpd_id, core, key, key_size)
          #print "Added: id:%s, core:%s, context:%s" % (str(cmpd_id),core,key)

        #set up the command line options


parser = OptionParser(description="Program to create an MMP db. ")
parser.add_option(
  '-p', '--prefix', action='store', dest='prefix', type='string',
  help='Prefix to use for the db file (and directory for SMARTS index). DEFAULT=mmp')
parser.add_option(
  '-m', '--maxsize', action='store', dest='maxsize', type='int',
  help="Maximum size of change (in heavy atoms) that is stored in the database. DEFAULT=15. \t\
                  Note: Any MMPs that involve a change greater than this value will not be stored in the database and hence not be identified in the searching."
)
parser.add_option(
  '-s', '--smarts', default=False, action='store_true', dest='sma', help=
  'Build SMARTS db so can perform SMARTS searching against db. Note: Will make the build process somewhat slower.'
)
#parse the command line options
(options, args) = parser.parse_args()

#this is the heavy atom limit
#of the size of change that goes in the db
max_size = 15
db_name = "mmp.db"
pre = "mmp"

#print options
if options.maxsize is not None:
  max_size = options.maxsize

if options.prefix is not None:
  pre = options.prefix
  db_name = "%s.db" % (pre)

con = sqlite3.connect(db_name)
cursor = con.cursor()

#these setting increase performance
cursor.execute('PRAGMA main.page_size = 4096;')
cursor.execute('PRAGMA main.cache_size=10000;')
cursor.execute('PRAGMA main.locking_mode=EXCLUSIVE;')
cursor.execute('PRAGMA main.synchronous=NORMAL;')
cursor.execute('PRAGMA main.journal_mode=WAL;')
cursor.execute('PRAGMA main.cache_size=5000;')
cursor.execute('PRAGMA main.temp_store = MEMORY;')

#tables
cursor.execute("DROP TABLE IF EXISTS context_table")
cursor.execute("""CREATE TABLE context_table
                  (
                    context_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    context_smi VARCHAR(1000) NOT NULL UNIQUE,
                    context_size INTEGER
                  )""")

cursor.execute("DROP TABLE IF EXISTS core_table")
cursor.execute("""CREATE TABLE core_table
                  (
                    context_id INTEGER NOT NULL,
                    cmpd_id VARCHAR(20) NOT NULL,
                    core_smi VARCHAR(1000) NOT NULL,
                    core_smi_ni VARCHAR(1000) NOT NULL
                  )""")

cursor.execute("DROP TABLE IF EXISTS cmpd_smisp")
cursor.execute("""CREATE TABLE cmpd_smisp
                  (
                    cmpd_id VARCHAR(20) NOT NULL UNIQUE,
                    smiles VARCHAR(1000),
                    cmpd_size INTEGER
                  )""")

for line in sys.stdin:

  line = line.rstrip()
  smi, id, core, context = line.split(',')

  cursor.execute("INSERT OR IGNORE INTO cmpd_smisp(smiles,cmpd_id) VALUES(?,?)", (smi, id))

  #check if heavy atom count for compound has been calculated
  cursor.execute("select cmpd_size from cmpd_smisp where cmpd_id = ?", (id, ))
  heavy_calculated = cursor.fetchone()
  #add heavy atom count if needed
  if (heavy_calculated[0] is None):
    cmpd_heavy_atoms = heavy_atom_count(smi)
    cursor.execute("update cmpd_smisp set cmpd_size = ? where cmpd_id = ? ", (cmpd_heavy_atoms, id))
  else:
    cmpd_heavy_atoms = heavy_calculated[0]

  #deal with cmpds that have not been core_smied
  if (len(core) == 0) and (len(context) == 0):
    continue

  #deal with single cuts
  if (len(core) == 0):
    side_chains = context.split('.')

    context = side_chains[0]
    core = side_chains[1]
    context_size = heavy_atom_count(context) - context.count("*")

    #limit what goes in db based on max_size
    if (cmpd_heavy_atoms - context_size <= max_size):
      #put in db
      add_to_db(id, core, context, context_size)

    context = side_chains[1]
    core = side_chains[0]
    context_size = heavy_atom_count(context) - context.count("*")

    #limit what goes in db based on max_size
    if (cmpd_heavy_atoms - context_size <= max_size):
      #put in db
      add_to_db(id, core, context, context_size)

  #double or triple cut
  else:
    context_size = heavy_atom_count(context) - context.count("*")
    if (cmpd_heavy_atoms - context_size <= max_size):
      #put in db
      add_to_db(id, core, context, context_size)

#index the H change
index_hydrogen_change()

#Note: Indices automatically built for primary key and unique columns
#so need to build the remaining indices that are needed
cursor.execute("DROP INDEX IF EXISTS context_smi_idx")
cursor.execute("CREATE INDEX context_smi_idx ON context_table (context_smi)")
cursor.execute("DROP INDEX IF EXISTS context_id_idx")
cursor.execute("CREATE INDEX context_id_idx ON core_table (context_id)")
cursor.execute("DROP INDEX IF EXISTS cmpd_id_idx")
cursor.execute("CREATE INDEX cmpd_id_idx ON core_table (cmpd_id)")
cursor.execute("DROP INDEX IF EXISTS core_smi_idx")
cursor.execute("CREATE INDEX core_smi_ni_idx ON core_table (core_smi_ni)")
cursor.execute("DROP INDEX IF EXISTS smiles_idx")
cursor.execute("CREATE INDEX smiles_idx ON cmpd_smisp (smiles)")

#add smarts searching capability if required
if (options.sma):
  #first dump out the distinct core_smi_ni from core_table
  #and write to a file
  temp_core_ni_file = 'temp_core_ni_file_%s' % (os.getpid())
  outfile = open(temp_core_ni_file, 'w')

  query_sql = """select distinct(core_smi_ni) from core_table"""

  cursor.execute(query_sql)
  results = cursor.fetchall()

  for r in results:
    #print "%s %s" % (r[0],r[0])
    outfile.write("%s %s\n" % (r[0], r[0]))

  outfile.close()

  #set os environment for rdkit to use sqllite
  os.environ['RD_USESQLLITE'] = '1'
  #use the DbCli utility in RDKit: http://code.google.com/p/rdkit/wiki/UsingTheDbCLI
  cmd = 'python $RDBASE/Projects/DbCLI/CreateDb.py --dbDir=%s_smarts --molFormat=smiles %s --noPairs --noFingerprints --noDescriptors --noProps --noMorganFps --noSmiles --silent' % (
    pre, temp_core_ni_file)
  subprocess.Popen(cmd, shell=True).wait()
  #remove temporary file
  os.unlink(temp_core_ni_file)
