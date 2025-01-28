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

from indexing import cansmirk, heavy_atom_count
from rfrag import fragment_mol


def cmpd_not_in_db_mmp_query(in_smi, cmpd_id):

  query_contexts = set()
  cmpd_frags = fragment_mol(in_smi, cmpd_id)
  for row in cmpd_frags:
    row = row.rstrip()
    row_fields = re.split(',', row)
    if (row_fields[3].count(".") == 1):
      a, b = row_fields[3].split(".")
      query_contexts.add(a)
      query_contexts.add(b)
    else:
      query_contexts.add(row_fields[3])

  q_string = "','".join(query_contexts)
  q_string = "'%s'" % (q_string)

  query_sql = """
    select  c.cmpd_id,
            c.core_smi,
            con.context_smi,
            con.context_size
    from    core_table c, context_table con
    where   c.context_id in (select context_id from context_table where context_smi in (%s))
            and c.context_id = con.context_id""" % (q_string)
  cursor.execute(query_sql)
  results = cursor.fetchall()

  cmpd_size = heavy_atom_count(in_smi)
  print_smallest_change_mmp(results, cmpd_id, cmpd_size)


def run_mmp_query(cmpd_id, cmpd_size):
  query_sql = """
    select  c.cmpd_id,
            c.core_smi,
            con.context_smi,
            con.context_size
    from    core_table c, context_table con
    where   c.context_id in (select context_id from core_table where cmpd_id = '%s')
            and c.context_id = con.context_id""" % (cmpd_id)
  cursor.execute(query_sql)
  results = cursor.fetchall()

  print_smallest_change_mmp(results, cmpd_id, cmpd_size)


def print_smallest_change_mmp(db_results, cmpd_id, query_size):

  uniq_list = {}
  for r in db_results:
    if (r[0] != cmpd_id):
      #print r
      #for each unique compound keep the largest one in common
      if (r[0] not in uniq_list):
        uniq_list[r[0]] = r
      elif (r[3] > uniq_list[r[0]][3]):
        uniq_list[r[0]] = r

  for key, value in uniq_list.items():
    size_of_change = query_size - value[3]
    #print "q_size: %s, Size od change: %s, Ratio: %s" % (query_size,size_of_change,float(size_of_change)/query_size)
    if (use_ratio):
      if (float(size_of_change) / query_size <= ratio):
        cursor.execute("SELECT smiles FROM cmpd_smisp WHERE cmpd_id = ?", (key, ))
        rsmi = cursor.fetchone()[0]
        print("%s,%s,%s,%s,%s,%s" % (smi, rsmi, id, value[0], value[1], value[2]))
    elif (size_of_change <= max_size):
      cursor.execute("SELECT smiles FROM cmpd_smisp WHERE cmpd_id = ?", (key, ))
      rsmi = cursor.fetchone()[0]
      print("%s,%s,%s,%s,%s,%s" % (search_string, rsmi, id, value[0], value[1], value[2]))


def run_subs_query(subs):

  query_sql = """
    select  lhs_smi.smiles,
            lhs.cmpd_id,
            lhs.core_smi,
            rhs_smi.smiles,
            rhs.cmpd_id,
            rhs.core_smi,
            context_table.context_smi,
            rhs_smi.cmpd_size-context_table.context_size
    from    (select cmpd_id,core_smi,context_id from core_table where core_smi_ni = '%s') lhs,
            core_table rhs,
            cmpd_smisp lhs_smi,
            cmpd_smisp rhs_smi,
            context_table
    where   lhs.context_id = rhs.context_id
            and context_table.context_id = rhs.context_id
            and lhs_smi.cmpd_id = lhs.cmpd_id
            and rhs_smi.cmpd_id = rhs.cmpd_id
            and lhs.cmpd_id != rhs.cmpd_id
            and rhs_smi.cmpd_size-context_table.context_size <= %s""" % (subs, max_size)
  cursor.execute(query_sql)
  results = cursor.fetchall()

  for r in results:
    #make sure it is not the same core on both sides
    if (r[2] != r[5]):
      #cansmirk
      smirks, context = cansmirk(str(r[2]), str(r[5]), str(r[6]))
      if (have_id):
        print("%s,%s,%s,%s,%s,%s,%s,%s" % (subs, id, r[0], r[3], r[1], r[4], smirks, context))
      else:
        print("%s,%s,%s,%s,%s,%s,%s" % (subs, r[0], r[3], r[1], r[4], smirks, context))


def run_subs_smarts_query(subs_smarts):

  #set os environment for rdkit to use sqllite
  os.environ['RD_USESQLLITE'] = '1'
  temp_core_ni_file = 'temp_core_ni_file_%s' % (os.getpid())
  cmd = "python $RDBASE/Projects/DbCLI/SearchDb.py --dbDir=%s_smarts --smarts='%s' --silent >%s" % (
    pre, subs_smarts, temp_core_ni_file)
  subprocess.Popen(cmd, shell=True).wait()

  infile = open(temp_core_ni_file, 'r')
  for row in infile:
    row = row.rstrip()

    query_sql = """
        select  lhs_smi.smiles,
                lhs.cmpd_id,
                lhs.core_smi,
                rhs_smi.smiles,
                rhs.cmpd_id,
                rhs.core_smi,
                context_table.context_smi,
                rhs_smi.cmpd_size-context_table.context_size
        from    (select cmpd_id,core_smi,context_id from core_table where core_smi_ni = '%s') lhs,
                core_table rhs,
                cmpd_smisp lhs_smi,
                cmpd_smisp rhs_smi,
                context_table
        where   lhs.context_id = rhs.context_id
                and context_table.context_id = rhs.context_id
                and lhs_smi.cmpd_id = lhs.cmpd_id
                and rhs_smi.cmpd_id = rhs.cmpd_id
                and lhs.cmpd_id != rhs.cmpd_id
                and rhs_smi.cmpd_size-context_table.context_size <= %s
                and lhs_smi.cmpd_size-context_table.context_size <= %s""" % (row, max_size,
                                                                             max_size)
    cursor.execute(query_sql)
    results = cursor.fetchall()

    for r in results:
      #cansmirk
      smirks, context = cansmirk(str(r[2]), str(r[5]), str(r[6]))
      if (have_id):
        print("%s,%s,%s,%s,%s,%s,%s" % (id, r[0], r[3], r[1], r[4], smirks, context))
      else:
        print("%s,%s,%s,%s,%s,%s" % (r[0], r[3], r[1], r[4], smirks, context))
  infile.close()
  #remove temporary files
  os.unlink(temp_core_ni_file)


def run_trans_smarts_query(transform):

  lhs, rhs = transform.split(">>")
  matching_lhs = []
  matching_rhs = []

  #set os environment for rdkit to use sqllite
  os.environ['RD_USESQLLITE'] = '1'

  cmd = "python $RDBASE/Projects/DbCLI/SearchDb.py --dbDir=%s_smarts --smarts='%s' --silent" % (pre,
                                                                                                lhs)
  p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
  output = p1.communicate()[0].decode().rstrip()
  matching_lhs = output.split("\n")
  #sys.stderr.write("rhs: %s\n" % (len(matching_lhs)) )

  cmd = "python $RDBASE/Projects/DbCLI/SearchDb.py --dbDir=%s_smarts --smarts='%s' --silent" % (pre,
                                                                                                rhs)
  p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
  output = p1.communicate()[0].decode().rstrip()
  matching_rhs = output.split("\n")
  #sys.stderr.write("rhs: %s\n" % (len(matching_rhs)) )

  #sys.stderr.write('SQLlite method\n')
  lhs_q_string = "','".join(matching_lhs)
  lhs_q_string = "'%s'" % (lhs_q_string)

  rhs_q_string = "','".join(matching_rhs)
  rhs_q_string = "'%s'" % (rhs_q_string)

  query_sql = """
    select  lhs_smi.smiles,
            lhs.cmpd_id,
            lhs.core_smi,
            rhs_smi.smiles,
            rhs.cmpd_id,
            rhs.core_smi,
            context_table.context_smi
    from    (select cmpd_id,core_smi,context_id from core_table where core_smi_ni in (%s) ) lhs,
            (select cmpd_id,core_smi,context_id from core_table where core_smi_ni in (%s) ) rhs,
            cmpd_smisp lhs_smi,
            cmpd_smisp rhs_smi,
            context_table
    where   lhs.context_id = rhs.context_id
            and context_table.context_id = rhs.context_id
            and lhs_smi.cmpd_id = lhs.cmpd_id
            and rhs_smi.cmpd_id = rhs.cmpd_id
            and lhs.cmpd_id != rhs.cmpd_id
            and rhs_smi.cmpd_size-context_table.context_size <= %s
            and lhs_smi.cmpd_size-context_table.context_size <= %s """ % (
    lhs_q_string, rhs_q_string, max_size, max_size)

  cursor.execute(query_sql)
  results = cursor.fetchall()

  for r in results:
    smirks, context = cansmirk(str(r[2]), str(r[5]), str(r[6]))
    if (have_id):
      print("%s,%s,%s,%s,%s,%s,%s,%s" % (transform, id, r[0], r[3], r[1], r[4], smirks, context))
    else:
      print("%s,%s,%s,%s,%s,%s,%s" % (transform, r[0], r[3], r[1], r[4], smirks, context))


def run_trans_query(transform):

  lhs, rhs = transform.split(">>")

  #remove connectivity info
  lhs_ni = remove_numbers(lhs)
  rhs_ni = remove_numbers(rhs)

  query_sql = """
    select  lhs_smi.smiles,
            lhs.cmpd_id,
            lhs.core_smi,
            rhs_smi.smiles,
            rhs.cmpd_id,
            rhs.core_smi,
            context_table.context_smi
    from    (select cmpd_id,core_smi,context_id from core_table where core_smi_ni = '%s') lhs,
            (select cmpd_id,core_smi,context_id from core_table where core_smi_ni = '%s') rhs,
            cmpd_smisp lhs_smi,
            cmpd_smisp rhs_smi,
            context_table
    where   lhs.context_id = rhs.context_id
            and context_table.context_id = rhs.context_id
            and lhs_smi.cmpd_id = lhs.cmpd_id
            and rhs_smi.cmpd_id = rhs.cmpd_id""" % (lhs_ni, rhs_ni)

  cursor.execute(query_sql)
  results = cursor.fetchall()

  for r in results:
    smirks, context = cansmirk(str(r[2]), str(r[5]), str(r[6]))
    #make sure connectivity is correct
    if (smirks == transform):
      if (have_id):
        print("%s,%s,%s,%s,%s,%s,%s" % (id, r[0], r[3], r[1], r[4], smirks, context))
      else:
        print("%s,%s,%s,%s,%s,%s" % (r[0], r[3], r[1], r[4], smirks, context))


def remove_numbers(in_string):

  out_string = re.sub(r'\[\*\:1\]', '[*]', in_string)
  out_string = re.sub(r'\[\*\:2\]', '[*]', out_string)
  out_string = re.sub(r'\[\*\:3\]', '[*]', out_string)

  return out_string


#quick class to help with the formatting of optparse
class MyParser(OptionParser):

  def format_description(self, formatter):
    return self.description


parser = MyParser(
  description="""Program to search MMP db. The types of searching that can be performed are as
follows:

mmp: Find all MMPs of a input/query compound to the compounds in the db

subs: Find all MMPs in the db where the LHS of the transform matches an input
substructure. Make sure the attached points are donated by an asterisk and the
input substructure has been canonicalised (eg. [*]c1ccccc1).

trans: Find all MMPs that match the input transform/SMIRKS. Make sure the input
SMIRKS has been canonicalised using the cansmirk.py program.

subs_smarts: Find all MMPs in the db where the LHS of the transform matches an
input SMARTS. The attachment points in the SMARTS can be donated by [#0] (eg.
[#0]c1ccccc1).

trans_smarts: Find all MMPs that match the LHS and RHS SMARTS of the input
transform. The transform SMARTS are input as LHS_SMARTS>>RHS_SMARTS (eg.
[#0]c1ccccc1>>[#0]c1ccncc1). Note: This search can take a long time to run if a
very general SMARTS expression is used.
""")
parser.add_option(
  '-t', '--type', action='store', dest='type', type='string',
  help='Type of search required. Options are: mmp, subs, trans, subs_smarts, trans_smarts')
parser.add_option(
  '-m', '--maxsize', action='store', dest='maxsize', type='int', help=
  'Maximum size of change (in heavy atoms) allowed in matched molecular pairs identified. DEFAULT=10. \
                  Note: This option overrides the ratio option if both are specified.')
parser.add_option(
  '-r', '--ratio', action='store', dest='ratio', type='float', help=
  'Only applicable with the mmp search type. Maximum ratio of change allowed in matched molecular pairs identified. The ratio is: size of change / \
                  size of cmpd (in terms of heavy atoms) for the QUERY MOLECULE. DEFAULT=0.3. Note: If this option is used with the maxsize option, the maxsize option will be used.'
)
parser.add_option('-p', '--prefix', action='store', dest='prefix', type='string',
                  help='Prefix for the db file. DEFAULT=mmp')

#parse the command line options
(options, args) = parser.parse_args()

#note max heavy atom count does not
#include the attachment points (*)
max_size = 10
ratio = 0.3
use_ratio = False
have_id = True

search_type = "mmp"
db_name = "mmp.db"
pre = "mmp"

if options.maxsize is not None:
  max_size = options.maxsize
elif options.ratio is not None:
  ratio = options.ratio
  if (ratio >= 1):
    print("Ratio specified: %s. Ratio needs to be less than 1.")
    sys.exit(1)
  use_ratio = True

if options.type is not None:
  if ((options.type == "mmp") or (options.type == "subs") or (options.type == "trans")
      or (options.type == "subs_smarts") or (options.type == "trans_smarts")):
    search_type = options.type
  else:
    print(
      "Unrecognised search type. Please choose from: mmp, subs, trans, subs_smarts, trans_smarts")
    sys.exit(1)
else:
  print(
    "Please specify search type. Please choose from: mmp, subs, trans, subs_smarts, trans_smarts")
  sys.exit(1)

if options.prefix is not None:
  pre = options.prefix
  db_name = "%s.db" % (pre)

#connect to db
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

#read the STDIN
for line in sys.stdin:

  line = line.rstrip()
  line_fields = re.split(r'\s|,', line)

  if (len(line_fields) == 1):
    id = line_fields[0]
    have_id = False
  else:
    id = line_fields[1]

  search_string = line_fields[0]

  if (search_type == "mmp"):
    #check smiles is in the database
    cursor.execute("SELECT cmpd_id,cmpd_size FROM cmpd_smisp WHERE smiles = ?", (search_string, ))
    d_res = cursor.fetchone()

    #cmpd in the db
    if (d_res):
      id_in_db, query_size = d_res
      run_mmp_query(id_in_db, query_size)
    else:
      #print "Not in db"
      cmpd_not_in_db_mmp_query(search_string, id)

  #if doing a subs query
  elif (search_type == "subs"):
    run_subs_query(search_string)

  elif (search_type == "trans"):
    run_trans_query(search_string)

  #smarts queries
  elif (search_type == "subs_smarts"):
    run_subs_smarts_query(search_string)

  elif (search_type == "trans_smarts"):
    run_trans_smarts_query(search_string)
