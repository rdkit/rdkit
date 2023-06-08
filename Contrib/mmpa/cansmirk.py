# Copyright (c) 2012, GlaxoSmithKline Research & Development Ltd.
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
# Created by Jameed Hussain, September 2012

import re
import sys

from indexing import cansmirk

from rdkit import Chem

if __name__ == '__main__':

  if (len(sys.argv) >= 2):
    print(
      "Program that canonicalises an input SMIRKS so its in same format as MMP identification program.\n"
    )
    print("USAGE: ./cansmirks.py <file_of_smirks\n")
    sys.exit(1)

  #read the STDIN
  for line in sys.stdin:

    line = line.rstrip()

    line_fields = re.split(r'\s|,', line)
    smirks = line_fields[0]

    if (len(line_fields) == 1):
      id = ""
    else:
      id = line_fields[1]

    lhs, rhs = smirks.split(">>")

    l = Chem.MolFromSmiles(lhs)
    if l is None:
      sys.stderr.write("Can't generate mol for: %s\n" % (lhs))
      continue

    r = Chem.MolFromSmiles(rhs)
    if r is None:
      sys.stderr.write("Can't generate mol for: %s\n" % (rhs))
      continue

    clhs = Chem.MolToSmiles(l, isomericSmiles=True)
    crhs = Chem.MolToSmiles(r, isomericSmiles=True)

    #just need to take care of [*H:1]
    if (clhs == '[*H:1]'):
      clhs = '[*:1][H]'

    if (crhs == '[*H:1]'):
      crhs = '[*:1][H]'

    #print clhs
    #print crhs

    csmirk, context = cansmirk(clhs, crhs, "")

    print("%s %s" % (csmirk, id))
