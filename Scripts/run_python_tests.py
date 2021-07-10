#
#  Copyright (C) 2018 Greg Landrum
#   All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import RDConfig
from rdkit import TestRunner
import os,time,sys

if __name__ == '__main__':
    script = 'test_list.py'
    os.chdir(RDConfig.RDCodeDir)
    t1 = time.time()
    failed,nTests = TestRunner.RunScript(script,doLongTests=False,
        verbose=True)
    t2 = time.time()
    TestRunner.ReportResults(script,failed,nTests,t2-t1,verbose=True,dest=sys.stderr)
    sys.exit(len(failed))
