#
#  Copyright (C) 2017 greg Landrum
#
#   @@ All Rights Reserved  @@
#
import os
import subprocess
import pytest

from rdkit import RDConfig


def test_Github1406():
    with open('data/simple.smi') as inf:
        p = subprocess.run(('python', 'rfrag.py'), stdin=inf, stdout=subprocess.PIPE)
    assert p.returncode == 0
    assert p.stdout == b'''c1ccccc1,benzene,,
Cc1ccccc1,toluene,,C[*:1].c1ccc(cc1)[*:1]
'''


if __name__ == '__main__':
    pytest.main()
