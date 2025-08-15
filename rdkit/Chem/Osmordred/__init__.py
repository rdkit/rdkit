# Copyright (c) 2020-2023 Greg Landrum and other RDKit contributors
#  All rights reserved.
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
"""
Osmordred descriptors for RDKit
"""
import sys

try:
    from . osmordred import CalcOsmordred, Calculate, GetDescriptorNames, CalcOsmordredFromMol
    HasOrmordred = True
except:
    sys.stderr.write("Osmordred descriptors not available, please ensure rdOsmordred is built")
    def CalcOsmordred(*a, **kw):
        return []
    def Calculate(*a, **kw):
        return []
    def GetDescriptorNames():
        return []
    def CalcOsmordredFromMol(*a, **kw):
        return []
    HasOrmordred = False
