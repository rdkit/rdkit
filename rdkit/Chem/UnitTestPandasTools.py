from __future__ import print_function
import os

from  rdkit.Chem import PandasTools

class TestLoadSDF(unittest.TestCase):
    def setUp(self):
        self.gz_filename = os.path.join(RDConfig.RDCodeDir, 'Chem/test_data', 'pandas_load.sdf.gz')
        
    # the doctest tests loading from a ".sdf" file
    
    def test_load_gzip_file(self):
        df = PandasTools.LoadSDF(self.gz_filename)
        self.assertEqual(len(df), 
