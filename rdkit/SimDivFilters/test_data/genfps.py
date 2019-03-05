

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Dbase import DbModule
from rdkit.Dbase.DbConnection import DbConnect
import pickle

if RDConfig.usePgSQL:
    dbName = "::RDTests"
else:
    dbName = "data.sqlt"

molTblName = 'simple_mols1'
fpTblName = 'simple_mols1_fp'
conn = DbConnect(dbName, molTblName)
conn.AddTable(fpTblName, 'id varchar(10),autofragmentfp %s' % DbModule.binaryTypeName)
d = conn.GetData()
for smi, ID in d:
    print(repr(ID), repr(smi))
    mol = Chem.MolFromSmiles(smi)
    fp = Chem.RDKFingerprint(mol)
    pkl = pickle.dumps(fp)
    conn.InsertData(fpTblName, (ID, DbModule.binaryHolder(pkl)))
conn.Commit()
