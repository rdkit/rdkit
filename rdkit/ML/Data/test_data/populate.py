import RDConfig
from Dbase import DbUtils
from io import StringIO

basic_2class = """ID,VAL,ACT
id-1,1.0,0
id-2,23.0,1
id-3,4.0,0
id-4,6.0,1
id-5,321.0,0
id-6,885.0,1
id-7,252.0,0
id-8,1351.1,1
id-9,3215.0,0
id-10,2585.0,1
id-11,55.0,0
id-12,12.0,1
"""
float_2class = """ID,VAL,ACT
id-1,1.0,-0.1
id-2,23.0,1.2
id-3,4.0,0.5
id-4,6.0,1.01
id-5,321.0,0.7
id-6,885.0,2.0
id-7,252.0,0.0
id-8,1351.1,1.1
id-9,3215.0,-0.9
id-10,2585.0,1.3
id-11,55.0,0.9
id-12,12.0,3.0
"""
io = StringIO(basic_2class)
DbUtils.TextFileToDatabase(RDConfig.RDTestDatabase, 'basic_2class', io)
io = StringIO(float_2class)
DbUtils.TextFileToDatabase(RDConfig.RDTestDatabase, 'float_2class', io)
