# $Id$
#
# Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
from pyRDKit import RDConfig

if hasattr(RDConfig,"usePgSQL") and RDConfig.usePgSQL:
  useGvib=0
  from pyPgSQL import PgSQL
  # as of this writing (March 2004), this results in a speedup in
  # getting results back from the wrapper:
  PgSQL.fetchReturnsList=1

  from pyPgSQL.PgSQL import *
  sqlTextTypes = [PG_CHAR,PG_BPCHAR,PG_TEXT,PG_VARCHAR,PG_NAME]
  sqlIntTypes = [PG_INT8,PG_INT2,PG_INT4]
  sqlFloatTypes = [PG_FLOAT4,PG_FLOAT8]
  sqlBinTypes = [PG_OID,PG_BLOB,PG_BYTEA]
  getTablesSql = """select tablename from pg_tables where schemaname='public'"""
  getTablesAndViewsSql = """SELECT c.relname as "Name"
  FROM pg_catalog.pg_class c
  LEFT JOIN pg_catalog.pg_user u ON u.usesysid = c.relowner
  LEFT JOIN pg_catalog.pg_namespace n ON n.oid = c.relnamespace
  WHERE c.relkind IN ('r','v','S','')
  AND n.nspname NOT IN ('pg_catalog', 'pg_toast')
  AND pg_catalog.pg_table_is_visible(c.oid)
                              
  """
  getDbSql = """ select datname from pg_database where datallowconn """
  fileWildcard=None
  placeHolder='%s'
  binaryTypeName="bytea"
  binaryHolder = PgBytea
  RDTestDatabase="::RDTests"
elif hasattr(RDConfig,"useSqlLite") and RDConfig.useSqlLite:
  useGvib=0
  try:
    import sqlite3 as sqlite
    #from sqlite3 import *
  except ImportError:
    from pysqlite2 import dbapi2 as sqlite
    #from pysqlite2 import *
  sqlTextTypes = []
  sqlIntTypes = []
  sqlFloatTypes = []
  sqlBinTypes = []
  getTablesSql = """select name from SQLite_Master where type='table'"""
  getTablesAndViewsSql = """select name from SQLite_Master where type in ('table','view')"""
  getDbSql = None
  dbFileWildcard='*.sqlt'
  placeHolder='?'
  binaryTypeName="blob"
  binaryHolder = buffer

  connect = lambda x,*args:sqlite.connect(x)
else:
  useGvib=1
  from gvib import *
  sqlTextTypes = [gvibTypes.SQL_TEXT,gvibTypes.SQL_VARYING]
  sqlIntTypes = [gvibTypes.SQL_SHORT,gvibTypes.SQL_LONG,
                 gvibTypes.SQL_INT64]
  sqlFloatTypes = [gvibTypes.SQL_FLOAT,gvibTypes.SQL_DOUBLE,
                   gvibTypes.SQL_D_FLOAT]
  sqlBinTypes = [gvibTypes.SQL_BLOB,gvibTypes.SQL_ARRAY]
  getTablesSql = """SELECT RDB$RELATION_NAME
    FROM RDB$RELATIONS
    WHERE ((RDB$SYSTEM_FLAG = 0) OR 
    (RDB$SYSTEM_FLAG IS NULL)) AND 
    (RDB$VIEW_SOURCE IS NULL)
    ORDER BY RDB$RELATION_NAME"""
  getTablesAndViewsSql="""SELECT RDB$RELATION_NAME
      FROM RDB$RELATIONS
      WHERE ((RDB$SYSTEM_FLAG = 0) OR 
      (RDB$SYSTEM_FLAG IS NULL))
      ORDER BY RDB$RELATION_NAME"""
  getDbSql = None
  dbFileWildcard='*.gdb'
  placeHolder='?'
  binaryTypeName="blob"
  binaryHolder = lambda x:x
