# $Id$
#
# Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import RDConfig

if hasattr(RDConfig, "usePgSQL") and RDConfig.usePgSQL:
  from pyPgSQL import PgSQL

  # as of this writing (March 2004), this results in a speedup in
  # getting results back from the wrapper:
  PgSQL.fetchReturnsList = 1

  from pyPgSQL.PgSQL import *
  sqlTextTypes = [PG_CHAR, PG_BPCHAR, PG_TEXT, PG_VARCHAR, PG_NAME]
  sqlIntTypes = [PG_INT8, PG_INT2, PG_INT4]
  sqlFloatTypes = [PG_FLOAT4, PG_FLOAT8]
  sqlBinTypes = [PG_OID, PG_BLOB, PG_BYTEA]
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
  fileWildcard = None
  placeHolder = '%s'
  binaryTypeName = "bytea"
  binaryHolder = PgBytea
  RDTestDatabase = "::RDTests"
elif hasattr(RDConfig, "useSqlLite") and RDConfig.useSqlLite:
  try:
    import sqlite3 as sqlite
  except ImportError:
    from pysqlite2 import dbapi2 as sqlite
  sqlTextTypes = []
  sqlIntTypes = []
  sqlFloatTypes = []
  sqlBinTypes = []
  getTablesSql = """select name from SQLite_Master where type='table'"""
  getTablesAndViewsSql = """select name from SQLite_Master where type in ('table','view')"""
  getDbSql = None
  dbFileWildcard = '*.sqlt'
  fileWildcard = dbFileWildcard
  placeHolder = '?'
  binaryTypeName = "blob"
  binaryHolder = memoryview

  def connect(x, *args):
    return sqlite.connect(x)
else:
  raise ImportError("Neither sqlite nor PgSQL support found.")
