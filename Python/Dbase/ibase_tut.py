# $Id$
#
#    Copyright (C) 2000-2006  greg Landrum
#
#   @@ All Rights Reserved  @@
#
"""
  A brief little demo/tutorial on use of the python <-> interbase adaptor
   gvib (http://www.zope.org/Members/RETierney/gvibDA), windows binary
     available from: http://clientes.netvisao.pt/luiforra/gvib/
"""

import gvib
import string

# default names for the database and table to be used
#  You'll need to change these to reflect your own system
defaultTestDBase = 'd:/glandrum/Python/Dbase.OLD/testdb.gdb'
defaultTestTable = 'blob_data2'

defaultUser='sysdba'
defaultPassword='masterkey'

#------------
#
# establish the connection and set up a cursor
#
#------------
cn = gvib.connect(defaultTestDBase,defaultUser,defaultPassword)
c = cn.cursor()

#------------
#
# clean up (remove) the table if it already exists in the database
#
#------------
try:
  # the SQL command to remove a table is 'drop table' (imagine that)
  c.execute('drop table %s'%(defaultTestTable))
except:
  print 'cannot drop'

#------------
#
# define a new table with two fields: a string and a BLOB
#  BLOBs can contain arbitrary data, but we're only putting
#    strings in for this demo.
#
#------------

# this is the definition of the columns
colDefs = 'name varchar(20), bdata blob'
# construct the SQL string to create the table
sqlStr = 'create table %s (%s)'%(defaultTestTable,colDefs)
# and execute it
c.execute(sqlStr)

#------------
#
# insert a new row into the table
#
#------------

# form a string of numbers to slap in the BLOB
numStr = string.join(map(lambda x:str(x),range(10)),' ')
dStr = numStr*6

#  This is an alternate set of code which inserts the contents of a file
#f = open('Gui.py','r')
#l = f.readlines()
#dStr = string.join(l,'')

# this is the SQL command to insert the data, the ? are placeholders
#   to tell gvib that data will follow
sqlStr = "insert into %s values (?, ?)"%(defaultTestTable)
# execute the SQL.  Pass in the data (for the ?s) as a *TUPLE*, even if
#  there's only one ? to fill.
c.execute(sqlStr,('hi',dStr))
# commit the change (otherwise the table is not updated)...
#  changes must be committed before the connection to the database is dropped or,
#  I think, before you query the database.
cn.commit()

#------------
#
# pull data out of the table
#
#------------

# select is used to grab data.  Here we just grab it all
c.execute('select * from %s'%(defaultTestTable))
# the data is now stored in the cursor.  Use fetchone() to pull out
#  one of the rows as a tuple
res= c.fetchone()
# and print the data.  This first form shows the internal representation
#  of a BLOB
print res
# and this shows that it can be printed out nicely
print res[1]

#------------
#
# Add a column to the table
#
#------------

# this adds a column called num which is an integer
c.execute('alter table %s add num integer'%(defaultTestTable))

# update the entry in the table
c.execute('update %s set num=12 where name=?'%(defaultTestTable),('hi',))
#  NOTE: there's no reason the 12 here couldn't be in the parameter
#  list at the end:
# c.execute('update %s set num=? where name=?'%(defaultTestTable),(12,'hi'))

# add a new row to the table
sqlStr = "insert into %s values (?, ?, ?)"%(defaultTestTable)
dStr = 'foo bar baz grn'
c.execute(sqlStr,('bye',dStr,14))
# commit the change
cn.commit()

# pull out all the data again
c.execute('select * from %s'%(defaultTestTable))

# print out the first entry to make sure it got updated
res = c.fetchone()
print res[0],res[1],res[2]
# print out the second entry to make sure it got updated
res = c.fetchone()
print res[0],res[1],res[2]

#------------
#
# Add some more rows to the table
#
#------------
sqlStr = "insert into %s values (?, ?, ?)"%(defaultTestTable)
# insert 4 rows using the executemany() command.  Here the SQL
#  command is run with each set of parameters.
c.executemany(sqlStr,[('foo',dStr,1),
                      ('bar',dStr,2),
                      ('baz',dStr,3),
                      ('grn',dStr,4)])
# don't forget to commit the change
cn.commit()

#------------
#
# Mildly fancier use of select
#
#------------

# pull out only the name and num fields from records where num < 10
c.execute('select name,num from %s where num<10'%(defaultTestTable))

# and print out all the results.  fetchall() gives back a list of them
print c.fetchall()

