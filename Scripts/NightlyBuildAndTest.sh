#!/bin/bash

DATESTRING=`date +%d_%m_%y`
DIRNAME="RDKitBuild_"$DATESTRING
TMPDIR=/scratch
export BASE=$TMPDIR/$DIRNAME

#MAILIT=/home/glandrum/RDKit/Scripts/MailResults.py
export RDBASE=$BASE/RDKit
export BOOSTBASE="boost_1_37_0"
export BOOSTHOME="/usr/local/src/$BOOSTBASE"
export BOOST_BUILD_PATH=$BOOSTHOME
export BJAM="/usr/local/bin/bjam"
export PYTHONPATH="$RDBASE"
export PATH="$RDBASE/bin:$PATH"
export LD_LIBRARY_PATH="$RDBASE/bin:/usr/local/lib"
LOGFILE=$BASE"/build_log_"$DATESTRING

export RD_USESQLLITE="1"
DBLOADER=sqlite3

# ------------------------- -------------------------
#
#               Setup
#
# ------------------------- -------------------------
rm -rf $TMPDIR/RDKitBuild_*`date +"_%m_%y"`

if ! rm -rf $BASE; then
    echo "CANNOT CLEANUP";
    exit -1
fi
if ! mkdir $BASE; then
    echo "CANNOT MAKE DIRECTORY";
    exit -1
fi
if ! cd $BASE; then
    echo "CANNOT CD";
    exit -1
fi

echo "Running in directory" $DIRNAME > $LOGFILE

# ------------------------- -------------------------
#
#               Pull
#
# ------------------------- -------------------------
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
echo "****  PULL  ****" >> $LOGFILE 2>&1
echo "****  PULL  ****"
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
SVNROOT=https://rdkit.svn.sourceforge.net/svnroot/rdkit/trunk
svn checkout  $SVNROOT RDKit &> /dev/null

if test -f $RDBASE/Data/RDTests.sqlt; then 
  rm $RDBASE/Data/RDTests.sqlt
fi
$DBLOADER $RDBASE/Data/RDTests.sqlt < $RDBASE/rdkit/Dbase/test_data/RDTests.sqlite
if test -f $RDBASE/Data/RDData.sqlt; then 
  rm $RDBASE/Data/RDData.sqlt
fi
$DBLOADER $RDBASE/Data/RDData.sqlt < $RDBASE/rdkit/Dbase/test_data/RDData.sqlite
# ------------------------- -------------------------
#
#               Build
#
# ------------------------- -------------------------
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
echo "****  BUILD  ****" >> $LOGFILE 2>&1
echo "****  BUILD  ****"
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
cd $RDBASE/Code
$BJAM -j2 >>$LOGFILE 2>&1

# ------------------------- -------------------------
#
#               Test
#
# ------------------------- -------------------------
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
echo "****  TEST  ****" >> $LOGFILE 2>&1
echo "****  TEST  ****"
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
cd $RDBASE/Code
python $RDBASE/rdkit/TestRunner.py test_list.py >> $LOGFILE 2>&1
cd $RDBASE/rdkit
python $RDBASE/rdkit/TestRunner.py test_list.py >> $LOGFILE 2>&1
cd $RDBASE/Projects
python $RDBASE/rdkit/TestRunner.py test_list.py >> $LOGFILE 2>&1


# ------------------------- -------------------------
#
#               Generate documentation
#
# ------------------------- -------------------------
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
echo "****  DOCS  ****" >> $LOGFILE 2>&1
echo "****  DOCS  ****"
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
cd $RDBASE/Code
/usr/bin/doxygen doxygen.config > /dev/null 2>&1
cd $RDBASE/rdkit
epydoc --config  epydoc.config > /dev/null 2>&1

# ------------------------- -------------------------
#
#               Report
#
# ------------------------- -------------------------
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
echo "****  TEST FAILURES  ****" >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
echo >> $LOGFILE 2>&1
grep -n "Failed [0-9]" $LOGFILE >> $LOGFILE.summary
grep -n "Failed [0-9]" $LOGFILE >> $LOGFILE
gzip -9 $LOGFILE
#$MAILIT $LOGFILE.gz $LOGFILE.summary
#$RDBASE/Scripts/AddIssue.py $LOGFILE.summary "user=NightlyBuild" "tracker_home=/home/roundup/trackers/RDTrack" "title=Test Failures: $DATESTRING" "priority=bug"









