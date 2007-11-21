#!/bin/bash

DATESTRING=`date +%d_%m_%y`
DIRNAME="RDKitBuild_"$DATESTRING
TMPDIR=/scratch
export BASE=$TMPDIR/$DIRNAME

#MAILIT=/home/glandrum/RDKit/Scripts/MailResults.py
export RDBASE=$BASE/RDKit
export RDOPTFLAGS="-O3"
export RDF77LIB="gfortran"
export PYTHON_ROOT="/usr"
export PYTHON_VERSION="2.5"
export GCCVERSION="41"
export BOOSTBASE="boost-1_34_1"
export BOOSTHOME="/usr/local"
export PYTHONPATH="$RDBASE/Python"
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
svn checkout  -N $SVNROOT RDKit &> /dev/null
cd $RDBASE
svn checkout  -N $SVNROOT/Data Data &> /dev/null
svn checkout  -N $SVNROOT/Data/NCI Data/NCI &> /dev/null
for foo in Code bin Python Projects; do
  svn checkout  $SVNROOT/$foo $foo &> /dev/null
done
svn checkout  -N $SVNROOT/External External &> /dev/null
for foo in Lapack++ svdlibc svdpackc vflib-2.0 cmim-1.0 HappyDoc-r1_3; do
  svn checkout  $SVNROOT/External/$foo External/$foo &> /dev/null
done
svn checkout  -N $SVNROOT/Scripts Scripts &> /dev/null



if test -f $RDBASE/Data/RDTests.sqlt; then 
  rm $RDBASE/Data/RDTests.sqlt
fi
$DBLOADER $RDBASE/Data/RDTests.sqlt < $RDBASE/Python/Dbase/testData/RDTests.sqlite
if test -f $RDBASE/Data/RDData.sqlt; then 
  rm $RDBASE/Data/RDData.sqlt
fi
$DBLOADER $RDBASE/Data/RDData.sqlt < $RDBASE/Python/Dbase/testData/RDData.sqlite 



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
cd $RDBASE/External
make >> $LOGFILE 2>&1
cd $RDBASE/Code
make all >> $LOGFILE 2>&1
make regrs >> $LOGFILE 2>&1


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
python $RDBASE/Python/TestRunner.py test_list.py >> $LOGFILE 2>&1
cd GraphMol
python $RDBASE/Python/TestRunner.py test_list.py >> $LOGFILE 2>&1
cd $RDBASE/Python
find . -name 'test_list.py' -exec python $RDBASE/Python/TestRunner.py \{\} >> $LOGFILE 2>&1 \; 
cd $RDBASE/Projects
python $RDBASE/Python/TestRunner.py test_list.py >> $LOGFILE 2>&1




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
make docs
cd $RDBASE/Python
epydoc --config  epydoc.config


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









