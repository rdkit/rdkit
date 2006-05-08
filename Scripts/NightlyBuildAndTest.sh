#!/bin/bash

DATESTRING=`date +%d_%m_%y`
DIRNAME="Build_"$DATESTRING
export BASE=/tmp/$DIRNAME
#OLD_DIRS=/tmp/BUILD_*`date +"_%m_%y"`

export RDBASE=$BASE/RD
export RDOPTFLAGS="-O3"
export RDF77LIB=""
export PYTHON_ROOT="/usr"
export PYTHON_VERSION="2.4"
export BOOSTBASE="boost-1_33"
export PYTHONPATH="$RDBASE/Python"
export PATH="$RDBASE/bin:$PATH"
export LD_LIBRARY_PATH="$RDBASE/bin:/usr/local/lib:$LD_LIBRARY_PATH"
LOGFILE=$BASE"/build_log_"$DATESTRING

# ------------------------- -------------------------
#
#               Setup
#
# ------------------------- -------------------------
rm -rf /tmp/Build_*`date +"_%m_%y"`

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
SVNROOT=https://svn.sourceforge.net/svnroot/rdkit/trunk
svn checkout  -N $SVNROOT RD &> /dev/null
cd RD
svn checkout  -N $SVNROOT/Data Data &> /dev/null
svn checkout  -N $SVNROOT/Data/NCI Data/NCI &> /dev/null
for foo in Code bin Python; do
  svn checkout  $SVNROOT/$foo $foo &> /dev/null
done
svn checkout  -N $SVNROOT/External External &> /dev/null
for foo in Lapack++ libsvm svdlibc svdpackc vflib-2.0 cmim-1.0 HappyDoc-r1_3; do
  svn checkout  $SVNROOT/External/$foo External/$foo &> /dev/null
done

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
cd $RDBASE/Code/RDGeneral
make install >> $LOGFILE 2>&1
cd $RDBASE/Code/RDBoost 
make >> $LOGFILE 2>&1
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
rm -rf /home/www/html/RDCodeDocs
mv docs /home/www/html/RDCodeDocs
rm -rf /home/www/html/RDPythonDocs
mv $RDBASE/Docs/Code /home/www/html/RDPythonDocs
cp $RDBASE/External/HappyDoc-r1_3/RD.css /home/www/html/RDPythonDocs

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
~glandrum/RD/trunk/Scripts/MailResults.py $LOGFILE.gz $LOGFILE.summary
~glandrum/RD/trunk/Scripts/AddIssue.py $LOGFILE.summary "user=NightlyBuild" "tracker_home=/home/roundup/trackers/RDTrack" "title=Test Failures: $DATESTRING" "priority=bug"









