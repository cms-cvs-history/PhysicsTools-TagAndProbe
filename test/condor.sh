#!/bin/bash

#
# variables from arguments string in jdl
#
# format:
#
# 1: condor cluster number
# 2: condor process number
# 3: CMSSW_DIR
# 4: RUN_DIR
# 5: PARAMETER_SET (full path, has to contain all needed files in PoolSource and filled following variables with keywords: maxEvents = CONDOR_MAXEVENTS, skipEvents = CONDOR_SKIPEVENTS, output fileName = CONDOR_OUTPUTFILENAME)
#

CONDOR_CLUSTER=$1
CONDOR_PROCESS=$2
CMSSW_DIR=$3
RUN_DIR=$4
PARAMETER_SET=$5

#
# header 
#

echo ""
echo "CMSSW on Condor"
echo ""

START_TIME=`/bin/date`
echo "started at $START_TIME"

echo ""
echo "parameter set:"
echo "CONDOR_CLUSTER: $CONDOR_CLUSTER"
echo "CONDOR_PROCESS: $CONDOR_PROCESS"
echo "CMSSW_DIR: $CMSSW_DIR"
echo "RUN_DIR: $RUN_DIR"
echo "PARAMETER_SET: $PARAMETER_SET"
#echo "NUM_EVENTS_PER_JOB: $NUM_EVENTS_PER_JOB"
#
# setup software environment at FNAL for the given CMSSW release
#
cd $CMSSW_DIR
source /uscmst1/prod/sw/cms/shrc uaf
export SCRAM_ARCH=slc4_ia32_gcc345
cd $CMSSW_DIR
eval `scramv1 runtime -sh`

#
# change to output directory
#
cd $RUN_DIR

#
# modify parameter-set
#

FINAL_PARAMETER_SET_NAME=`echo batch_${CONDOR_CLUSTER}_${CONDOR_PROCESS}`
FINAL_PARAMETER_SET=`echo $FINAL_PARAMETER_SET_NAME.py`
FINAL_LOG=`echo $FINAL_PARAMETER_SET_NAME.log`
FINAL_FILENAME=`echo $FINAL_PARAMETER_SET_NAME.root`
echo ""
echo "Writing final parameter-set: $FINAL_PARAMETER_SET to RUN_DIR: $RUN_DIR"
echo ""

# the following line leaves some room for ariphmetic operations
let "skip = $CONDOR_PROCESS"
cat $PARAMETER_SET | sed -e s/CONDOR_PROCESS/$skip/ | sed -e s/CONDOR_OUTPUTFILENAME/$FINAL_FILENAME/ > $FINAL_PARAMETER_SET

#
# run cmssw
#

echo "run: cmsRun $FINAL_PARAMETER_SET > $FINAL_LOG 2>&1"
cmsRun $FINAL_PARAMETER_SET >> $FINAL_LOG 2>&1
exitcode=$?

#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode

