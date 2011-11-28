#!/bin/bash

status=0;
../../bin/subsetBaseQCStats --inStats testFiles/statsBaseQCSum.txt --regionList testFiles/regions.txt --outStats results/statsBaseQCSum.txt 2> results/statsBaseQCSum.log
let "status |= $?"
diff results/statsBaseQCSum.txt expected/statsBaseQCSum.txt
let "status |= $?"
diff results/statsBaseQCSum.log expected/statsBaseQCSum.log
let "status |= $?"

if [ $status != 0 ]
then
  echo failed subsetStats test.
  exit 1
fi

