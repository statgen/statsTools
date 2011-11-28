#!/bin/bash

status=0;
../../bin/mergeBaseQCSumStats --out results/mergeBaseQCSum.stats testFiles/test1.stats testFiles/test2.stats testFiles/test3.stats testFiles/test4.stats testFiles/test5.stats 2> results/mergeBaseQCSum.log
let "status |= $?"
diff results/mergeBaseQCSum.stats expected/mergeBaseQCSum.stats
let "status |= $?"
diff results/mergeBaseQCSum.log expected/mergeBaseQCSum.log
let "status |= $?"

if [ $status != 0 ]
then
  echo failed mergeBaseQCSum test.
  exit 1
fi

