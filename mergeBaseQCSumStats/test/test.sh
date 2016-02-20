#!/bin/bash

status=0;
../../bin/mergeBaseQCSumStats --out results/mergeBaseQCSum.stats testFiles/test1.stats testFiles/test2.stats testFiles/test3.stats testFiles/test4.stats testFiles/test5.stats 2> results/mergeBaseQCSum.log
let "status |= $?"
diff results/mergeBaseQCSum.stats expected/mergeBaseQCSum.stats
let "status |= $?"
diff results/mergeBaseQCSum.log expected/mergeBaseQCSum.log
let "status |= $?"

../../bin/mergeBaseQCSumStats --out results/mergeBaseQCSumFromFile.stats --chrList testFiles/chrList.txt testFiles/test1.stats testFiles/test2.stats testFiles/test3.stats testFiles/test4.stats testFiles/test5.stats 2> results/mergeBaseQCSumFromFile.log
let "status |= $?"
diff results/mergeBaseQCSumFromFile.stats expected/mergeBaseQCSum.stats
let "status |= $?"
diff results/mergeBaseQCSumFromFile.log expected/mergeBaseQCSumFromFile.log
let "status |= $?"

../../bin/mergeBaseQCSumStats --out results/mergeBaseQCSumShort.stats testFiles/test1short.stats testFiles/test2short.stats testFiles/test3short.stats testFiles/test4short.stats testFiles/test5short.stats 2> results/mergeBaseQCSumShort.log
let "status |= $?"
diff results/mergeBaseQCSumShort.stats expected/mergeBaseQCSumShort.stats
let "status |= $?"
diff results/mergeBaseQCSumShort.log expected/mergeBaseQCSumShort.log
let "status |= $?"

../../bin/mergeBaseQCSumStats --out results/mergeBaseQCSumFromFileShort.stats --chrList testFiles/chrList.txt testFiles/test1short.stats testFiles/test2short.stats testFiles/test3short.stats testFiles/test4short.stats testFiles/test5short.stats 2> results/mergeBaseQCSumFromFileShort.log
let "status |= $?"
diff results/mergeBaseQCSumFromFileShort.stats expected/mergeBaseQCSumShort.stats
let "status |= $?"
diff results/mergeBaseQCSumFromFileShort.log expected/mergeBaseQCSumFromFileShort.log
let "status |= $?"

if [ $status != 0 ]
then
  echo failed mergeBaseQCSum test.
  exit 1
fi

