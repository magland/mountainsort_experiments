#!/bin/bash
algname=$1
dsname=$2
nodejs ../boxsort/boxsort2.node.js $algname $dsname --outpath=$PWD/output --alglist=$PWD/alglist.txt --dslist=$PWD/dslist.txt "$@"