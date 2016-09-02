#!/bin/bash
algname=$1
dsname=$2
firings_name=$3

nodejs ../boxsort/tabulate_results.node.js $algname $dsname --firings_name=$firings_name --outpath=$PWD/output --alglist=$PWD/alglist.txt --dslist=$PWD/dslist.txt  "$@"