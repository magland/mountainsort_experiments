#!/bin/bash

print_usage()
{
	echo "Usage: [].sh [algnames|all] [dsnames|all]"
	echo "Example: [].sh ms example1,example2"
}

alg=$1
ds=$2

if [ -z $alg ] || [ -z $ds ] ; then
	print_usage
	exit -1
fi;

algname=$1
dsname=$2
../boxsort/boxsort2.node.js $algname $dsname --outpath=$PWD/output --alglist=$PWD/alglist.txt --dslist=$PWD/dslist.txt "$@"