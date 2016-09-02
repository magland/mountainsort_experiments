#!/bin/bash

print_usage()
{
	echo "Usage: [].sh [algnames|all] [dsnames|all] [firings.mda|firings.curated.mda|optional]"
	echo "Example: [].sh ms example1,example2"
}

alg=$1
ds=$2
firings_name=${3:-firings.mda}

if [ -z $alg ] || [ -z $ds ] || [ -z $firings_name ] ; then
	print_usage
	exit -1
fi;

../boxsort/tabulate_results.node.js $alg $ds --firings_name=$firings_name --outpath=$PWD/output --alglist=$PWD/alglist.txt --dslist=$PWD/dslist.txt  "$@"