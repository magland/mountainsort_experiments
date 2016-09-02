#!/bin/bash

print_usage()
{
	echo "Usage: [].sh [algname|all] [dsname|all] [firings.mda|firings.curated.mda|optional]"
	echo "Example: [].sh ms example1"
}

alg=$1
ds=$2
firings_name=${3:-firings.mda}

if [ -z $alg ] || [ -z $ds ] || [ -z $firings_name ] ; then
	print_usage
	exit -1
fi;

outpath=output/$alg-$ds
mountainview --firings=$outpath/$firings_name --raw=$outpath/pre0.mda.prv --filt=$outpath/pre1b.mda.prv --pre=$outpath/pre2.mda.prv --samplerate=30000


