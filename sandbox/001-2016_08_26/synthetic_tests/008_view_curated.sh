#!/bin/bash

print_usage()
{
	echo "Usage: [].sh ms example1"
	echo "Usage: [].sh ms example2"
	echo "etc"
	echo "Usage: [].sh ms example1 firings.annotated.mda"
}

alg=$1
ds=$2
firings_name=${3:-firings.mda}

if [ -z $alg ]; then
	print_usage
	exit -1
fi;

outpath=output/$alg-$ds
mountainview --firings=$outpath/$firings_name --raw=$outpath/pre0.mda.prv --filt=$outpath/pre1b.mda.prv --pre=$outpath/pre2.mda.prv --samplerate=30000


