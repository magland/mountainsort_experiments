#!/bin/bash

alg=${1:-ms80}
ds=${2:-h1}

outpath=output/$alg-$ds

firings_true=$outpath/firings_true.mda.prv
if [ ! -f $firings_true ]; then
	firings_true=$outpath/firings_true.mda
fi

mountaincompare --firings1=$firings_true --firings2=$outpath/firings.mda --raw=$outpath/pre0.mda.prv --filt=$outpath/pre1b.mda.prv --pre=$outpath/pre2.mda.prv --samplerate=30000
