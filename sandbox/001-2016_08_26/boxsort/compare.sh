#!/bin/bash

alg=${1:-ms80}
ds=${2:-h1}

outpath=output/$alg-$ds

mountaincompare --firings1=$outpath/firings_true.mda.prv --firings2=$outpath/firings.mda --raw=$outpath/pre0.mda.prv --filt=$outpath/pre1b.mda.prv --pre=$outpath/pre2.mda.prv --samplerate=30000
