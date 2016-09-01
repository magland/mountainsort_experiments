#!/bin/bash

alg1=${1:-ms80}
ds1=${2:-example1}
alg2=${3:-ms80_noise}
ds2=${4:-example2}

outpath1=output/$alg1-$ds1
outpath2=output/$alg2-$ds2

mountaincompare --firings1=$outpath1/firings.mda --firings2=$outpath2/firings.mda --raw=$outpath1/pre0.mda.prv --filt=$outpath1/pre1b.mda.prv --pre=$outpath1/pre2.mda.prv --samplerate=30000
