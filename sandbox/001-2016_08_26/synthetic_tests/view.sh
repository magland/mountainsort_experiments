#!/bin/bash

alg=${1:-ms}
ds=${2:-example1}

outpath=output/$alg-$ds

echo mountainview --firings=$outpath/firings.mda --raw=$outpath/pre0.mda.prv --filt=$outpath/pre1b.mda.prv --pre=$outpath/pre2.mda.prv --samplerate=30000
mountainview --firings=$outpath/firings.mda --raw=$outpath/pre0.mda.prv --filt=$outpath/pre1b.mda.prv --pre=$outpath/pre2.mda.prv --samplerate=30000
