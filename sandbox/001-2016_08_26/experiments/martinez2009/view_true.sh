#!/bin/bash

num=${1:-1}
raw=../../raw/martinez2009/$num

mountainview --firings=$raw/firings_true.mda.prv --raw=output-$num/pre0.mda.prv --filt=output-$num/pre1b.mda.prv --pre=output-$num/pre2.mda.prv --samplerate=10000 --geom=geom.csv
