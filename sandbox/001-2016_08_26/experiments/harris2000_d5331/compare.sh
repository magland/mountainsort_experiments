#!/bin/bash

raw=../../raw/harris2000_d5331
mountaincompare --firings1=$raw/firings_true.mda.prv --firings2=output/firings.mda --raw=output/pre0.mda.prv --filt=output/pre1b.mda.prv --pre=output/pre2.mda.prv --samplerate=30000 --geom=geom.csv
