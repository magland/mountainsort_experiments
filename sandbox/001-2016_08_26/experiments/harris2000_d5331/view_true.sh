#!/bin/bash

raw=../../raw/harris2000_d5331
mountainview --firings=$raw/firings_true.mda.prv --raw=output/pre0.mda.prv --filt=output/pre1b.mda.prv --pre=output/pre2.mda.prv --samplerate=10000 --geom=geom.csv
