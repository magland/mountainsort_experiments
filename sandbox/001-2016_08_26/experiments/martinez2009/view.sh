#!/bin/bash

num=${1:-1}

mountainview --firings=output-$num/firings.mda --raw=output-$num/pre0.mda.prv --filt=output-$num/pre1b.mda.prv --pre=output-$num/pre2.mda.prv --samplerate=24000
