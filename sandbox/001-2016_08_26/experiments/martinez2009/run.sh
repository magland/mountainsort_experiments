#!/bin/bash

num=${1:-1}
raw=../../raw/martinez2009/$num

mkdir output-$num
mountainprocess run-script --_nodaemon alg_scda_009.js $raw/params.json --raw=$raw/raw.mda.prv --outpath=output-$num --mask_threshold=0 | tee sort.log
