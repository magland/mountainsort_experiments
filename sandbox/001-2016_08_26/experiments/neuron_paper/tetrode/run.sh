#!/bin/bash

mkdir output
mountainprocess run-script --_nodaemon alg_scda_009.js --samplerate=30000 --raw=raw.mda.prv --outpath=output | tee sort.log
