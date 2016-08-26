#!/bin/bash

raw=../../raw/FL_ms11d45_18chan

mkdir output
mountainprocess run-script --_nodaemon alg_scda_009.js $raw/params.json --raw=$raw/raw.mda.prv --geom=$raw/geom.csv --outpath=output | tee sort.log
