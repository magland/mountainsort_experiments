#!/bin/bash

raw=../../raw/harris2000_d5331

mkdir output
mountainprocess run-script --_nodaemon alg_scda_009.js $raw/params.json --raw=$raw/raw.mda.prv --outpath=output | tee sort.log
