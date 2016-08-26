#!/bin/bash

raw=../../raw/FL_20160426_r1_nt16

mkdir output
mountainprocess run-script --_nodaemon alg_scda_009.js $raw/params.json --raw=$raw/raw.mda.prv --outpath=output | tee sort.log
