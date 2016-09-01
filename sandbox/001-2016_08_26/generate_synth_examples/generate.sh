#!/bin/bash

basepath=../../../../mountainsort_experiments

opts="--basepath=$basepath --firing_rate_min=0.5 --firing_rate_max=3 --M=5 --duration=1800"

nodejs generate_synth_examples.node.js $opts --dsname=example1 --K=16 --noise_level=1 
nodejs generate_synth_examples.node.js $opts --dsname=example2 --K=16 --noise_level=2 
nodejs generate_synth_examples.node.js $opts --dsname=example3 --K=16 --noise_level=2 --amp_variation_min=0.2