#!/bin/bash

opts="--firing_rate_min=0.5 --firing_rate_max=3 --M=5 --duration=1800"

nodejs generate_synth_examples.node.js $opts --dsname=example1 --K=16 --noise_level=1 
nodejs generate_synth_examples.node.js $opts --dsname=example2 --K=16 --noise_level=2 
nodejs generate_synth_examples.node.js $opts --dsname=example3 --K=16 --noise_level=1 --amp_variation_min=0.2
nodejs generate_synth_examples.node.js $opts --dsname=example4 --K=16 --noise_level=2 --amp_variation_min=0.2
nodejs generate_synth_examples.node.js $opts --dsname=example5 --K=40 --noise_level=1 --amp_variation_min=0.2
