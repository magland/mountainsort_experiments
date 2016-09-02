#!/bin/bash

firings_name=${3:-firings.curated.mda}
./003_tabulate_results.sh $1 $2 $firings_name
