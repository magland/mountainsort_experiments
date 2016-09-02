#!/bin/bash

print_usage()
{
	echo "Usage: [].sh ms example1"
	echo "Usage: [].sh ms example2"
	echo "etc"
}

if [ -z $1 ]; then
	print_usage
	exit -1
fi;

./003_compare_with_ground_truth.sh "$@" firings.annotated.mda


