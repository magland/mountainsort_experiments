#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
nodejs $DIR/boxsort.node.js "$@"
