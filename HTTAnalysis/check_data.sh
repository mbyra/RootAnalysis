#!/bin/bash

# This script searches all outputs in directories of 1...T threads, each one with 1...N tests and prints number of analized data.
# USAGE: bash check_data.sh N T

N_TESTS=1 #number of tests for each number of threads passed by user
T_THREADS=10 #maximal number of threads to use during each test

usage() {
	echo "Usage: bash ./check_data N T"
	echo "(Searches all outputs in directories of 1...T threads, each one with 1...N tests and prints number of analized data)"
}

parse_arguments() {
	if [[ $# -ne 2 ]]; then
		echo "Illegal number of parameters"
		usage
		exit 1
	else 
		N_TESTS=$1
		T_THREADS=$2
	fi;
}

grep_outputs() {
	for (( i=1; $i <= $T_THREADS; i++ )); do
		for (( j=1; $j <= $N_TESTS; j++ )); do
			echo -n "Threads: $i. TestL $j. Analyzed data:"
			grep 'Data = ' tests/$i/$j/output.txt	
			echo ""
		done
	done
}

parse_arguments $@ #check if user entered number of tests to run, otherwise print usage info
grep_outputs