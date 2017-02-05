#!/bin/bash

# This script runs test equal to "./test htt.ini" N times, each time for 1, 2, 3...T threads.
# Creates (if necessary) directory ./tests, and subdirectories for each number of threads with subdirectories for each test.
# Saves output of each test and resulting png_jpg images to proper directories.
# Measures time of each test to compare time on different schedules. To change directory in tests/ which will keep results of tests
# of current schedule from "dynamic", add third argument - name of directory e.g. "bash test_all.sh 2 10 guided"
# USAGE: bash test_all.sh N T    OR    bash test_all.sh N T SCHEDULE_DIRECTORY_NAME 

N_TESTS=1 #number of tests for each number of threads passed by user
T_THREADS=20 #maximal number of threads to use during each test
SCHEDULE="dynamic"

usage() {
	echo "Usage: bash ./test_all N T"
	echo "(Run N times ''./test htt.ini'', each time for 1,2...T threads, save output and png_jpg to ./tests/T/N"
}

parse_arguments() {
	if [[ $# -ne 2 ]] && [[ $# -ne 3 ]]; then
		echo "Illegal number of parameters"
		usage
		exit 1
	else 
		N_TESTS=$1
		T_THREADS=$2

	if [[ $# -eq 3 ]]; then
		SCHEDULE=$3
	fi;
}

create_directories() {
	mkdir -p tests #create tests directory if not exist
	for (( i=1; $i <= $T_THREADS; i++ )); do
		path="tests/$SCHEDULE/$i/" #path to htt.ini containing "threads = i"
		mkdir -p $path #create tests/thread_no directory if not exist
		cp htt.ini $path/ #copy main htt.ini to this directory
		echo "%s/threads = 1/threads = $i/g
		w
		q
		" | ex $path/htt.ini #change "threads = 1" to "threads = i"
		for (( j=1; $j <= $N_TESTS; j++ )); do
			mkdir -p $path/test_$j #create subdirectories for each of N_TESTS tests
		done
	done
}

run_tests() {
	for (( i=1; $i <= $T_THREADS; i++ )); do
			path="tests/$SCHEDULE/$i/" #path to main directory of given schedule and currently iterated amount of threads
		for (( j=1; $j <= $N_TESTS; j++ )); do
			echo -n "./test $path/htt.ini   schedule = $SCHEDULE, test no = $j"
			echo ""
			#./test tests/$i/htt.ini 2>&1 | tee tests/$i/$j/output.txt #run test with output redirected also to file in test directory
			{ time ./test $path/htt.ini 2>&1 $path/test_$j/output.txt ; } 2>&1 $path/test_$j/time.txt #run test with output redirected also to file in test directory
			cp -R fig_png $path/test_$j/ #copy directory containing images to proper directory
		done
	done
}

grep_outputs() {
	for (( i=1; $i <= $T_THREADS; i++ )); do
		path="tests/$SCHEDULE/$i/" #path to main directory of given schedule and currently iterated amount of threads
		for (( j=1; $j <= $N_TESTS; j++ )); do
			echo -n "schedule $SCHEDULE threads $i test $j. Analyzed data:"
			grep 'Data: ' $path/test_$j/output.txt
			echo ""
			grep 'real' $path/test_$j/time.txt
			echo ""
			grep 'user' $path/test_$j/time.txt
			echo ""
			grep 'sys' $path/test_$j/time.txt
			echo ""
			echo ""
		done
	done
}

parse_arguments $@ #check if user entered number of tests to run, otherwise print usage info
create_directories #create necessary directories to run and store test results
run_tests
grep_outputs #print number of processed data, and time of each evaluation