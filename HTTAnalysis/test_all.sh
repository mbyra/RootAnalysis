#!/bin/bash

#This script runs test equal to "./test htt.ini" N times, each time for 1, 2, 3...T threads.
#Creates if necessary directory ./tests, and subdirectories for each number of threads with subdirectories for each test.
#Saves output of each test and resulting png_jpg images to proper directories.
# USAGE: bash test_all.sh N T 

N_TESTS=1 #number of tests for each number of threads passed by user
T_THREADS=10 #maximal number of threads to use during each test

usage() {
	echo "Usage: bash ./test_all N T"
	echo "(Run N times ''./test htt.ini'', each time for 1,2...T threads, save output and png_jpg to ./tests/T/N"
}

parse_argument() {
	if [[ $# -ne 2 ]]; then
		echo "Illegal number of parameters"
		usage
		exit 1
	else 
		N_TESTS=$1
		T_THREADS=$2
	fi;
}

create_directories() {
	mkdir -p tests #create tests directory if not exist
	for (( i=1; $i <= $T_THREADS; i++ )); do
		mkdir -p tests/$i #create tests/thread_no directory if not exist
		cp htt.ini tests/$i/ #copy main htt.ini to this directory
		echo "%s/threads = 1/threads = $i/g
		w
		q
		" | ex tests/$i/htt.ini #change "threads = 1" to "threads = i"
		for (( j=1; $j <= $N_TESTS; j++ )); do
			mkdir -p tests/$i/$j #create subdirectories for each of N_TESTS tests
		done
	done
}

run_tests() {
	for (( i=1; $i <= $T_THREADS; i++ )); do
			htt_ini_path="tests/$i/htt.ini" #path to htt.ini containing "threads = i"
		for (( j=1; $j <= $N_TESTS; j++ )); do
			echo -n "./test $htt_ini_path"
			echo ""
			./test tests/$i/htt.ini 2>&1 | tee tests/$i/$j/output.txt #run test with output redirected also to file in test directory
			cp -R fig_jpg tests/$i/$j #copy directory containing images to proper directory
		done
	done
}

parse_argument $@ #check if user entered number of tests to run, otherwise print usage info
create_directories #create necessary directories to run and store test results
run_tests