#!/bin/bash

scl enable devtoolset-3 bash
source /opt/soft/root-6.06.00/bin/thisroot.sh
cd /home/mb347056/RootAnalysis
git checkout IIpracownia_A
libtoolize
aclocal
automake --add-missing
autoreconf
./configure 
make