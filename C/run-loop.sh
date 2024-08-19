#!/bin/bash

# File paths
DAT_FILE1="../datafiles/m_function_regions_numerical_lorenz_"
DAT_FILE2="../datafiles/rpd_numerical_lorenz_"
READ_FILE="../datafiles/reinjection_region_lorenz_"
CONFIG_FILE="config/config2.ini"
CSCRIPT1="m_function_lorenz"
CSCRIPT2="rpd_regions_numerical_lorenz"

for i in {0..9}
do
	# Define the output datafile
	DATAFILE1="${DAT_FILE1}${i}.dat"
	DATAFILE2="${DAT_FILE2}${i}.dat"
	READFILE="${READ_FILE}${i}.dat"

	# Run the C script
	./run.sh $CSCRIPT1 $DATAFILE1 $READFILE
	echo "=============================================================================="
	echo "=============================================================================="
	./run.sh $CSCRIPT2 $DATAFILE2 $READFILE
	echo "=============================================================================="
	echo "=============================================================================="

done
