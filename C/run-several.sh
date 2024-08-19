#!/bin/bash

# File paths
DAT_FILE="../datafiles/reinject_regions.dat"
CONFIG_FILE="config/config2.ini"
CSCRIPT="reinjection_region_lorenz_gsl_omp"

# Read the .dat file line by line
i=1
while IFS= read -r line; do
	# Extract xmin and xmax from line
	xmin=$(echo $line | cut -d':' -f1)
	xmax=$(echo $line | cut -d':' -f2)
	rtarget=$(echo $line | cut -d':' -f3)
	
	# Update the config.ini file
	sed -i "s/^xmin = .*/xmin = $xmin/" $CONFIG_FILE
	sed -i "s/^xmax = .*/xmax = $xmax/" $CONFIG_FILE
	sed -i "s/^rtarget = .*/rtarget = $rtarget/" $CONFIG_FILE
	
	# Define the output datafile
	DATAFILE="../datafiles/reinjection_region_lorenz_gsl_omp_$i.dat"
	
	# Run the C script
	./run.sh $CSCRIPT $DATAFILE $CONFIG_FILE
	
	# Increment the counter
	i=$((i + 1))
done < "$DAT_FILE"
