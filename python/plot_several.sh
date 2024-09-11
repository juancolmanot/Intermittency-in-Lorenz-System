#!/bin/bash

BASE_FILENAME_1="reinjection_region_lorenz_"
BASE_FILENAME_2="m_function_regions_numerical_lorenz_"
BASE_FILENAME_3="rpd_regions_numerical_lorenz_"
CONFIG_FILE="plot_config5.txt"
PLOT_PREFIX="reinj_m_rpd_regions_numerical_lorenz_"

for i in {0..9}
do
	filename1="../datafiles/${BASE_FILENAME_1}${i}.dat"
	filename2="../datafiles/${BASE_FILENAME_2}${i}.dat"
	filename3="../datafiles/${BASE_FILENAME_3}${i}.dat"
	plotname="../plots/${PLOT_PREFIX}${i}.pdf"
	
	sed -i "s|^\([[:space:]]*file: \).*${BASE_FILENAME_1}[0-9]*\.dat|\1${filename1}|" "$CONFIG_FILE"
	sed -i "s|^\([[:space:]]*file: \).*${BASE_FILENAME_2}[0-9]*\.dat|\1${filename2}|" "$CONFIG_FILE"
	sed -i "s|^\([[:space:]]*file: \).*${BASE_FILENAME_3}[0-9]*\.dat|\1${filename3}|" "$CONFIG_FILE"
	sed -i "s|^save_file: .*|save_file: $plotname|" "$CONFIG_FILE"
	python3 plot_file.py "$CONFIG_FILE"
	echo "Plotted in ${plotname}"
done

