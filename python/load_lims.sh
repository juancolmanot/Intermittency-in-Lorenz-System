#!/bin/bash

DATA_FILE=$1
CONFIG_FILE="plot_config5.txt"

xmin=$(awk '{print $1}' $DATA_FILE | sort -n | head -1)
xmax=$(awk '{print $1}' $DATA_FILE | sort -n | tail -1)
ymin=$(awk '{print $2}' $DATA_FILE | sort -n | head -1)
ymax=$(awk '{print $2}' $DATA_FILE | sort -n | tail -1)

sed -i "s/^[[:space:]]*xlim: .*/\txlim: [$xmin, $xmax]/" $CONFIG_FILE
sed -i "s/^[[:space:]]*ylim: .*/\tylim: [$ymin, $ymax]/" $CONFIG_FILE

echo "Replaced xlim and ylim in $CONFIG_FILE with lims: [$xmin, $xmax], [$ymin, $ymax]"
