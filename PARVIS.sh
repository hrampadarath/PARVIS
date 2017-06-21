#!/usr/bin/env bash


NOW=$(date +"%d-%b-%y:%T")

simulator="PARVIS.py"
echo "executing $simulator"
dir="./logs/"
LOGFILE="$dir$simulator-$NOW.log"
echo "Saving logs to $LOGFILE"

python $simulator -f simulation.inputs > $LOGFILE

echo "$simulator completed"
