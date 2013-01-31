#!/bin/sh
#$ -o /srv/GT/analysis/masaomi/workflow_manager/out.log
#$ -e /srv/GT/analysis/masaomi/workflow_manager/err.log
#PARAMETER SLEEP_TIME
#PARAMETER ECHO_STRING

SLEEP_TIME=10
ECHO_STRING=hogehoge

echo "START"
echo $ECHO_STRING
ls
pwd
sleep $SLEEP_TIME
pwd
sleep $SLEEP_TIME
pwd
sleep $SLEEP_TIME
pwd
ls
echo 'END'
