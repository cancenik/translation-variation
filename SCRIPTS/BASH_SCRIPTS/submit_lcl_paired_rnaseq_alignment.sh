#!/bin/bash

LCLRNA=/srv/gs1/projects/snyder/ccenik/LCL_RNASEQ
cd $LCLRNA
for dir in $(find . -maxdepth 2 -mindepth 1 -name 'GM*' -type d)
do  
if [ ! -e $LCLRIBO/$dir/LOG* ]
then
cd $LCLRNA/$dir/
FILE1=$(ls *_1*)
FILE2=$(ls *_2*)
pwd
echo "qsub -v first=$FILE1,second=$FILE2 ~/SCRIPTS/BASH_SCRIPTS/alignment_strategy_paired.sh"
qsub -v first=$FILE1,second=$FILE2 ~/SCRIPTS/BASH_SCRIPTS/alignment_strategy_paired.sh
fi
done