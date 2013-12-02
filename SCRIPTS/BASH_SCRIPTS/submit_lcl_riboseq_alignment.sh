#!/bin/bash

LCLRIBO=/srv/gs1/projects/snyder/ccenik/LCL_RIBOSEQ
cd $LCLRIBO
for dir in $(find . -maxdepth 1 -mindepth 1 -name 'GM*' -type d -printf '%f\n')
do
if [ ! -e $LCLRIBO/$dir/LOG* ]
then
cd $LCLRIBO/$dir/
FILE=$(ls)
echo $LCLRIBO/$dir/
qsub -v adapter=AGATCGGAAGAGCACACGTCT,file=$FILE ~/SCRIPTS/BASH_SCRIPTS/alignment_strategy_qsub.sh
fi
done