#!/bin/bash

set -euo pipefail

X_SEQ=$1
FOLDER_Y_SEQS=$2

mkdir -p hprc_pafs
for Y_SEQ in $(ls $FOLDER_Y_SEQS); do
    echo $X_SEQ
    echo $FOLDER_Y_SEQS/$Y_SEQ
    ./scripts/align_parent_to_child.sh $X_SEQ $FOLDER_Y_SEQS/$Y_SEQ 50000 > hprc_pafs/refto_${Y_SEQ}_50000.paf
    #echo "Done aligning $X_SEQ to $Y_SEQ"
done
