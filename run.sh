#!/bin/bash

SEED=42

for CELL_LINE in 'GM12878' 'HUVEC' 'HeLa-S3' 'K562'; do
    for CROSS_CELL_LINE in 'GM12878' 'HUVEC' 'HeLa-S3' 'K562'; do
        printf "\n\n\nPREPARING DATA FOR '$CELL_LINE' AND '$CROSS_CELL_LINE' \n\n\n"
        python prepare_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --balanced

        printf "\n\n\nSPLITTING DATA FOR '$CELL_LINE' AND '$CROSS_CELL_LINE' \n\n\n"
        python split_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --seed=$SEED

        printf "\n\n\nTRAIN_TEST FOR '$CELL_LINE' AND '$CROSS_CELL_LINE' \n\n\n"
        python train_test.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --seed=$SEED
    done
done


printf "\n\n\nPREPARING DATA FOR combined \n\n\n"
python prepare_data.py --cell_line='combined' --balanced

printf "\n\n\nSPLITTING DATA FOR combined \n\n\n"
python split_data.py --cell_line='combined' --seed=$SEED

printf "\n\n\nTRAIN_TEST FOR combined \n\n\n"
python train_test.py --cell_line='combined' --seed=$SEED
