#!/bin/bash

SEED=$1

# for CELL_LINE in 'GM12878' 'K562' 'HUVEC' 'HeLa-S3' 'IMR90' 'NHEK'; do
#     python prepare_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CELL_LINE" --k_mer=5 --label_rate=0.2 --label=$LABEL --seed=$SEED --balanced
#     python split_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CELL_LINE" --k_mer=5 --label_rate=0.2 --seed=$SEED
#     python train_test.py --cell_line="$CELL_LINE" --cross_cell_line="$CELL_LINE" --k_mer=5 --label_rate=0.2 --label=$LABEL --seed=$SEED
# done

for CELL_LINE in 'GM12878' 'K562' 'HUVEC' 'HeLa-S3' 'NHEK' 'IMR90'; do
    for CROSS_CELL_LINE in 'GM12878' 'K562' 'HUVEC' 'HeLa-S3' 'NHEK' 'IMR90'; do
        printf "\n\n\nPREPARING DATA FOR $CELL_LINE AND $CROSS_CELL_LINE (5-MER)\n\n\n"
        python prepare_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=5 --label=1 --seed=$SEED --label_rate=0.2 --balanced

        printf "\n\n\nSPLITTING DATA FOR $CELL_LINE AND $CROSS_CELL_LINE (5-MER)\n\n\n"
        python split_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=5 --seed=$SEED --label_rate=0.2

        printf "\n\n\nTRAIN_TEST FOR $CELL_LINE AND $CROSS_CELL_LINE (5-MER)\n\n\n"
        python train_test.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=5 --label=1 --seed=$SEED --label_rate=0.2
    done
done

