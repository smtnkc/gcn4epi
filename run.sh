#!/bin/bash

SEED=42 

for CELL_LINE in 'GM12878' 'HUVEC' 'HeLa-S3' 'K562' 'combined'; do
    printf "\n\n\n\n\n PREPARING SAME-CELL-LINE DATA FOR '$CELL_LINE' \n\n\n\n\n"
    python prepare_gcn_data.py --cell_line="$CELL_LINE" --from_scratch --balanced --seed=$SEED
    printf "\n\n\n\n\n RUNNING SAME-CELL-LINE TESTS FOR '$CELL_LINE' \n\n\n\n\n"
    python train.py --cell_line="$CELL_LINE" --seed=$SEED
done

for CELL_LINE in 'GM12878' 'HUVEC' 'HeLa-S3' 'K562'; do
    for CROSS_CELL_LINE in 'GM12878' 'HUVEC' 'HeLa-S3' 'K562'; do
        printf "\n\n\n\n\n PREPARING CROSS-CELL-LINE DATA FOR '$CELL_LINE' AND '$CROSS_CELL_LINE' \n\n\n\n\n"
        python prepare_gcn_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --from_scratch --balanced --seed=$SEED
        printf "\n\n\n\n\n RUNNING CROSS-CELL-LINE TESTS FOR '$CELL_LINE' AND '$CROSS_CELL_LINE' \n\n\n\n\n"
        python train.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --seed=$SEED
    done
done
