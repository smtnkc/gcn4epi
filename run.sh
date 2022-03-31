#!/bin/bash

SEED=$1
PREPARE=$2
K_MER=$3

# for K_MER in 3 4 5 6; do
    # for CELL_LINE in 'combined'; do
    #     PREPARE=1
    #     for SEED in 42 91 123; do
    #         if [ "$PREPARE" = 1 ] ; then
    #             printf "\n\n\nPREPARING DATA FOR $CELL_LINE (K-MER=$K_MER SEED=$SEED)\n\n\n"
    #             python prepare_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CELL_LINE" --k_mer=$K_MER --balanced
    #         fi

    #         printf "\n\n\nSPLITTING DATA FOR $CELL_LINE (K-MER=$K_MER SEED=$SEED)\n\n\n"
    #         python split_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CELL_LINE" --k_mer=$K_MER --seed=$SEED

    #         printf "\n\n\nTRAIN_TEST FOR $CELL_LINE (K-MER=$K_MER SEED=$SEED)\n\n\n"
    #         python train_test.py --cell_line="$CELL_LINE" --cross_cell_line="$CELL_LINE" --k_mer=$K_MER --seed=$SEED
    #         PREPARE=0
    #     done
    # done
# done


# for CELL_LINE in 'GM12878' 'HUVEC' 'HeLa-S3' 'K562'; do
#     for CROSS_CELL_LINE in 'GM12878' 'HUVEC' 'HeLa-S3' 'K562'; do
#         if [ "$PREPARE" = 1 ] ; then
#             printf "\n\n\nPREPARING DATA FOR $CELL_LINE AND $CROSS_CELL_LINE ($K_MER-MER)\n\n\n"
#             python prepare_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=$K_MER --balanced
#         fi

#         printf "\n\n\nSPLITTING DATA FOR $CELL_LINE AND $CROSS_CELL_LINE ($K_MER-MER)\n\n\n"
#         python split_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=$K_MER --seed=$SEED

#         printf "\n\n\nTRAIN_TEST FOR $CELL_LINE AND $CROSS_CELL_LINE ($K_MER-MER)\n\n\n"
#         python train_test.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=$K_MER --seed=$SEED
#     done
# done

# if [ "$PREPARE" = 1 ] ; then
#     printf "\n\n\nPREPARING DATA FOR combined ($K_MER-MER)\n\n\n"
#     python prepare_data.py --cell_line='combined' --k_mer=$K_MER --balanced
# fi

# printf "\n\n\nSPLITTING DATA FOR combined ($K_MER-MER)\n\n\n"
# python split_data.py --cell_line='combined' --k_mer=$K_MER --seed=$SEED

# printf "\n\n\nTRAIN_TEST FOR combined ($K_MER-MER)\n\n\n"
# python train_test.py --cell_line='combined' --k_mer=$K_MER --seed=$SEED


for CELL_LINE in 'NHEK' 'IMR90'; do
    for CROSS_CELL_LINE in 'NHEK' 'IMR90'; do
        if [ "$PREPARE" = 1 ] ; then
            printf "\n\n\nPREPARING DATA FOR $CELL_LINE AND $CROSS_CELL_LINE ($K_MER-MER)\n\n\n"
            python prepare_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=$K_MER
        fi

        printf "\n\n\nSPLITTING DATA FOR $CELL_LINE AND $CROSS_CELL_LINE ($K_MER-MER)\n\n\n"
        python split_data.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=$K_MER --seed=$SEED

        printf "\n\n\nTRAIN_TEST FOR $CELL_LINE AND $CROSS_CELL_LINE ($K_MER-MER)\n\n\n"
        python train_test.py --cell_line="$CELL_LINE" --cross_cell_line="$CROSS_CELL_LINE" --k_mer=$K_MER --seed=$SEED
    done
done
