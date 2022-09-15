import os
import random
import pickle as pkl
import argparse
import numpy as np
from sklearn.model_selection import train_test_split
from prepare_data import getTuples, getLabels

def trainTestSplit(cell_line, cross_cell_line, id_dict, cross_begin_id, labels, label_rate, seed):

    def getIdPortions(cell_line, cross_cell_line, id_dict, cross_begin_id, labels, seed):

        """
            Returns ID portions for train, test, validation split.

            Label rate is the number of labeled nodes (x) that are used
            for training divided by the total number of nodes in dataset.

            Example for 0.2 label rate:
            20% tx (Fixed 20% of the same or cross cell-line)
            80% lx + ux + vx
                16% vx (Fixed 20% of the lx + ux + vx data)
                12.8% lx (20% of the lx + ux data since label rate is 0.2)
                51.2% ux (80% of the lx + ux data)
        """

        np.random.seed(seed)
        idx_lx_ux_vx_tx = list(id_dict.values())[0:cross_begin_id] # do not include cross cell-line elements
        scl_labels = labels[0:cross_begin_id] # do not include cross cell-line elements (scl = same cell line)
        idx_lx_ux_vx, idx_tx = train_test_split(idx_lx_ux_vx_tx, test_size=0.2, random_state=seed, stratify=scl_labels)
        idx_lx_ux, idx_vx = train_test_split(idx_lx_ux_vx, test_size=0.2, random_state=seed)
        idx_lx, idx_ux = train_test_split(idx_lx_ux, test_size=1-label_rate, random_state=seed)

        if cell_line == cross_cell_line or cross_cell_line == None:
            # No cross cell-line specified. Use same cell-line for testings.
            print('SAME CELL-LINE TESTING:\n{} labeled training \n{} validation \n{} test ({}) \n{} unlabeled training'
                .format(len(idx_lx), len(idx_vx), len(idx_tx), cell_line, len(idx_ux)))
        else:
            # Use cross cell-line for testing. Overwrite tx.
            cross_idx_lx_ux_vx_tx = list(id_dict.values())[cross_begin_id:]
            ccl_labels = labels[cross_begin_id:] # ccl = cross cell line
            _, idx_tx = train_test_split(cross_idx_lx_ux_vx_tx, test_size=0.2, random_state=seed, stratify=ccl_labels)
            print('CROSS CELL-LINE TESTING:\n{} labeled training \n{} validation \n{} test ({}) \n{} unlabeled training'
                .format(len(idx_lx), len(idx_vx), len(idx_tx), cross_cell_line, len(idx_ux)))

        return idx_lx, idx_ux, idx_vx, idx_tx


    # TRAIN / TEST / VALIDATION SPLIT
    idx_lx, idx_ux, idx_vx, idx_tx = getIdPortions(cell_line, cross_cell_line, id_dict, cross_begin_id, labels, seed)
    print('Writing index files for train/test/validation split...')

    if (args.cross_cell_line != None) and (args.cross_cell_line != args.cell_line):
        dump_dir = 'data/{}/'.format(cell_line + '_' + cross_cell_line)
    else:
        dump_dir = 'data/{}/'.format(cell_line)

    if not os.path.exists(dump_dir):
        os.makedirs(dump_dir)

    lr = '{:.2f}'.format(label_rate).split('.')[1]

    lx_file = open('{}/lx_{}.index'.format(dump_dir, lr), "wb")
    pkl.dump(idx_lx, lx_file)
    lx_file.close()

    ux_file = open('{}/ux_{}.index'.format(dump_dir, lr), "wb")
    pkl.dump(idx_ux, ux_file)
    ux_file.close()

    vx_file = open('{}/vx_{}.index'.format(dump_dir, lr), "wb")
    pkl.dump(idx_vx, vx_file)
    vx_file.close()

    tx_file = open('{}/tx_{}.index'.format(dump_dir, lr), "wb")
    pkl.dump(idx_tx, tx_file)
    tx_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='gcn4epi')
    parser.add_argument('--cell_line', default='GM12878', type=str)
    parser.add_argument('--cross_cell_line', default=None, type=str) # set to run cross cell-line testing
    parser.add_argument('--k_mer', default=5, type=int)
    parser.add_argument('--seed', default=42, type=int)
    parser.add_argument('--label_rate', default=0.2, type=float) # [0.2, 0.1, 0.05]
    parser.add_argument('--frag_len', default=200, type=int) # set 0 to disable fragmentation and use full sequences
    args = parser.parse_args()
    random.seed(args.seed)

    df_ep, id_dict, cross_begin_id = getTuples(args.cell_line, args.cross_cell_line, args.k_mer)  # requires successful run of prepare_gcn_data.py
    labels = getLabels(df_ep, len(id_dict))

    trainTestSplit(args.cell_line, args.cross_cell_line, id_dict, cross_begin_id, labels, args.label_rate, args.seed)
