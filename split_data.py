import os
import random
import pickle as pkl
import argparse
from sklearn.model_selection import train_test_split
from prepare_data import getTuples

def trainTestSplit(cell_line, cross_cell_line, id_dict, cross_begin_id, label_rate, seed):

    def getIdPortions(cell_line, cross_cell_line, id_dict, cross_begin_id, seed):

        """
            Returns ID portions for train, test, validation split.

            Label rate is the number of labeled nodes (x) that are used
            for training divided by the total number of nodes in dataset.

            Example: Label rate = 0.1
            10% labeled training (x)
            60% unlabaled training (ux)
            10% validation (vx)
            20% test (tx) !!! 20% of the same or cross cell-line !!!

            allx = x + ux + vx
        """

        idx = list(id_dict.values())[0:cross_begin_id] # do not include cross cell-line elements
        idx_allx, idx_tx = train_test_split(idx, test_size=0.2, random_state=seed)
        idx_x_vx, idx_ux = train_test_split(idx_allx, test_size=1-(label_rate*2/0.8), random_state=seed)
        idx_x, idx_vx = train_test_split(idx_x_vx, test_size=0.5, random_state=seed)

        if cross_begin_id == len(id_dict):
            # No cross cell-line specified. Use same cell-line for testings.
            print('SAME CELL-LINE TESTING:\n {} labeled training \n {} validation \n {} test ({}) \n{} unlabeled training'
                .format(len(idx_x), len(idx_vx), len(idx_tx), cell_line, len(idx_ux)))
        else:
            # Use cross cell-line for testing. Overwrite idx_tx.
            cross_idx = list(id_dict.values())[cross_begin_id:]
            _, idx_tx = train_test_split(cross_idx, test_size=0.2, random_state=seed)
            print('CROSS CELL-LINE TESTING:\n {} labeled training \n {} validation \n {} test ({}) \n{} unlabeled training'
                .format(len(idx_x), len(idx_vx), len(idx_tx), cross_cell_line, len(idx_ux)))

        return idx_x, idx_ux, idx_vx, idx_tx


    # TRAIN / TEST / VALIDATION SPLIT
    idx_x, idx_ux, idx_vx, idx_tx = getIdPortions(cell_line, cross_cell_line, id_dict, cross_begin_id, seed)
    print('Writing index files for train/test/validation split...')

    if (args.cross_cell_line != None) and (args.cross_cell_line != args.cell_line):
        dump_dir = 'data/{}/'.format(cell_line + '_' + cross_cell_line)
    else:
        dump_dir = 'data/{}/'.format(cell_line)

    if not os.path.exists(dump_dir):
        os.makedirs(dump_dir)

    lr = '{:.2f}'.format(label_rate).split('.')[1]

    idx_x_file = open('{}/x_{}.index'.format(dump_dir, lr), "wb")
    pkl.dump(idx_x, idx_x_file)
    idx_x_file.close()

    idx_ux_file = open('{}/ux_{}.index'.format(dump_dir, lr), "wb")
    pkl.dump(idx_ux, idx_ux_file)
    idx_ux_file.close()

    idx_vx_file = open('{}/vx_{}.index'.format(dump_dir, lr), "wb")
    pkl.dump(idx_vx, idx_vx_file)
    idx_vx_file.close()

    idx_tx_file = open('{}/tx_{}.index'.format(dump_dir, lr), "wb")
    pkl.dump(idx_tx, idx_tx_file)
    idx_tx_file.close()

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

    _, id_dict, cross_begin_id = getTuples(args.cell_line, args.cross_cell_line, args.k_mer)  # requires successful run of prepare_gcn_data.py

    trainTestSplit(args.cell_line, args.cross_cell_line, id_dict, cross_begin_id, args.label_rate, args.seed)
