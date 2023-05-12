import pickle as pkl
import os
import numpy as np
import argparse

def loader(cell_line='GM12878', cross_cell_line='', label_rate=0.2, k_mer=5):

    if (cross_cell_line != None) and (cross_cell_line != cell_line):
        read_dir = 'data/{}_{}/'.format(cell_line, cross_cell_line)
    else:
        read_dir = 'data/{}/'.format(cell_line)

    # STEP 1: Load all feature vectors, class labels and graph
    features_file = open('{}/features_{}mer'.format(read_dir, k_mer), "rb")
    features = pkl.load(features_file)
    features_file.close()

    labels_file = open('{}/labels'.format(read_dir), "rb")
    labels = pkl.load(labels_file)
    labels_file.close()

    graph_file = open('{}/graph'.format(read_dir), "rb")
    graph = pkl.load(graph_file)
    graph_file.close()

    # STEP 2: Load IDs of labeled_train/unlabeled_train/validation/test nodes
    lr = txt = '{:.2f}'.format(label_rate).split('.')[1]

    idx_lx_file = open('{}/lx_{}.index'.format(read_dir, lr), "rb")
    idx_lx = pkl.load(idx_lx_file)
    idx_lx_file.close()

    idx_ux_file = open('{}/ux_{}.index'.format(read_dir, lr), "rb")
    idx_ux = pkl.load(idx_ux_file)
    idx_ux_file.close()

    idx_vx_file = open('{}/vx_{}.index'.format(read_dir, lr), "rb")
    idx_vx = pkl.load(idx_vx_file)
    idx_vx_file.close()

    idx_tx_file = open('{}/tx_{}.index'.format(read_dir, lr), "rb")
    idx_tx = pkl.load(idx_tx_file)
    idx_tx_file.close()

    # STEP 3: Take subsets from loaded features and class labels using loaded IDs
    idx_lx = np.sort(idx_lx)
    x = features[idx_lx]
    y = labels[idx_lx]

    idx_allx = np.concatenate((idx_lx, idx_ux, idx_vx))
    idx_allx = np.sort(idx_allx)

    allx = features[idx_allx]
    ally = labels[idx_allx]

    idx_tx = np.sort(idx_tx)
    tx = features[idx_tx]
    ty = labels[idx_tx]

    print("x={} allx={}".format(x.shape[0], allx.shape[0]))

    return graph, x, y, allx, ally, tx, ty, list(idx_tx)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='gcn4epi')
    parser.add_argument('--cell_line', default='GM12878', type=str)
    parser.add_argument('--cross_cell_line', default=None, type=str) # set to run cross cell-line testing
    parser.add_argument('--k_mer', default=5, type=int)
    parser.add_argument('--label_rate', default=0.2, type=float) # [0.2, 0.1, 0.05]
    args = parser.parse_args()


    graph, x, y, allx, ally, tx, ty, test_index = loader(args.cell_line, args.cross_cell_line, args.label_rate, args.k_mer)

    dump_dir = 'explanations/{}/Cora/raw/'.format(args.cell_line)

    if not os.path.exists(dump_dir):
        os.makedirs(dump_dir)

    file_names = ['x', 'y', 'tx', 'ty', 'allx', 'ally', 'graph']
    data = [x, y, tx, ty, allx, ally, graph]

    for f, d in zip(file_names, data):
        with open('{}/ind.{}.{}'.format(dump_dir, 'cora',f), "wb") as f:
            pkl.dump(d, f)
    
    with open('{}/ind.{}.test.index'.format(dump_dir, 'cora'), "w") as f:
        for item in test_index:
            f.write(str(item) + "\n")
