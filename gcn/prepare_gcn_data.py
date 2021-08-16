import sys
import random
import numpy as np
import pandas as pd
import pickle as pkl
import scipy.sparse as sp
import argparse
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split

def getNodeById(df_ep, node_id):
    for row in range(len(df_ep)):
        enh = df_ep['enhancer'][row]
        pro = df_ep['promoter'][row]
        if enh[0] == node_id:
            return enh
        elif pro[0] == node_id:
            return pro

def getTuples(args):
    """
    Returns a new DF where each element is a tuple of 3 elements: (id, name, sequence)
    """
    ep_sentences = pd.read_csv('data/{}/ep_sentences_{}mer.csv'.format(args.cell_line, args.k_mer))

    e_list = []
    p_list = []

    for i in range(len(ep_sentences)):
        e_list.append((ep_sentences['enhancer_name'][i],
                       ep_sentences['enhancer_sentence'][i]))

        p_list.append((ep_sentences['promoter_name'][i],
                       ep_sentences['promoter_sentence'][i]))

    ep_list = sorted(list(set(list(ep_sentences['enhancer_name']) + \
                              list(ep_sentences['promoter_name']))))

    # CREATE ID_DICT
    id_dict = {}
    chr_id = 0
    for ep in ep_list:
        id_dict[ep] = chr_id
        chr_id += 1

    # DUMP ID_DICT
    nodes_file = open('data/{}/nodes'.format(args.cell_line), "wb")
    pkl.dump(id_dict, nodes_file)
    nodes_file.close()

    for i in range(len(e_list)):
        e_list[i] = (id_dict[e_list[i][0]], ) + e_list[i]
        
    for i in range(len(p_list)):
        p_list[i] = (id_dict[p_list[i][0]], ) + p_list[i]

    df_ep = pd.DataFrame({'enhancer': e_list, 'promoter': p_list})
    return df_ep, id_dict

def getAdjMatrix(df_ep, node_count):
    adj = sp.csr_matrix((node_count, node_count), dtype=np.int32)

    for i in range(len(df_ep)):
        x = df_ep['enhancer'][i][0]
        y = df_ep['promoter'][i][0]
        adj[x,y] = 1
        adj[y,x] = 1

    return adj

def getFeatureVectors(df_ep):
    merged_list = list(set(list(df_ep['enhancer']) + list(df_ep['promoter'])))
    merged_list = sorted(merged_list) # sort by first element (id)

    corpus = []
    for t in merged_list:
        corpus.append(t[2])

    vectorizer = TfidfVectorizer()
    features = vectorizer.fit_transform(corpus)
    return features

def getLabels(df_ep, node_count):
    labels = np.zeros(shape=(node_count,2), dtype=np.int8) # values from -128 to 127

    for i in range(len(df_ep)):
        eid = df_ep['enhancer'][i][0]
        pid = df_ep['promoter'][i][0]
        labels[eid] = [1,0] # enhancer class
        labels[pid] = [0,1] # promoter class

    return labels

def getIdPortions(id_dict, args):

    """
        Returns ID portions for train, test, validation split.

        Label rate is the number of labeled nodes (x) that are used
        for training divided by the total number of nodes in dataset.

        Example: Label rate = 0.1
        10% labeled training (x)
        60% unlabaled training (ux)
        10% validation (vx)
        20% test (tx)

        allx = x + ux + vx
    """

    idx = list(id_dict.values())
    idx_allx, idx_tx = train_test_split(idx, test_size=0.2, random_state=args.seed)
    idx_x_vx, idx_ux = train_test_split(idx_allx, test_size=1-(args.label_rate*2/0.8),
                                        random_state=args.seed)
    idx_x, idx_vx = train_test_split(idx_x_vx, test_size=0.5, random_state=args.seed)

    print(' {} labeled training \n {} validation \n {} test \n{} unlabeled training'
        .format(len(idx_x), len(idx_vx), len(idx_tx), len(idx_ux)))

    return idx_x, idx_ux, idx_vx, idx_tx

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='gcn4epi')
    parser.add_argument('--cell_line', default='GM12878', type=str)
    parser.add_argument('--k_mer', default=5, type=int)
    parser.add_argument('--seed', default=42, type=int)
    parser.add_argument('--label_rate', default=0.2, type=float) # [0.2, 0.1, 0.05]
    args = parser.parse_args()
    random.seed(args.seed)

    df_ep, id_dict = getTuples(args)

    adj = getAdjMatrix(df_ep, node_count=len(id_dict))
    print('Writing adjacency matrix...')
    graph = {i: np.nonzero(row)[1].tolist() for i,row in enumerate(adj)}
    graph_file = open('data/{}/graph'.format(args.cell_line), "wb")
    pkl.dump(graph, graph_file)
    graph_file.close()

    features = getFeatureVectors(df_ep)
    print('Writing feature vectors...')
    features_file = open('data/{}/features'.format(args.cell_line), "wb")
    pkl.dump(features, features_file)
    features_file.close()

    labels = getLabels(df_ep, len(id_dict))
    print('Writing binary class labels...')
    labels_file = open('data/{}/labels'.format(args.cell_line), "wb")
    pkl.dump(labels, labels_file)
    labels_file.close()

    idx_x, idx_ux, idx_vx, idx_tx = getIdPortions(id_dict, args)
    print('Writing index files for train/test/validation split...')

    idx_x_file = open('data/{}/x.index'.format(args.cell_line), "wb")
    pkl.dump(idx_x, idx_x_file)
    idx_x_file.close()

    idx_ux_file = open('data/{}/ux.index'.format(args.cell_line), "wb")
    pkl.dump(idx_ux, idx_ux_file)
    idx_ux_file.close()

    idx_vx_file = open('data/{}/vx.index'.format(args.cell_line), "wb")
    pkl.dump(idx_vx, idx_vx_file)
    idx_vx_file.close()

    idx_tx_file = open('data/{}/tx.index'.format(args.cell_line), "wb")
    pkl.dump(idx_tx, idx_tx_file)
    idx_tx_file.close()
