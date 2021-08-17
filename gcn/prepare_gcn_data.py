import os
import sys
import random
import numpy as np
import pandas as pd
import pickle as pkl
import scipy.sparse as sp
import argparse
from Bio import SeqIO
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split


def getPairs(cell_line):
    """
        If your cell line is not available in TargetFinder repo,
        Place your ep_pairs.csv file manually under your cell line directory.
    """
    available_cell_lines = ['GM12878', 'HUVEC', 'HeLa-S3', 'IMR90', 'K562', 'NHEK', 'combined']

    if cell_line not in available_cell_lines:
        print('{} cell line is not in available.\nSelect one of {}\n' \
              'Or manually create data/{}/ep_pairs.csv'.format(cell_line, available_cell_lines, cell_line))
        return None

    if os.path.isfile('data/{}/ep_pairs.csv'.format(cell_line)):
        print('Reading pairs from local file...')
        ep_pairs = pd.read_csv('data/{}/ep_pairs.csv'.format(cell_line))
    else:
        print('Reading pairs from remote github repo...')
        ep_pairs = pd.read_csv('https://raw.githubusercontent.com/shwhalen/' \
                               'targetfinder/master/paper/targetfinder/{}/' \
                               'output-ep/pairs.csv'.format(cell_line))
        if not os.path.isdir('data/{}'.format(args.cell_line)):
            print('Creating directory for {} cell line...'.format(args.cell_line))
            os.makedirs('data/{}'.format(args.cell_line))
        print('Writing pairs to data/{}/ep_pairs.csv'.format(args.cell_line))
        ep_pairs.to_csv('data/{}/ep_pairs.csv'.format(args.cell_line), index=False)
    return ep_pairs

def getSequences(ep_pairs):
    # DOWNLOAD HUMAN GENOME v37 (3.2 Gb)
    # Older version but compatible with genomic coordinates of TargetFinder dataset
    # https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml
    # https://github.com/shwhalen/targetfinder/tree/master/paper/targetfinder

    print('Parsing GRCh37 genome...')
    hg37 = SeqIO.to_dict(SeqIO.parse('data/GRCh37_latest_genomic.fna', 'fasta'))

    RefSeqIDs = []

    for k in hg37.keys():
        if k.startswith('NC_0000'):
            RefSeqIDs.append(hg37[k].id)

    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
               'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', \
               'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    RefSeqDict = {chromosomes[i]: RefSeqIDs[i] for i in range(len(chromosomes))}

    enhancer_sequences = []
    promoter_sequences = []
    n = len(ep_pairs)

    print('Getting DNA sequences for {} EP pairs...'.format(n))

    for i in range(n):
        enhancer_seq_id = ep_pairs['enhancer_chrom'][i]
        enhancer_seq_start = ep_pairs['enhancer_start'][i] - 1
        enhancer_seq_end = ep_pairs['enhancer_end'][i]

        promoter_seq_id = ep_pairs['promoter_chrom'][i]
        promoter_seq_start = ep_pairs['promoter_start'][i] - 1
        promoter_seq_end = ep_pairs['promoter_end'][i]
        
        enhancer_sequences.append(str(hg37[RefSeqDict[enhancer_seq_id]]
                                    .seq[enhancer_seq_start:enhancer_seq_end]).upper())

        promoter_sequences.append(str(hg37[RefSeqDict[promoter_seq_id]]
                                    .seq[promoter_seq_start:promoter_seq_end]).upper())

    ep_sequences = pd.DataFrame({'enhancer_name': ep_pairs['enhancer_name'][0:n],
                                 'promoter_name': ep_pairs['promoter_name'][0:n],
                                 'enhancer_seq': enhancer_sequences,
                                 'promoter_seq': promoter_sequences})
    return ep_sequences

def DNA2Sentence(dna, K, clean=False):
    if clean:
        dna = dna.replace("N", "")

    sentence = ""
    length = len(dna)

    for i in range(length - K + 1):
        sentence += dna[i: i + K] + " "

    # remove spaces
    sentence = sentence[0 : len(sentence) - 1]
    return sentence

def getSentences(ep_sequences, k_mer):
    enhancer_sentences = []
    promoter_sentences = []
    n = len(ep_sequences)

    print('Creating {}-mer sentences for {} EP pairs...'.format(k_mer, n))

    for i in range(len(ep_sequences)):
        enhancer_sentences.append(DNA2Sentence(ep_sequences['enhancer_seq'][i], k_mer))
        promoter_sentences.append(DNA2Sentence(ep_sequences['promoter_seq'][i], k_mer))

    ep_sentences = pd.DataFrame({'enhancer_name': ep_sequences['enhancer_name'][0:n],
                                 'promoter_name': ep_sequences['promoter_name'][0:n],
                                 'enhancer_sentence': enhancer_sentences,
                                 'promoter_sentence': promoter_sentences})
    return ep_sentences

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

    ep_pairs = getPairs(args.cell_line)
    if (ep_pairs is None):
        sys.exit()

    print('{} EP pairs have been read.'.format(len(ep_pairs)))

    ep_pairs = ep_pairs[ep_pairs['label'] == 1] # Keep only the interacting pairs
    print('{} EP pairs are labeled as 1.'.format(len(ep_pairs)))

    ep_sequences = getSequences(ep_pairs)
    ep_sequences.to_csv('data/{}/ep_sequences.csv'.format(args.cell_line), index=False)
    print('EP sequences are written!')

    ep_sentences = getSentences(ep_sequences, args.k_mer)
    ep_sentences.to_csv('data/{}/ep_sentences_{}mer.csv'.format(
        args.cell_line, args.k_mer), index=False)
    print('EP sentences are written!')

    df_ep, id_dict = getTuples(args)

    adj = getAdjMatrix(df_ep, node_count=len(id_dict))
    print('Writing adjacency matrix...')
    graph = {i: np.nonzero(row)[1].tolist() for i,row in enumerate(adj)}
    graph_file = open('data/{}/graph'.format(args.cell_line), "wb")
    pkl.dump(graph, graph_file)
    graph_file.close()

    features = getFeatureVectors(df_ep)
    print('Writing feature vectors...')
    features_file = open('data/{}/features_{}mer'.format(args.cell_line, args.k_mer), "wb")
    pkl.dump(features, features_file)
    features_file.close()

    labels = getLabels(df_ep, len(id_dict))
    print('Writing binary class labels...')
    labels_file = open('data/{}/labels'.format(args.cell_line), "wb")
    pkl.dump(labels, labels_file)
    labels_file.close()

    # TRAIN / TEST / VALIDATION SPLIT
    idx_x, idx_ux, idx_vx, idx_tx = getIdPortions(id_dict, args)
    print('Writing index files for train/test/validation split...')

    lr = txt = '0' + '{:.2f}'.format(args.label_rate).split('.')[1]

    idx_x_file = open('data/{}/x_{}.index'.format(args.cell_line, lr), "wb")
    pkl.dump(idx_x, idx_x_file)
    idx_x_file.close()

    idx_ux_file = open('data/{}/ux_{}.index'.format(args.cell_line, lr), "wb")
    pkl.dump(idx_ux, idx_ux_file)
    idx_ux_file.close()

    idx_vx_file = open('data/{}/vx_{}.index'.format(args.cell_line, lr), "wb")
    pkl.dump(idx_vx, idx_vx_file)
    idx_vx_file.close()

    idx_tx_file = open('data/{}/tx_{}.index'.format(args.cell_line, lr), "wb")
    pkl.dump(idx_tx, idx_tx_file)
    idx_tx_file.close()
