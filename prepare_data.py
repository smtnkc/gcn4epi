import os
import sys
import random
import numpy as np
import pandas as pd
import pickle as pkl
import scipy.sparse as sp
import warnings
import argparse
import pcdhit
from Bio import SeqIO
from sklearn.feature_extraction.text import TfidfVectorizer
from collections import Counter
from progress.bar import Bar


def getSentences(cell_line, ep_sequences, k_mer):

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

    ep_sentences.to_csv('data/{}/ep_sentences_{}mer.csv'.format(cell_line, k_mer), index=False)
    print('EP sentences are written!')
    return ep_sentences


def getNodeById(df_ep, node_id):
    for row in range(len(df_ep)):
        enh = df_ep['enhancer'][row]
        pro = df_ep['promoter'][row]
        if enh[0] == node_id:
            return enh
        elif pro[0] == node_id:
            return pro


def getTuples(cell_line, cross_cell_line, k_mer):
    """
    Returns a new DF where each element is a tuple of 3 elements: (id, name, sequence)
    """
    ep_sentences = pd.read_csv('data/{}/ep_sentences_{}mer.csv'.format(cell_line, k_mer))
    if (cross_cell_line != None) and (cross_cell_line != cell_line):
        cross_ep_sentences = pd.read_csv('data/{}/ep_sentences_{}mer.csv'.format(cross_cell_line, k_mer))

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

    cross_begin_id = chr_id

    if (cross_cell_line != None) and (cross_cell_line != cell_line):
        for i in range(len(cross_ep_sentences)):
            e_list.append((cross_ep_sentences['enhancer_name'][i],
                           cross_ep_sentences['enhancer_sentence'][i]))
            p_list.append((cross_ep_sentences['promoter_name'][i],
                           cross_ep_sentences['promoter_sentence'][i]))

        cross_ep_list = sorted(list(set(list(cross_ep_sentences['enhancer_name']) + \
                                        list(cross_ep_sentences['promoter_name']))))

        # ADD CROSS CELL-LINE ENHANCERS AND PROMOTERS INTO ID_DICT
        for ep in cross_ep_list:
            id_dict[ep] = chr_id
            chr_id += 1

    for i in range(len(e_list)):
        e_list[i] = (id_dict[e_list[i][0]], ) + e_list[i]
        
    for i in range(len(p_list)):
        p_list[i] = (id_dict[p_list[i][0]], ) + p_list[i]

    df_ep = pd.DataFrame({'enhancer': e_list, 'promoter': p_list})

    return df_ep, id_dict, cross_begin_id


def getAdjMatrix(df_ep, node_count):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=sp.SparseEfficiencyWarning)
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


def getSequences(cell_line):

    def fetchPairs(cell_line):
        """
            If your cell line is not available in TargetFinder repo,
            Place your ep_pairs.csv file manually under your cell line directory.
        """
        available_cell_lines = ['GM12878', 'HUVEC', 'HeLa-S3', 'IMR90', 'K562', 'NHEK', 'combined']
        pairs_file = 'data/{}/ep_pairs.csv'.format(cell_line)

        if cell_line not in available_cell_lines:
            print('{} cell line is not in available.\nSelect one of {}\n' \
                'Or manually create {}'.format(cell_line, available_cell_lines, pairs_file))
            return None

        if os.path.isfile(pairs_file):
            print('Reading pairs from {}...'.format(pairs_file))
            ep_pairs = pd.read_csv(pairs_file)
        else:
            print('Reading pairs from remote github repo...')
            ep_pairs = pd.read_csv('https://raw.githubusercontent.com/shwhalen/' \
                                'targetfinder/master/paper/targetfinder/{}/' \
                                'output-ep/pairs.csv'.format(cell_line))
            if not os.path.isdir('data/{}'.format(cell_line)):
                print('Creating directory for {} cell line...'.format(cell_line))
                os.makedirs('data/{}'.format(cell_line))
            print('Writing pairs to {}'.format(pairs_file))
            ep_pairs.to_csv(pairs_file, index=False)
        return ep_pairs

    def fetchSequences(ep_pairs):
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


    seq_file = 'data/{}/ep_sequences.csv'.format(cell_line)
    if os.path.isfile(seq_file):
        print('Reading existing sequences from {}...'.format(seq_file))
        ep_sequences = pd.read_csv(seq_file)
    else:
        ep_pairs = fetchPairs(cell_line)
        if (ep_pairs is None):
            sys.exit()
        print('{} EP pairs have been read.'.format(len(ep_pairs)))
        ep_pairs = ep_pairs[ep_pairs['label'] == 1].reset_index() # Keep only the interacting pairs
        print('{} EP pairs are labeled as 1.'.format(len(ep_pairs)))
        ep_sequences = fetchSequences(ep_pairs)
        ep_sequences.to_csv(seq_file, index=False)
        print('EP sequences are written!')
    return ep_sequences


def getFragments(cell_line, frag_len, balanced):

    def generateFragments(ep_sequences, frag_len):
        enh_names = []
        enh_frag_names = []
        enh_frag_seqs = []
        for i in range(len(ep_sequences)):
            seq = ep_sequences['enhancer_seq'][i]
            name = ep_sequences['enhancer_name'][i]
            coordinates = name.split(':')[1]
            coor_start = int(coordinates.split('-')[0])
            coor_end = coor_start + frag_len
            while len(seq) >= frag_len:
                fragment = str(coor_start) + '-' + str(coor_end)
                enh_names.append(name)
                enh_frag_names.append(name.split(':')[0] + ':' + fragment)
                enh_frag_seqs.append(seq[:frag_len])
                seq = seq[frag_len:]
                coor_start = coor_end
                coor_end = coor_start + frag_len

        pro_names = []
        pro_frag_names = []
        pro_frag_seqs = []
        for i in range(len(ep_sequences)):
            seq = ep_sequences['promoter_seq'][i]
            name = ep_sequences['promoter_name'][i]
            coordinates = name.split(':')[1]
            coor_start = int(coordinates.split('-')[0])
            coor_end = coor_start + frag_len
            while len(seq) >= frag_len:
                fragment = str(coor_start) + '-' + str(coor_end)
                pro_names.append(name)
                pro_frag_names.append(name.split(':')[0] + ':' + fragment)
                pro_frag_seqs.append(seq[:frag_len])
                seq = seq[frag_len:]
                coor_start = coor_end
                coor_end = coor_start + frag_len

        df_enh_fragments = pd.DataFrame({'enhancer_name': enh_names, 'enhancer_frag_name': enh_frag_names, 'enhancer_frag_seq': enh_frag_seqs})
        df_pro_fragments = pd.DataFrame({'promoter_name': pro_names, 'promoter_frag_name': pro_frag_names, 'promoter_frag_seq': pro_frag_seqs})

        df_enh_fragments = df_enh_fragments.drop_duplicates(subset=['enhancer_frag_name']).reset_index(drop=True)
        df_pro_fragments = df_pro_fragments.drop_duplicates(subset=['promoter_frag_name']).reset_index(drop=True)
        return df_enh_fragments, df_pro_fragments

    def filterFragments(df_frags, threshold):
        '''
            Filters out the fragments with similarity higher than a specified threshold
        '''
        filtered_frags = list(pcdhit.filter(list(zip(df_frags.iloc[:,1], df_frags.iloc[:,2])), threshold=threshold))
        df_ff = df_frags[df_frags.iloc[:,2].isin([e[1] for e in filtered_frags])].reset_index(drop=True)
        return df_ff

    def mergeFragments(ep_sequences, df_fef, df_fpf):
        col_names =  ['enhancer_name', 'enhancer_frag_name', 'enhancer_frag_seq',
                'promoter_name', 'promoter_frag_name', 'promoter_frag_seq']

        merged_df = pd.DataFrame(columns = col_names)

        with Bar('Processing', max=len(ep_sequences)) as bar:
            for i in range(len(ep_sequences)):
                enh_frags = df_fef[df_fef['enhancer_name'] == ep_sequences['enhancer_name'][i]]
                pro_frags = df_fpf[df_fpf['promoter_name'] == ep_sequences['promoter_name'][i]]

                for e in range(len(enh_frags)):
                    for p in range(len(pro_frags)):
                        e_row = enh_frags[e:e+1].reset_index(drop=True)
                        p_row = pro_frags[p:p+1].reset_index(drop=True)
                        merged_row = pd.concat([e_row, p_row], axis=1)
                        merged_df = pd.concat([merged_df, merged_row])
                bar.next()

        return merged_df.reset_index(drop=True)

    def getBalancedDf(df, cell_line):
        # To balance the fragments, we use most frequent promoters or enhancers
        # For example, 3189 is selected for GM12878 after several trials

        balance_cutoffs = {
            'GM12878': 3189,
            'HUVEC': 3522,
            'HeLa-S3': 1771,
            'IMR90': 218,
            'K562': 1277,
            'NHEK': 32,
            'combined': 9903
        }

        n_enh = len(set(df['enhancer_frag_name']))
        n_pro = len(set(df['promoter_frag_name']))

        if n_pro > n_enh:
            most_freq_promoters = [p[0] for p in Counter(df['promoter_frag_name']).most_common(balance_cutoffs[cell_line])]
            df_balanced = df[df['promoter_frag_name'].isin(most_freq_promoters)].reset_index(drop=True)
        else:
            most_freq_enhancers = [p[0] for p in Counter(df['enhancer_frag_name']).most_common(balance_cutoffs[cell_line])]
            df_balanced = df[df['enhancer_frag_name'].isin(most_freq_enhancers)].reset_index(drop=True)

        return df_balanced

    frag_path = 'data/{}/frag_pairs{}_{}.csv'.format(cell_line, '_balanced' if balanced else '', frag_len)

    if os.path.isfile(frag_path):
        print('Reading fragments from {}...'.format(frag_path))
        ep_frags = pd.read_csv(frag_path)
        ep_frags = ep_frags[['enhancer_frag_name', 'enhancer_frag_seq', 'promoter_frag_name', 'promoter_frag_seq']]
        ep_frags.columns = ['enhancer_name', 'enhancer_seq', 'promoter_name', 'promoter_seq']
        print('{} enhancer fragments.'.format(len(set(ep_frags['enhancer_name']))))
        print('{} promoter fragments.'.format(len(set(ep_frags['promoter_name']))))
        print('{} interactions between EP fragments.'.format(len(ep_frags)))
    else:
        print('Generating fragments from scratch...')
        ep_sequences = getSequences(cell_line)

        print('Removing sequences shorter than {}bp...'.format(frag_len))
        ep_sequences = ep_sequences[
            ep_sequences['enhancer_seq'].apply(lambda x: len(x)>=frag_len) &
            ep_sequences['promoter_seq'].apply(lambda x: len(x)>=frag_len)].reset_index(drop=True)

        print('{} enhancers with length >= {}'.format(len(set(ep_sequences['enhancer_name'])), frag_len))
        print('{} promoters with length >= {}'.format(len(set(ep_sequences['promoter_name'])), frag_len))
        print('SPLITTING INTO FRAGMENTS...')
        df_enh_frags, df_pro_frags = generateFragments(ep_sequences, frag_len)
        print('{} fragments from {} enhancers.'.format(len(df_enh_frags), len(set(df_enh_frags['enhancer_name']))))
        print('{} fragments from {} promoters.'.format(len(df_pro_frags), len(set(df_pro_frags['promoter_name']))))
        df_fef = filterFragments(df_enh_frags, 0.8)  # filter out if similarity is higher than 80%
        df_fpf = filterFragments(df_pro_frags, 0.8)  # filter out if similarity is higher than 80%

        df_merged_frags = mergeFragments(ep_sequences, df_fef, df_fpf)
        df_merged_frags.to_csv('data/{}/frag_pairs_{}.csv'.format(cell_line, frag_len), index=False)

        df_merged_frags_balanced = getBalancedDf(df_merged_frags, cell_line)
        print('{} enhancer fragments with low similarity.'.format(len(set(df_merged_frags_balanced['enhancer_frag_name']))))
        print('{} promoter fragments with low similarity.'.format(len(set(df_merged_frags_balanced['promoter_frag_name']))))

        df_merged_frags_balanced.to_csv('data/{}/frag_pairs_balanced_{}.csv'.format(cell_line, frag_len), index=False)

        if balanced:
            ep_frags = df_merged_frags_balanced
        else:
            ep_frags = df_merged_frags

        ep_frags = ep_frags[['enhancer_frag_name', 'enhancer_frag_seq', 'promoter_frag_name', 'promoter_frag_seq']]
        ep_frags.columns = ['enhancer_name', 'enhancer_seq', 'promoter_name', 'promoter_seq']
    return ep_frags


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='gcn4epi')
    parser.add_argument('--cell_line', default='GM12878', type=str)
    parser.add_argument('--cross_cell_line', default=None, type=str) # set to run cross cell-line testing
    parser.add_argument('--k_mer', default=5, type=int)
    parser.add_argument('--label_rate', default=0.2, type=float) # [0.2, 0.1, 0.05]
    parser.add_argument('--frag_len', default=200, type=int) # set 0 to disable fragmentation and use full sequences
    parser.add_argument('--balanced', action='store_true') # set to balance enhancers and promoters
    args = parser.parse_args()

    if args.frag_len > 0:
        # Use fix-sized fragments (not full sequences)
        ep_sequences = getFragments(args.cell_line, args.frag_len, args.balanced)
        if (args.cross_cell_line != None) and (args.cross_cell_line != args.cell_line):
            cross_ep_sequences = getFragments(args.cross_cell_line, args.frag_len, args.balanced)
    else:
        # Use full sequences (not fragments)
        ep_sequences = getSequences(args.cell_line)
        if (args.cross_cell_line != None) and (args.cross_cell_line != args.cell_line):
            cross_ep_sequences = getSequences(args.cross_cell_line)

    ep_sentences = getSentences(args.cell_line, ep_sequences, args.k_mer)  # also writes EP sentences to files
    if (args.cross_cell_line != None) and (args.cross_cell_line != args.cell_line):
        cross_ep_sentences = getSentences(args.cross_cell_line, cross_ep_sequences, args.k_mer)  # also writes EP sentences to files

    df_ep, id_dict, cross_begin_id = getTuples(args.cell_line, args.cross_cell_line, args.k_mer)  # requires successful run of getSentences()

    if (args.cross_cell_line != None) and (args.cross_cell_line != args.cell_line):
        dump_dir = 'data/{}/'.format(args.cell_line + '_' + args.cross_cell_line)
    else:
        dump_dir = 'data/{}/'.format(args.cell_line)

    if not os.path.exists(dump_dir):
        os.makedirs(dump_dir)

    nodes_file = open('{}/nodes'.format(dump_dir), "wb")
    pkl.dump(id_dict, nodes_file)
    nodes_file.close()

    adj = getAdjMatrix(df_ep, node_count=len(id_dict))

    print('Writing adjacency matrix...')
    graph = {i: np.nonzero(row)[1].tolist() for i,row in enumerate(adj)}
    graph_file = open('{}/graph'.format(dump_dir), "wb")
    pkl.dump(graph, graph_file)
    graph_file.close()

    features = getFeatureVectors(df_ep)
    print('Writing feature vectors...')
    features_file = open('{}/features_{}mer'.format(dump_dir, args.k_mer), "wb")
    pkl.dump(features, features_file)
    features_file.close()

    labels = getLabels(df_ep, len(id_dict))
    print('Writing binary class labels...')
    labels_file = open('{}/labels'.format(dump_dir), "wb")
    pkl.dump(labels, labels_file)
    labels_file.close()
