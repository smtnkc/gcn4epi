from Bio import SeqIO
import pandas as pd
import os
import argparse


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
    return ep_pairs

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='gcn4epi')
    parser.add_argument('--cell_line', default='GM12878', type=str)
    parser.add_argument('--k_mer', default=5, type=int)
    args = parser.parse_args()

    ep_pairs = getPairs(args.cell_line)
    if (ep_pairs is not None):
        print('{} EP pairs have been read.'.format(len(ep_pairs)))
        if not os.path.isdir(args.cell_line):
            print('Creating directory for {} cell line...'.format(args.cell_line))
            os.makedirs(args.cell_line)
        ep_pairs.to_csv('data/{}/ep_pairs.csv'.format(args.cell_line), index=False)

        ep_pairs = ep_pairs[ep_pairs['label'] == 1]
        print('{} EP pairs are labeled as 1.'.format(len(ep_pairs)))

        ep_sequences = getSequences(ep_pairs)
        ep_sequences.to_csv('data/{}/ep_sequences.csv'.format(args.cell_line), index=False)
        print('EP sequences are written!')

        ep_sentences = getSentences(ep_sequences, args.k_mer)
        ep_sentences.to_csv('data/{}/ep_sentences_{}mer.csv'.format(
            args.cell_line, args.k_mer), index=False)
        print('EP sentences are written!')
