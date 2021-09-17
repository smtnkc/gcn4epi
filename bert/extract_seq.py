import os
import re

full_path = os.path.realpath(__file__)
os.chdir(os.path.dirname(full_path))

print(f"Change CWD to: {os.path.dirname(full_path)}")

def extract_seq(file_path, dir_path='seq', seq_length=200):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    nseq = 0
    data = re.split(r'(^>.*)', ''.join(open(file_path).readlines()), flags=re.M)

    for i in range(2, len(data), 2):
        fid = data[i-1][1:].split(':')
        nseq = nseq + 1
        fasta = data[i].replace('\n', '').replace(' ', '')
        seq = ' '.join(fasta)
        ffas = open(f"{dir_path}/{fid[0] + '_' + fid[1]}.seq", "w")
        ffas.write(seq)

    print(f"Number of sequences for {dir_path} -> {nseq}")

extract_seq('data/enhancer.cv.txt', 'cv_enh')
extract_seq('data/promoter.cv.txt', 'cv_pro')

extract_seq('data/enhancer.ind.txt', 'ind_enh')
extract_seq('data/promoter.ind.txt', 'ind_pro')

extract_seq('data/enhancer.all.txt', 'all_enh')
extract_seq('data/promoter.all.txt', 'all_pro')
