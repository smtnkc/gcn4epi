import os
import pandas as pd


files = []
seed_kmer = []


for sub_dir in os.listdir('results'):
    if os.path.isdir('results/' + sub_dir) and not sub_dir.startswith('.'):
        seed_kmer.append(sub_dir)


for file in os.listdir('results/' + seed_kmer[0]):
    if file.endswith(".txt"):
        files.append(file)

files = sorted(files)
print(files)

f1_scores = [[0.0 for i in range(len(files))] for j in range(len(seed_kmer))]

for i in range(len(seed_kmer)):
    train_cell_lines = []
    test_cell_lines = []
    for j in range(len(files)):
        f_path = 'results/' + seed_kmer[i] + '/' + files[j]
        f = open(f_path, "r")
        f1 = f.read().split('Test F1')[1].split('\n')[0].split('= ')[1]
        print(f_path, '=', f1)
        f1 = float(f1)
        f.close()

        cell_lines = files[j].split('.')[0].split('_')
        train_cell_lines.append(cell_lines[0])
        if len(cell_lines) == 1:
            test_cell_lines.append(cell_lines[0])
        else:
            test_cell_lines.append(cell_lines[1])
        f1_scores[i][j] = float(f1)

data = {
    'train_cell_line': train_cell_lines,
    'test_cell_line': test_cell_lines
}

for i in range(len(seed_kmer)):
    f1_key = 'f1_' + seed_kmer[i]
    data[f1_key] = f1_scores[i]
    print('Seed {} f1 count {}'.format(seed_kmer[i], len(data[f1_key])))


df = pd.DataFrame(data)
df.to_csv('results/results.csv', index=False)
