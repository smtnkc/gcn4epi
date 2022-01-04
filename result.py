import os
import pandas as pd


files = []
seeds = []


for seed in os.listdir('results'):
    if os.path.isdir('results/' + seed) and not seed.startswith('.'):
        seeds.append(seed)


for file in os.listdir('results/' + seeds[0]):
    if file.endswith(".txt"):
        files.append(file)

files = sorted(files)

aucs = [[0.0 for i in range(len(files))] for j in range(len(seeds))]

for i in range(len(seeds)):
    train_cell_lines = []
    test_cell_lines = []
    for j in range(len(files)):
        f = open('results/' + seeds[i] + '/' + files[j], "r")
        auc = float(f.read().split('\n')[-4].split('= ')[1])
        f.close()

        cell_lines = files[j].split('.')[0].split('_')
        train_cell_lines.append(cell_lines[0])
        if len(cell_lines) == 1:
            test_cell_lines.append(cell_lines[0])
        else:
            test_cell_lines.append(cell_lines[1])
        aucs[i][j] = float(auc)

data = {
    'train_cell_line': train_cell_lines,
    'test_cell_line': test_cell_lines
}

for i in range(len(seeds)):
    auc_key = 'auc_' + seeds[i]
    data[auc_key] = aucs[i]
    print(len(data[auc_key]))

print(len(data['train_cell_line']))
print(len(data['test_cell_line']))


df = pd.DataFrame(data)
df.to_csv('results/results.csv', index=False)