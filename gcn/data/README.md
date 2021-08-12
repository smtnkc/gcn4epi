# Requirements

1. Download Human Genome **GRCh37** from [Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml) and place it under this directory. **Example:** `data/GRCh37_latest_genomic.fna`

2. Download `pairs.csv` file for desired cell line from [TargetFinder Repo](https://github.com/shwhalen/targetfinder/tree/master/paper/targetfinder). Rename it as `ep_pairs.csv`. Place it under a subdirectory with the same name as its cell line. **Example:** `data/GM12878/ep_pairs.csv`

3. Run `fetch_and_parse` and `prepare_gcn_data` modules respectively to prepare data files required by GCN.

![S1](./s1.png)

| **File Path** | **Description** |
| :-- | :-- |
| GM12878/x.index  | the indices of labeled train instances as list object |
| GM12878/ux.index | the indices of unlabeled train instances as list object |
| GM12878/vx.index | the indices of validation instances as list object |
| GM12878/tx.index | the indices of test instances as list object |
| GM12878/labels   | the one-hot labels of **all** instances as numpy.ndarray object |
| GM12878/features | the feature vectors of **all** instances as scipy.sparse.csr.csr_matrix object |
| GM12878/nodes    | a dict in the format **{id: chromosome_name}** as collections.defaultdict object |
| GM12878/graph | a dict in the format **{id: [id_of_neighbor_nodes]}** as collections.defaultdict object |
