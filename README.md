# INSIDER #

This repository contains scripts for detecting foreign DNA sequences in genomes.

## Requirements ##

python >= 3.7.0, pyspark >= 3.0.0, scikit-learn >= 0.24.0, scipy >= 1.6.0, statsmodels >= 0.12.0

## Usage ##

### Calculate K-mer frequencies

```
python calculate_kmer_frequencies.py \
    split \
    -f test_file.fa \
    -k 2 \
    -n \
    -o test_2mer
```

For each sequence, extract 2-mers and count their frequencies. Ambigous bases (i.e., N's) are ignored)

### Cluster K-mer frequencies

```
python cluster_kmer_frequencies.py \
    consensus \
    --obs test_2mer \
    --params params.json \
    -o test_2mer_cIds.txt
```

Cluster sequences based on their K-mers. Hyperparameters can be specified in the JSON file.

### Analyse K-mer frequencies

```
python analyse_kmer_frequencies.py
    main \
    --obs test_2mer \
    --cid test_2mer_cIds.txt \
    -o test_2mer_output.txt 
```
Assess the similarity between each cluster and the genome based on their K-mer frequencies.
