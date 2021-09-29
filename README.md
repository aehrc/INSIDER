# INSIDER #
This repository contains scripts for detecting foreign DNA sequences in genomes.

## Requirements ##
python >= 3.7.0, pyspark >= 3.0.0, scikit-learn >= 0.24.0, scipy >= 1.6.0, statsmodels >= 0.12.0

Alternatively, create the conda environment: `conda create env -f environment.yml`

## Quick start ##
To run the INSIDER pipeline: `sh INSIDER_Pipeline.sh`

## Usage ##

### Calculate K-mer frequencies

```
python bin/calculate_kmer_frequencies.py \
    split \
    -f test_file.fa \
    -k 2 \
    -n \
    -o test_2mer
```

For each sequence, extract 2-mers and count their frequencies. Ambigous bases (i.e., N's) are ignored)

### Cluster K-mer frequencies

```
python insider_cluster.py \
    consensus \
    --freqDir test_2mer \
    --params params.json \
    -o test_2mer_cIds.txt
```

Cluster sequences based on their K-mers. Hyperparameters can be specified in the JSON file.

### Analyse K-mer frequencies

```
python insider_analyse.py
    main \
    --freqDir test_2mer \
    --cIdFile test_2mer_cIds.txt \
    -o test_2mer_output.txt 
```
Assess the similarity between each cluster and the genome based on their K-mer frequencies.

## Reference ##
For more information, please refer to the following article:

INSIDER: alignment-free detection of foreign DNA sequences

Aidan P. Tay, Brendan Hosking, Cameron Hosking, Denis C. Bauer, and Laurence O.W. Wilson

Computational and Structural Biotechnology Journal, 2021, 19, 3810-3816

