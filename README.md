# Kmer Analyzer

# What This Does
This script reads DNA sequences from a file and counts how often each
k-mer (substring of length k) appears, and what character comes after it.

# How to Run It
python kmer_analyzer.py sequences.txt 2 output.txt

- sequences.txt = your input file with one DNA sequence per line
- 2 = the value of k (length of substrings)
- output.txt = where results will be saved

# How to Run the Tests
pytest test_kmer_analyzer.py

# Requirements
- Python 3.x
- pytest (pip install pytest)
