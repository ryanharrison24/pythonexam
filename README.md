# Python Kmer Analyzer

## What it does
This script takes a file of DNA sequences and counts how often each kmer 
shows up and what character comes after it. You can set k to whatever length 
you want.

## How to run it
python pythonbadcode.py sequences.txt 2 output.txt

sequences.txt is your input file with one DNA sequence per line
2 is the value of k
output.txt is where the results get saved

## How to run the tests
pytest test_kmer_analyzer.py

## What you need
- Python 3
- pytest (pip install pytest)

## AI use
I used AI to help me understand the bugs in the original code and 
to help me figure out how to write the tests.
