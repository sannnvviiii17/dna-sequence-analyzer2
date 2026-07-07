# dna-sequence-analyzer2
A command-line Python tool that analyzes DNA sequences from FASTA files — built as a learning project combining Linux, Python, and basic bioinformatics.

## Features

- Parses one or more sequences from a FASTA file
- Reports sequence length and base composition (A/T/G/C counts)
- Calculates GC content
- Computes the reverse complement strand
- Transcribes DNA into mRNA
- Translates DNA into protein (stops at the first stop codon)
- Searches for a motif/subsequence and reports its position(s)

## Requirements

- Python 3.7+ (no external libraries needed)

## Usage

Basic analysis:
python3 dna_analyzer.py --file example.fasta

Include protein translation:
python3 dna_analyzer.py --file example.fasta --translate

Search for a motif:
python3 dna_analyzer.py --file example.fasta --motif ATG
Combine flags:
python3 dna_analyzer.py --file example.fasta --motif ATG --translate
 What I learned

This project helped me practice: reading files in Python, dictionaries, string slicing, loops with custom step sizes, command-line argument parsing with `argparse`, and basic Linux/git workflow.


