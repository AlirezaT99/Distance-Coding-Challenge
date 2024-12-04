# Author: Alireza Tajmirriahi
# Date: 2024-12-04
# GitHub: https://github.com/AlirezaT99/Distance-Coding-Challenge

import os
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

from Bio import SeqIO  # For reading FASTA files
import mmh3  # For Murmur hashing

# Constants
FILE_PATHS = [
    "R6.fa",
    "TIGR4.fa",
    "14412_3#82.contigs_velvet.fa",
    "14412_3#84.contigs_velvet.fa"
]
DEFAULT_KMER_SIZE = 14
DEFAULT_SKETCH_SIZE = 1000
COMPLEMENTS = {"A": "T", "T": "A", "C": "G", "G": "C"}

DEBUG = True
output = ""

def main(args) -> None:    
    # TODO
    pass

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Distance Coding Challenge - Alireza Tajmirriahi")
    parser.add_argument("--seq-pair", "-p", type=str, help="Sequence pairs to analyze", choices=["references", "draft"], default="references")
    parser.add_argument("--data-dir", "-d", type=str, help="Data directory path", default="s_pneumoniae_genomes")
    parser.add_argument("--kmer-size", "-k", type=int, help="K-mer size", default=DEFAULT_KMER_SIZE)
    parser.add_argument("--sketch-size", "-s", type=int, help="Sketch size", default=DEFAULT_SKETCH_SIZE)
    parser.add_argument("--output-dir", "-o", type=str, help="Output directory path", default="output")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    main(args)
