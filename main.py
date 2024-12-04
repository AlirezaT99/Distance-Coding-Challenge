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


def dprint(msg) -> None:
    output += msg + "\n"
    if DEBUG:
        print(msg)

def get_kmers(sequences: list, kmer_size: int) -> pd.Series:
    """
    Returns a list of all kmers of size kmer_size in the given sequence.
    """
    result = []
    for sequence in sequences:
        result.extend([sequence.seq[i:i + kmer_size] for i in tqdm(range(len(sequence) - kmer_size + 1))])
    return pd.Series(result)

def calculate_hash(kmers: pd.Series) -> np.array:
    """
    Returns the Murmur hash of the given kmers by assigning the minimum value
    between the kmer and its reverse complement as the canonical kmer hash.
    """
    def reverse_complement(sequence: str) -> str:
        """
        Returns the reverse complement of the given sequence.
        """
        return "".join(COMPLEMENTS[base] for base in sequence[::-1])
    
    def hash_kmer(kmer: str) -> int:
        return mmh3.hash(kmer, seed=0, signed=False)

    rev_kmers = kmers.apply(reverse_complement)
    kmers_hash, rev_kmers_hash = kmers.apply(hash_kmer), rev_kmers.apply(hash_kmer)
    return np.min([kmers_hash, rev_kmers_hash], axis=0)

def test_hash_consistency(sequences: pd.Series, hash_values, attempts=1) -> None:
    """
    Tests the consistency of the hash function by comparing the hash values
    of the given sequences to the hash values of the same sequences.
    """
    for _ in range(attempts):
        attempt_values = calculate_hash(sequences)
        assert (hash_values == attempt_values).all(), "Hash values are inconsistent."

def make_sketch(kmers: set, sketch_size: int) -> np.array:
    """
    Returns a MinHash sketch of the given kmers with the given sketch size.
    """
    hash_values = calculate_hash(kmers)
    return np.sort(hash_values)[:sketch_size]

def jaccard_distance(kmers1, kmers2) -> float:
    """
    Returns the Jaccard distance between the kmer profiles of two sequences.
    The presence of a kmer/hash value is considered rather than the count.
    """
    return 1 - len(set(kmers1) & set(kmers2)) / len(set(kmers1) | set(kmers2))

def read_fasta(file_path) -> list:
    """
    Reads the FASTA file at the given path and returns a SeqRecord iterator.
    """
    with open(file_path, "r") as file:
        return list(SeqIO.parse(file, "fasta"))

def main(args) -> None:    
    file_pair = (0, 1) if args.seq_pair == "references" else (2, 3)
    
    dprint("Step 1 - Read FASTA files and calculate kmer counts")
    files_kmer_counts = dict()  # Storing the kmer counts dict for each file
    file_names = FILE_PATHS[file_pair[0]], FILE_PATHS[file_pair[1]]
    for file_path in FILE_PATHS:
        sequences = read_fasta(f"{args.data_dir}/{file_path}")
        kmers = get_kmers(sequences, args.kmer_size)
        kmers_count = kmers.value_counts().to_dict()
        files_kmer_counts[file_path] = kmers_count
    
    # Step 2 - Calculate Jaccard distance between the two files
    distance = jaccard_distance(files_kmer_counts[file_names[0]].keys(), files_kmer_counts[file_names[1]].keys())
    output = f"Jaccard distance between {file_names[0]} and {file_names[1]}: {distance}"
    
    # Step 3 - Calculate MinHash signatures
    for file_path in FILE_PATHS:
        hash_values = calculate_hash(kmers)
        test_hash_consistency(kmers, hash_values)
    
    # Step 4 - Create a Sketch from input sequences
    for file_path in FILE_PATHS:
        sketch = make_sketch(files_kmer_counts[file_path].keys(), args.sketch_size)
        np.save(f"{args.output_dir}/{file_path}_sketch_{args.sketch_size}.npy", sketch)
    
    # Step 5 - Calculate Jaccard distance between the two sketches
    for pair in [(0, 1), (2, 3)]:
        sketch1 = np.load(f"{args.output_dir}/{FILE_PATHS[pair[0]]}_sketch_{args.sketch_size}.npy")
        sketch2 = np.load(f"{args.output_dir}/{FILE_PATHS[pair[1]]}_sketch_{args.sketch_size}.npy")
        distance = jaccard_distance(sketch1, sketch2)
        dprint(f"\nJaccard distance between {FILE_PATHS[pair[0]]} and {FILE_PATHS[pair[1]]} sketches: {distance}")
    # TODO Comparing Isolates and Neighbor Joining
    

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
