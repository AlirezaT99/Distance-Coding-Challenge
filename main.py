# Author: Alireza Tajmirriahi
# Date: 2024-12-04
# GitHub: https://github.com/AlirezaT99/Distance-Coding-Challenge

import os
import re
import argparse
import itertools
import numpy as np
import pandas as pd

import mmh3  # For Murmur hashing
from Bio import SeqIO, Phylo  # For reading FASTA files and drawing trees
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from tqdm import tqdm

tqdm.pandas()  # Enables progress_apply

# Constants
FILE_PATHS = [
    "R6.fa",
    "TIGR4.fa",
    "14412_3#82.contigs_velvet.fa",
    "14412_3#84.contigs_velvet.fa"
]
DEBUG = True
DEFAULT_KMER_SIZE = 14
DEFAULT_SKETCH_SIZE = 1000
COMPLEMENTS = {"A": "T", "T": "A", "C": "G", "G": "C"}


class Logger:
    """
    A logger class that writes messages to a file
    and optionally prints them to the console.
    """

    def __init__(self, file_path, output="", debug=False):
        self.file_path = file_path
        self.output = output
        self.debug = debug

    def save_output(self) -> None:
        with open(self.file_path, "w") as f:
            f.write(self.output)

    def dprint(self, msg) -> None:
        self.output += msg + "\n"
        if self.debug:
            print(msg)


def extract_kmers(sequences: list, kmer_size: int) -> pd.Series:
    """
    Returns a list of all valid kmers of size kmer_size in
    the given sequence, after removing ambiguous base chunks.
    """
    sequences = filter(
        lambda fragment: len(fragment) >= kmer_size,
        [fragment for sequence in sequences for fragment in re.split("N+", str(sequence.seq))]
    )  # Naive way of cutting out ambiguous base chunks
    result = []
    for sequence in sequences:
        result.extend([sequence[i:i + kmer_size] for i in tqdm(range(len(sequence) - kmer_size + 1))])
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

    rev_kmers = kmers.progress_apply(reverse_complement)
    kmers_hash, rev_kmers_hash = kmers.progress_apply(hash_kmer), rev_kmers.progress_apply(hash_kmer)
    return np.min([kmers_hash, rev_kmers_hash], axis=0)


def test_hash_consistency(sequences: pd.Series, hash_values, attempts=1) -> None:
    """
    Tests the consistency of the hash function by comparing the hash values
    of the given sequences to the hash values of the same sequences.
    """
    for _ in range(attempts):
        attempt_values = calculate_hash(sequences)
        assert (hash_values == attempt_values).all(), "Hash values are inconsistent."
    dprint("\tHash values are consistent.")


def make_sketch(kmers: pd.Series, sketch_size: int) -> np.array:
    """
    Returns a MinHash sketch of the given kmers with the given sketch size.
    """
    hash_values = calculate_hash(kmers.drop_duplicates())
    return np.sort(hash_values)[:sketch_size]


def jaccard_distance(kmers1, kmers2) -> float:
    """
    Returns the Jaccard distance between the kmer profiles of two sequences.
    The presence of a kmer/hash value is considered rather than the count.
    """
    return 1 - len(set(kmers1) & set(kmers2)) / len(set(kmers1) | set(kmers2))


def evaluate_neighbour_joining(dist_matrix: np.array, names: list, output_path: str, msg: str) -> None:
    """
    Constructs a phylogenetic tree from the given distance matrix and names,
    and saves it to the given output path followed by the given message.
    """
    dist_matrix = [
        [0],
        [dist_matrix[0, 1], 0],
        [dist_matrix[0, 2], dist_matrix[1, 2], 0],
        [dist_matrix[0, 3], dist_matrix[1, 3], dist_matrix[2, 3], 0]
    ]
    dist_matrix = DistanceMatrix(names=names, matrix=dist_matrix)
    tree = DistanceTreeConstructor().nj(dist_matrix)

    with open(output_path, "a") as f:
        f.write(f"{msg}\n")
        Phylo.draw_ascii(tree, file=f)
        f.write("\n\n")


def read_fasta(file_path) -> list:
    """
    Reads the FASTA file at the given path and returns a SeqRecord iterator.
    """
    with open(file_path, "r") as f:
        return list(SeqIO.parse(f, "fasta"))


def main(args) -> None:
    dprint(f"Analysis started. K={args.kmer_size}, Sketch Size={args.sketch_size}")
    dprint("Step 1 - Reading FASTA files and calculating kmer counts")
    kmer_counts = dict()  # Storing the kmer counts dict for each file
    kmers_dict = dict()
    for file_path in FILE_PATHS:
        dprint(f"\tReading {file_path}...")
        sequences = read_fasta(f"{args.data_dir}/{file_path}")

        dprint(f"\tExtracting kmers from {file_path}...")
        kmers = extract_kmers(sequences, args.kmer_size)
        kmers_count = kmers.value_counts().to_dict()
        kmers_dict[file_path] = kmers
        kmer_counts[file_path] = kmers_count

    dprint("Step 2 - Calculating Jaccard distance between the two pairs")
    full_distance_matrix = np.zeros((4, 4))
    for pair in itertools.combinations(range(4), 2):  # Indices in the paths list   
        distance = jaccard_distance(kmer_counts[FILE_PATHS[pair[0]]].keys(), kmer_counts[FILE_PATHS[pair[1]]].keys())
        full_distance_matrix[pair[0], pair[1]] = distance
    for pair in [(0, 1), (2, 3)]:  # References and drafts
        dprint(f"\tJaccard distance between {FILE_PATHS[pair[0]]} and {FILE_PATHS[pair[1]]}: "
               + f"{full_distance_matrix[pair[0], pair[1]]}\n")

    dprint("Step 3 - Calculating MinHash signatures")
    for file_path in FILE_PATHS:
        dprint(f"\tCalculating Hash for {file_path} kmers (three steps)...")
        hash_values = calculate_hash(kmers_dict[file_path])
        dprint(f"\tTesting Hash consistency...")
        test_hash_consistency(kmers_dict[file_path], hash_values)

    dprint("Step 4 - Creating a Sketch from input sequences")
    for file_path in FILE_PATHS:
        dprint(f"\tMaking sketch for {file_path} with size {args.sketch_size}")
        sketch = make_sketch(kmers_dict[file_path], args.sketch_size)
        output_sketch_path = f"{args.output_dir}/{file_path}_sketch_{args.sketch_size}.npy"
        np.save(output_sketch_path, sketch)
        dprint(f"\tSaved sketch to {output_sketch_path}.")

    dprint("Step 5 - Calculating Jaccard distance between the sketches")
    jaccard_distance_matrix = np.zeros((4, 4))
    for pair in itertools.combinations(range(4), 2):  # Indices in the paths list
        sketch1 = np.load(f"{args.output_dir}/{FILE_PATHS[pair[0]]}_sketch_{args.sketch_size}.npy")
        sketch2 = np.load(f"{args.output_dir}/{FILE_PATHS[pair[1]]}_sketch_{args.sketch_size}.npy")
        distance = jaccard_distance(sketch1, sketch2)
        jaccard_distance_matrix[pair[0], pair[1]] = distance
        dprint(f"\tJaccard distance between {FILE_PATHS[pair[0]]} and {FILE_PATHS[pair[1]]} sketches: {distance}")

    dprint("(Bonus) Step 6 - Constructing a phylogenetic tree based on the distances")
    output_tree_path = f"{args.output_dir}/NJ_tree_k{args.kmer_size}_sketch{args.sketch_size}.txt"
    evaluate_neighbour_joining(jaccard_distance_matrix, FILE_PATHS, output_tree_path, "Sketch Distance")
    evaluate_neighbour_joining(full_distance_matrix, FILE_PATHS, output_tree_path, "Full Distance")
    dprint(f"\tPhylogenetic tree saved to {output_tree_path}.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Distance Coding Challenge - Alireza Tajmirriahi")
    parser.add_argument("--data-dir", "-d", type=str, help="Data directory path", default="s_pneumoniae_genomes")
    parser.add_argument("--kmer-size", "-k", type=int, help="K-mer size", default=DEFAULT_KMER_SIZE)
    parser.add_argument("--sketch-size", "-s", type=int, help="Sketch size", default=DEFAULT_SKETCH_SIZE)
    parser.add_argument("--output-dir", "-o", type=str, help="Output directory path", default="output")
    parser.add_argument("--debug", "-v", type=bool, help="Enable verbose output", default=DEBUG)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    output_log_path = os.path.join(args.output_dir, f"log-{pd.Timestamp.now().strftime('%H-%M-%S')}.txt")
    logger = Logger(output_log_path, debug=args.debug)
    dprint = logger.dprint
    try:
        main(args)
    except Exception as e:
        dprint(f"An error occurred: {e}")
    finally:
        logger.save_output()
        print(f"Log written to {output_log_path}.")
