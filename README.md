# Distance Coding Challenge
My attempt at the Distance Coding Challenge.

## How to Run
This code uses `mmh3` (for Murmur hashing) and `Biopython` (to read FASTA files and draw phylogenetic trees). You can install all dependencies by running:

```bash
pip install -r requirements.txt
```

To run the analysis (located in [`main.py`](main.py)), use the following command:

```bash
# `python main.py` if you use the default parameters, otherwise:
python main.py --data-dir "s_pneumoniae_genomes" -k 14 -s 1000
```

## Observations
### Comparing Isolates
- Full distances provide the most accurate estimation of Jaccard similarity by considering the entire k-mer profile; whereas MinHash distances approximate this similarity while decreasing computational and memory cost.
- As demonstrated in a [sample log](./output/log-00-18-03.txt) (Steps 2 & 5), introducing sketches significantly decreases the estimated distance between the reference pair, correcting the unrealistically high distance observed from the full profiles. Furthermore, using sketches enhances the separation between the reference and draft pairs, making the distances more distinguishable (see the next section for details).
- Such estimation methods are particularly useful in scenarios with limited computational resources or for applications where adjustable estimations suffice.
- Enlarging the sketch size will result in a more stable estimation of the Jaccard full distances.
- Increasing the sketch size results in more stable (and more costly) approximations of the full Jaccard distances.
### Neighbour-Joining Tree
- The neighbor-joining (NJ) tree construction, implemented using Biopython, is included as the [6th step](main.py#L191-L195) in the code. Trees are generated using both full Jaccard distances and MinHash distances, and the outputs are visualized in the same file (sample output [here](./output/NJ_tree_k14_sketch1000.txt)).
- Based on the sample results, the MinHash distances provide a more reasonable separation between the reference and draft pairs. In other words, Distances within each reference or draft pair are reduced, reflecting closer relationships, while distances between the reference and draft pairs are increased, indicating clearer separation.
- This adjustment aligns better with biological expectations, justifying our choice of approach and input parameters.
