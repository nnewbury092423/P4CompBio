[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/iI4DdsQR)
# ECE 208 Homework 4: Tree-based Multiple Sequence Alignment 

You have learned in the lecture that the complexity of dynamic programming (DP) for pairwise and threeway sequence alignment are O(n<sup>2</sup>) and O(n<sup>3</sup>), respectively, where n is maximum length of the input sequences. 
Generalizing this approach, the cost to align k sequences is O(n<sup>k</sup>), which means the algorithm's complexity grows exponentially with the numer of species. 
Despite the intensive computational cost, multiple sequence alignment (MSA) has a large number of applications in bioinformatics. 
Therefore, MSA is always an active research in the field. 

In this part of the assignment, you will design an algorithm to align multiple sequences (e.g. up to 50 species) given a *guidance tree*. The guidance tree gives information about the evolutionary relationships of the sequences: its topology shows how the species related to each other and the length of each branch shows the expected number of substitutions occured on that branch. 
You will use the guidance tree to come up with an algorithm that can infer an alignment more accurately and more efficiently than the DP algorithm. 

As in homework 2, you are given the BLOSUM62 score matrix, the homologous probabilies, and the indel-rate. In addition, you are also given the guidance tree (i.e. a binary rooted tree that has non-negative branch lengths).
You will implement your algorithm using Python3 and analyze the complexity and accuracy of your algorithm on a simulated dataset. 

### Input:
A rooted binary tree, a set of taxon - sequence pairs, where each taxon is mapped exactly to one sequence, and an indel rate.

### Output:
MSA of the input sequences.

### Starting code
* **[compute_MSA.py](compute_MSA.py):** Given to you. This script handles input/output. Please **DO NOT** modify it. 

```
python3 compute_MSA.py -i [INPUT_FILE.fas] -t [TREE_FILE.newick] -o [OUTPUT_FILE.fas] -r [INDEL_RATE]
```

* **[msa/todo.py](msa/todo.py):** You will complete the TODO in function ```MSA```.
* **[autocheck.py](autocheck.py):** Given to you. Please **DO NOT** modify. As always, this script automatically runs some basic tests for your program. To run autocheck, you will need to install `func-timeout`:

```
python3 -m pip install func-timeout
```

then simply type

```
python3 autocheck.py
```

## Deliverables:
* [msa/todo.py](msa/todo.py)
* **writeup_<YOUR_PID>.pdf**: description, proof of correctness, and complexity analyses of the algorithms in part 1 and part 2.

## Grade Breakdown (100 Points)
   * Accuracy (70 Points)
     * 5 tests for low indel rate (0.01), 7 pts each. Full grade: SP-Score 0.95 or above. Half grade: SP-Score between 0.8 and 0.95
     * 5 tests for high indel rate (0.1), 7 pts each. Full grade: SP-Score 0.85 or above. Half grade: SP-Score between 0.8 and 0.85
     * Time limit per test: 600 seconds
   * Algorithm description and analysis (30 Points): describe your algorithm to compute the MSA and analyze the time complexity with respect to n and k, where n is the (largest) sequence length and k is the number of species.
   
## Restrictions on the usage of external software/library
* You can use any of the built-in functions in TreeSwift. Other than that, you are NOT allowed to use any other external libaries other than the built-in of Python 3.
