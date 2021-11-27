"""This module generate trees based on a multi sequences alignment or \
does perform the alignment in case it has not been previously done.

See the documentation of the sequences_alignment module for the alignment.

NJ and UPGMA algorithms are both supported for the tree generation.
"""

import os

from Bio import AlignIO, Phylo
import Bio.Phylo.TreeConstruction as Tree

import matplotlib.pyplot as plt
import pandas as pd

# Custom alignement module, see its documentation
from sequences_alignment import align_seq


def draw_tree(tree):
    """Save a phylogenetic tree as a svg file.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        The phylogenetic tree to save.
    """
    # draw a matplotlib figure containing the tree
    fig = plt.figure(figsize=(10, 20), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)
    
    # save tree with current timestamp as name
    now = str(pd.to_datetime("now")).replace(" ", "_").rsplit(".")[0]
    plt.savefig(now+".svg")


def tree_constructor(seq, spike708=False, svg=False):
    """Construct a tree using the NJ or UPGMA method.

    Parameters
    ----------
    seq : Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
    spike708 : boolean (default False)
        If the data provided are the original spike data.
    svg : boolean (default False)
        Wether or not the tree should be saved as a svg file.
    """
    calculator = Tree.DistanceCalculator('blosum90')

    # if the sequences provided are the original ones
    # and the distance matrix exist
    if spike708 and os.path.isfile("data/distance_matrix.csv"):
        dm_df = pd.read_csv("data/distance_matrix.csv",
                            index_col=0)

        # keep only the lower part of the matrix
        lower_triangle = []
        k = 1
        for i in dm_df.values:
            lower_triangle.append(list(map(float, i[:k])))
            k += 1

        # construct a DistanceMatrix object from the dataframe
        distance_matrix = Tree._DistanceMatrix(names = list(dm_df.index),
                                               matrix=lower_triangle)

    else:
        distance_matrix = calculator.get_distance(seq)

        # save the distance matrix
        # distances=[]
        # for dis_lis in distance_matrix:
        #     distances.append(dis_lis)
        #     distancedf = pd.DataFrame(distances, index=distance_matrix.names)

        #     # save the distance matrix with current timestamp as name
        #     now = str(pd.to_datetime("now")).replace(" ", "_").rsplit(".")[0]
        #     distancedf.to_csv("data/dm_"+now+".csv")

    constructor = Tree.DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # save tree with current timestamp as name
    now = str(pd.to_datetime("now")).replace(" ", "_").rsplit(".")[0]
    Phylo.write(tree, now+"_tree.nwk", "newick")

    if svg:
        draw_tree(tree)


def tree_manager(seq, aligned=False, align_method="clustalw",
                 align_out=False, svg=False):
    """Construct a tree from a multi sequence alignement \
    or from a multi sequence file that will get aligned.

    Parameters
    ----------
    seq : Bio.Align.MultipleSeqAlignment or str
        A multiple sequences alignement or the path the a fasta file.
    aligned : boolean (default False)
        Wether or not the seq specified is already aligned.
    align_method : str (default clustalw)
        The alignment algorithm, either clustalw or muscle.
    align_out : boolean (default False)
        Wether or not an output file should be generated from the alignment.
    svg : boolean (default False)
        Wether or not the tree should be saved as a svg file.
    """
    # if the sequence provided is a path test if it exist
    if isinstance(seq, str):
        if os.path.isfile(seq):
            if aligned:
                sequence = AlignIO.read(seq, "fasta")
                og = True if seq == "data/spike_data_708_aligned.fasta" \
                    else False
            else:
                sequence = align_seq(seq, align_method, align_out)
                og = True if seq == "data/spike_data_708.fasta" \
                    else False
        else:
            print("Error: Not a valid file.")

    # else test if it is a multiple sequences alignment
    elif isinstance(seq, AlignIO.MultipleSeqAlignment):
        sequence = seq

        og = True if seq == AlignIO.read("data/spike_data_708_aligned.fasta") \
            else False
    else:
        print("Error: Input sequence need to be a fasta file \
              or an multi sequences alignment.")
        return

    tree_constructor(sequence, og, svg)


if __name__ == "__main__":
    SEQ_FILE = "data/spike_data_708.fasta"
    aligned_seq = AlignIO.read(SEQ_FILE.replace(".fasta", "_aligned.fasta"),
                               "fasta")
    tree_manager(aligned_seq, aligned=True, svg=True)
