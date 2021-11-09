"""This module generate trees based on a multi sequences alignment or \
does perform the alignment in case it has not been preiviously done.

See the documentation of the sequences_alignment module for the alignment.

NJ and UPGMA algorithms are both supported for the tree generation.
"""

import os
import matplotlib.pyplot as plt
import Bio.Phylo.TreeConstruction as Tree
from Bio import AlignIO, Phylo

# Custom alignement module, see its documentation
from sequences_alignment import align_seq


def draw_tree(tree, method):
    """Save a phylogenetic tree as a svg file.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        The phylogenetic tree to save.
    method : str
        The method used to generate the tree, used for the file name.
    """
    fig = plt.figure(figsize=(10, 20), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)
    plt.savefig(method+".svg")


def tree_constructor(seq, method="nj", svg=False):
    """Construct a tree using the NJ or UPGMA method.

    Parameters
    ----------
    seq :Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
    method : str (default nj)
        Indicates which method to use, NJ or UPGMA.
    svg : boolean (default False)
        Wether or not the tree should be saved as a svg file.

    Returns
    -------
    Bio.Phylo.BaseTree.Tree
        The generated UPGMA tree
    """
    calculator = Tree.DistanceCalculator('blosum62')
    distance_matrix = calculator.get_distance(seq)
    constructor = Tree.DistanceTreeConstructor()

    if method == "nj":
        tree = constructor.nj(distance_matrix)
    if method == "upgma":
        tree = constructor.upgma(distance_matrix)

    Phylo.write(tree, method+".nwk", "newick")

    if svg:
        draw_tree(tree, method)
    return tree


def parsimony_constructor(tree, seq, svg=False):
    """Construct a tree using the parsimony method.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        A model tree to use as comparison.
    seq : Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
    svg : boolean (default False)
        Wether or not the tree should be saved as a svg file.
    """
    scorer = Tree.ParsimonyScorer()
    searcher = Tree.NNITreeSearcher(scorer)
    constructor = Tree.ParsimonyTreeConstructor(searcher, tree)

    pars_tree = constructor.build_tree(seq)
    Phylo.write(pars_tree, "parsimony.nwk", "newick")

    if svg:
        draw_tree(pars_tree, "parsimony")


def tree_manager(seq, aligned=False, align_method="clustalw", align_out=False,
                 tree_method="nj", svg=False, parsimony=False):
    """Construct a tree from a multi sequence alignement \
    or from a multi sequence file that will get aligned.

    Parameters
    ----------
    seq : Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
    aligned : boolean (default False)
        Wether or not the seq specified is already aligned.
    align_method : str (default clustalw)
        The alignment algorithm, either clustalw or muscle.
    align_out : boolean (default False)
        Wether or not an output file should be generated from the alignment.
    tree_method : str (default nj)
        The method to use for the tree generation, either NJ or UPGMA.
    svg : boolean (default False)
        Wether or not the tree should be saved as a svg file.
    parsimony : boolean (default False)
        Wether or not a parsimony tree should be generated.
    """
    if isinstance(seq, str):
        if os.path.isfile(seq):
            if aligned:
                sequence = AlignIO.read(seq, "fasta")
            else:
                sequence = align_seq(seq, align_method, align_out)
        else:
            print("Error: Not a valid file.")
    elif isinstance(seq, AlignIO.MultipleSeqAlignment):
        sequence = seq
    else:
        print("Error: Input sequence need to be a fasta file \
              or an multi sequences alignment.")

    if tree_method in ('nj', 'upgma'):
        tree = tree_constructor(sequence, tree_method, svg)
    else:
        print("Error: method should either be 'nj' or 'upgma'. \
              Leave blank for nj method.")

    if parsimony:
        parsimony_constructor(tree, svg)


if __name__ == "__main__":
    SEQ_FILE = "data/spike_data_708.fasta"
    aligned_seq = AlignIO.read(SEQ_FILE.replace(".fasta", "_aligned.fasta"),
                               "fasta")
    tree_manager(SEQ_FILE, aligned=True, svg=True)
