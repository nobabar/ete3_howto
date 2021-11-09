import os
import matplotlib.pyplot as plt
import Bio.Phylo.TreeConstruction as Tree
from Bio import AlignIO, Phylo

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
    plt.savefig(seq_file.replace(".fasta", "_"+method+".svg"))


def tree_constructor(seq, method, svg):
    """Construct a tree using the NJ or UPGMA method.
    
    Parameters
    ----------
    seq :Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
    method : str
        Indicates which method to use, NJ or UPGMA.
    svg : boolean
        Wether or not the tree should be saved as a svg file.
    
    Returns
    -------
    Bio.Phylo.BaseTree.Tree
        The generated UPGMA tree
    """
    calculator = Tree.DistanceCalculator('blosum62')
    distance_matrix = calculator.get_distance(seq)
    constructor = Tree.DistanceTreeConstructor()
    
    if method == "NJ":
        tree = constructor.nj(distance_matrix)
    if method == "UPGMA":
        tree = constructor.upgma(distance_matrix)
    
    Phylo.write(tree, seq_file.replace(".fasta", "_"+method+".nwk"), "newick")
    
    if svg:
        draw_tree(tree, method)
    return tree


def parsimony_constructor(tree, seq, svg):
    """Construct a tree using the parsimony method.
    
    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        A model tree to use as comparison.
    seq : Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
    svg : boolean
        Wether or not the tree should be saved as a svg file. 
    """
    scorer = Tree.ParsimonyScorer()
    searcher = Tree.NNITreeSearcher(scorer)
    constructor = Tree.ParsimonyTreeConstructor(searcher, tree)
    
    pars_tree = constructor.build_tree(seq)
    Phylo.write(pars_tree,
                seq_file.replace(".fasta", "_parsimony.nwk"),
                "newick")
    
    if svg:
            draw_tree(pars_tree, "parsimony")


def tree_manager(seq, aligned=False, method="nj", svg=False, parsimony=False):
    """Construct a tree from a multi sequence alignement
    or from a multi sequence file that will get aligned.
    
    Parameters
    ----------
    seq : Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
    aligned : boolean
        Wether or not the seq specified is already aligned.
    method : str
        The method to use for the tree generation, either NJ or UPGMA.
    svg : boolean
        Wether or not the tree should be saved as a svg file. 
    parsimony : boolean
        Wether or not a parsimony tree should be generated. 
    """
    if isinstance(seq, str):
        if os.path.isfile(seq):
            if aligned:
                sequence = AlignIO.read(seq)
            else:
                sequence = align_seq(seq)
        else :
            print("Error: Not a valid file.")
    elif isinstance(seq, AlignIO.MultipleSeqAlignment):
        sequence = seq
    else:
        print("Error: Input sequence need to be a fasta file \
              or an multi sequences alignment.")
    
    if method == "nj" or method == "upgma":   
        tree = tree_constructor(sequence, method, svg)
    else:
        print("Error: method should either be 'nj' or 'upgma'. \
              Leave blank for nj method.")

    if parsimony:
        parsimony_constructor(tree, seq, svg)

if __name__ == "__main__":
    seq_file = "data/spike_data_708.fasta"
    aligned_seq = AlignIO.read(seq_file.replace(".fasta", "_aligned.fasta"),
                               "fasta")
    tree_manager(aligned_seq, True, svg=True)
