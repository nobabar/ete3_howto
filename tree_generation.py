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


def upgma_constructor(seq, svg):
    """Construct a tree using the UPGMA method.
    
    Parameters
    ----------
    seq :Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
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
    upgma_tree = constructor.upgma(distance_matrix)
    
    Phylo.write(upgma_tree, seq_file.replace(".fasta", "_upgma.nwk"), "newick")
    
    if svg:
        draw_tree(upgma_tree, "upgma")
    return upgma_tree

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
    Phylo.write(pars_tree, seq_file.replace(".fasta", "_parsimony.nwk"), "newick")
    
    if svg:
            draw_tree(pars_tree, "parsimony")


def tree_constructor(seq, aligned=False, svg=False, parsimony=False):
    """Construct a tree using the UPGMA method.
    
    Parameters
    ----------
    seq : Bio.Align.MultipleSeqAlignment
        A multiple sequences alignement.
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
    
    upgma_tree = upgma_constructor(sequence, svg)

    if parsimony:
        parsimony_constructor(upgma_tree, seq, svg)

if __name__ == "__main__":
    seq_file = "og_variant.fasta"
    #aligned_seq = AlignIO.read(seq_file.replace(".fasta", "_aligned.fasta"), "fasta")
    tree_constructor(seq_file)
