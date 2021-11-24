# Custom modules, see their documentation
from sequences_alignment import align_seq
from tree_generation import tree_manager

from Bio import SeqIO, AlignIO
from shutil import copyfile


def add_variant(variant_file, align_file, out_file):
    variant_seq = SeqIO.read(variant_file, "fasta")
    
    copyfile(align_file, out_file)
    
    with open(out_file, "a") as handle:
        SeqIO.write(variant_seq, handle, "fasta")





if __name__ == "__main__":
    OUTGROUP = "data/outgroup.fasta"
    ALIGNMENT = "data/spike_old_aligned.fasta"
    ALIGNMENT_VARIANT = "data/spike_variant.fasta"
    
    add_variant(OUTGROUP, ALIGNMENT, ALIGNMENT_VARIANT)
    # compiled_alignment = align_seq(ALIGNMENT_VARIANT, out_file=True)
    compiled_alignment = AlignIO.read(ALIGNMENT_VARIANT.replace(".fasta", "_aligned.fasta"), "fasta")
    tree_manager(compiled_alignment, aligned=True,
                 tree_method="nj", svg=False, parsimony=False)