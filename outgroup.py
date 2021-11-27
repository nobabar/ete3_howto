from shutil import copyfile

from Bio import SeqIO, AlignIO

# Custom tree module, see its documentation
from tree_generation import tree_manager


def add_variant(variant_file, align_file, out_file):
    """Add an outgroup to a sequences file.

    Parameters
    ----------
    variant_file : str
        The directory to the file containing the outgroup sequence.
    align_file : str
        The directory to the file containing the sequences.
    out_file : str
        The output file directory.
    """
    # read the outgroup sequence
    variant_seq = SeqIO.read(variant_file, "fasta")

    # copy the content of the sequences file
    copyfile(align_file, out_file)

    # append the outgroup
    with open(out_file, "a") as handle:
        SeqIO.write(variant_seq, handle, "fasta")


if __name__ == "__main__":
    # sequence to add, will constitute the outgroup
    OUTGROUP = "data/outgroup.fasta"
    # the original sequences file
    ALIGNMENT = "data/spike_old_aligned.fasta"
    # the output file
    ALIGNMENT_VARIANT = "data/spike_variant.fasta"

    add_variant(OUTGROUP, ALIGNMENT, ALIGNMENT_VARIANT)

    # TODO: align output file using module of command line
    # generate tree
    compiled = AlignIO.read(ALIGNMENT_VARIANT.replace(".fasta",
                                                      "_aligned.fasta"),
                                      "fasta")
    tree_manager(compiled, aligned=True)
