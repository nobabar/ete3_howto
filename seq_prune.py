from Bio import SeqIO, AlignIO
import pandas as pd

# Custom tree module, see its documentation
from tree_generation import tree_manager

def country_prune(seq_file, country, out_file):
    """Prune sequence given a country name.

    Parameters
    ----------
    seq_file : str
        The directory to the file containing the sequences.
    country : str
        The country we want the sequences from.
    out_file : str
        The output file directory.

    Returns
    -------
    dict
        The pruned dictionnary of sequences.
    """
    keyword_seq={}

    # read the sequences and store them as a dictionnary
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

    for seq_key in seq_dict:
        # keep the sequence if its country match the wanted country
        if country in seq_dict[seq_key].description:
            keyword_seq[seq_key] = seq_dict[seq_key]
    # write the pruned dictionnary in the output file
    SeqIO.write(keyword_seq.values(), out_file, "fasta")
    return keyword_seq


def date_prune(seq_file, date, out_file):
    """Prune sequence given a date.

    Parameters
    ----------
    seq_file : str
        The directory to the file containing the sequences.
    date : list
        A year and month, all prior sequences will be kept. (e.g. [2020, 4])
    out_file : str
        The output file directory.

    Returns
    -------
    dict
        The pruned dictionnary of sequences.
    """
    keyword_seq={}

    # read the sequences and store them as a dictionnary
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

    for seq_key in seq_dict:
        # parse the sequence date
        seq_date = seq_dict[seq_key].description.rsplit("|", 1)[1].strip()
        seq_period = pd.Period(seq_date)
        # keep the sequence if its country match the wanted country
        if seq_period.year < date[0] and seq_period.month < date[1]:
            keyword_seq[seq_key] = seq_dict[seq_key]
    # write the pruned dictionnary in the output file
    SeqIO.write(keyword_seq.values(), out_file, "fasta")
    return keyword_seq


if __name__ == "__main__":
    # main sequence file
    SEQ_FILE = "data/spike_data_708.fasta"

    # keep only sequences from Australia
    COUNTRY = "Australia"
    # output file
    PRUNED_FILE = "data/spike_"+COUNTRY+".fasta"
    country_prune(SEQ_FILE, COUNTRY, PRUNED_FILE)
    # align output file using module of command line
    # generate tree
    tree_manager(AlignIO.read(PRUNED_FILE.replace(".fasta", "_aligned.fasta"), "fasta"))

    # output file
    PRUNED_OLD = "data/spike_old.fasta"
    date_prune(SEQ_FILE, [2020, 4], PRUNED_OLD)
    # TODO: align output file using module of command line
    # generate tree
    tree_manager(AlignIO.read(PRUNED_OLD.replace(".fasta", "_aligned.fasta"), "fasta"))
