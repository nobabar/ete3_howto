from Bio import SeqIO, AlignIO
from datetime import date
import pandas as pd

# Custom tree module, see its documentation
from tree_generation import tree_manager

def country_prune(seq_file, country, out_file):
    keyword_seq={}
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    for seq_key in seq_dict:
        if country in seq_dict[seq_key].description:
            keyword_seq[seq_key] = seq_dict[seq_key]
    SeqIO.write(keyword_seq.values(), out_file, "fasta")
    return keyword_seq


def date_prune(seq_file, out_file):
    keyword_seq={}
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    
    for seq_key in seq_dict:
        seq_date = seq_dict[seq_key].description.rsplit("|", 1)[1].strip()
        seq_period = pd.Period(seq_date)
        if seq_period.year == 2020 and seq_period.month < 4:
            keyword_seq[seq_key] = seq_dict[seq_key]
    SeqIO.write(keyword_seq.values(), out_file, "fasta")
    return keyword_seq


if __name__ == "__main__":
    SEQ_FILE = "data/spike_data_708.fasta"
    
    COUNTRY = "Australia"
    PRUNED_FILE = "data/spike_"+COUNTRY+".fasta"
    # country_prune(SEQ_FILE, COUNTRY, PRUNED_FILE)
    # tree_manager(AlignIO.read(PRUNED_FILE.replace(".fasta", "_aligned.fasta"), "fasta"))
    
    PRUNED_OLD = "data/spike_old.fasta"
    date_prune(SEQ_FILE, PRUNED_OLD)
