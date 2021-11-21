from Bio import SeqIO, AlignIO

# Custom tree module, see its documentation
from tree_generation import tree_manager

def country_prune(seq_file, keyword, out_file):
    keyword_seq={}
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    for seq_key in seq_dict:
        if keyword in seq_dict[seq_key].description:
            keyword_seq[seq_key] = seq_dict[seq_key]
    SeqIO.write(keyword_seq.values(), out_file, "fasta")
    return keyword_seq

if __name__ == "__main__":
    SEQ_FILE = "data/spike_data_708.fasta"
    COUNTRY = "Australia"
    PRUNED_FILE = "data/spike_"+COUNTRY+".fasta"
    
    country_prune(SEQ_FILE, COUNTRY, PRUNED_FILE)
    # tree_manager(AlignIO.read(PRUNED_FILE.replace(".fasta", "_aligned.fasta"), "fasta"))