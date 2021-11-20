from Bio import SeqIO

def country_prune(seq_file, keyword):
    keyword_seq={}
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    for seq_key in seq_dict:
        if keyword in seq_dict[seq_key].description:
            keyword_seq[seq_key] = seq_dict[seq_key]
    return keyword_seq

if __name__ == "__main__":
    SEQ_FILE = "data/spike_data_708.fasta"
    COUNTRY = "Australia"
    pruned_sequences = country_prune(SEQ_FILE, COUNTRY)
    SeqIO.write(pruned_sequences.values(), "data/spike_"+COUNTRY+".fasta", "fasta")



