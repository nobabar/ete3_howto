from io import StringIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline


def align_seq(in_dir, method="clustalw", out_file=False):
    """Align multiple sequences using the MUSCLE algorithm.
    
    Parameters
    ----------
    in_dir : dir
        The directory to the file containing the sequences.
    out_dir : boolean
        Wether or not a output file should be generated.
    
    Returns
    -------
    align
        The aligned sequences.
    """
    clustalo_exe = ".\clustal-omega\clustalo.exe"
    # muscle_exe = ".\muscle3.8.31_i86win32.exe"
    
    if out_file:
        out_dir = in_dir.replace(".fasta", "_aligned.fasta")
        open(out_dir, "w+").close()
        clustalw_cline = ClustalOmegaCommandline(clustalo_exe,
                                                 infile=in_dir,
                                                 outfile=out_dir,
                                                 force=True)
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_dir, "fasta")

    else:
        clustalw_cline = ClustalOmegaCommandline(clustalo_exe,
                                                 infile=in_dir,
                                                 auto=True)
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(StringIO(stdout), "fasta")
    return align

if __name__ == "__main__":
    seq_file = "data/spike_data_708.fasta"
    sequence = align_seq(seq_file, "clustalw", True)