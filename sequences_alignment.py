"""This module align a multi sequences fasta file.

ClustalW and Muscle algorithms are both supported for the alignment.
"""

from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline
from io import StringIO


def align_seq(in_dir, method="clustalw", out_file=False):
    """Align multiple sequences using the MUSCLE algorithm.

    Parameters
    ----------
    in_dir : dir
        The directory to the file containing the sequences.
    method : str (default clustalw)
        The alignment algorithm, either clustalw or muscle.
    out_dir : boolean (default False)
        Wether or not a output file should be generated.

    Returns
    -------
    Bio.Align.MultipleSeqAlignment
        The aligned sequences.
    """
    clustalo_exe = "./clustal-omega/clustalo.exe"
    muscle_exe = "./muscle3.8.31_i86win32.exe"

    if out_file:
        out_dir = in_dir.replace(".fasta", "_aligned.fasta")
        open(out_dir, "w+").close()
        parameters = {"infile": in_dir, "outfile": out_dir, "force": True}
        if method == "clustalw":
            cline = ClustalOmegaCommandline(clustalo_exe, parameters)
        elif method == "muscle":
            cline = MuscleCommandline(muscle_exe, parameters)
        stdout, stderr = cline()
        align = AlignIO.read(out_dir, "fasta")

    else:
        parameters = {"infile": in_dir, "auto": True}
        if method == "clustalw":
            cline = ClustalOmegaCommandline(clustalo_exe, parameters)
        elif method == "muscle":
            cline = MuscleCommandline(muscle_exe, parameters)
        stdout, stderr = cline()
        align = AlignIO.read(StringIO(stdout), "fasta")
    return align


if __name__ == "__main__":
    SEQ_FILE = "data/spike_data_708.fasta"
    sequence = align_seq(SEQ_FILE, "clustalw", True)
