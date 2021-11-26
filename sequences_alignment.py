"""This module align a multi sequences fasta file.

ClustalW and Muscle algorithms are both supported for the alignment.
"""


from io import StringIO

from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline


def align_seq(in_dir, method="clustalw", out_file=False):
    """Align multiple sequences.

    Parameters
    ----------
    in_dir : str
        The directory to the file containing the sequences.
    method : str (default clustalw)
        The alignment algorithm, either clustalw or muscle.
    out_dir : boolean (default False)
        Wether or not an output file should be generated.

    Returns
    -------
    Bio.Align.MultipleSeqAlignment
        The aligned sequences.
    """
    if out_file:
        # create the output file or make sure it is empty
        out_dir = in_dir.replace(".fasta", "_aligned.fasta")
        open(out_dir, "w+").close()

        # call the correct alignement method via command line
        if method == "clustalw":
            cline =  ClustalOmegaCommandline(infile=in_dir,
                                             outfile=out_file,
                                             force=True)
        elif method == "muscle":
            cline = MuscleCommandline(input=in_dir,
                                      out=out_file)
        _ = cline()

        # read the output file
        align = AlignIO.read(out_dir, "fasta")

    else:
        # call the correct alignement method via command line
        if method == "clustalw":
            cline =  ClustalOmegaCommandline(infile=in_dir, auto=True)
        elif method == "muscle":
            cline = MuscleCommandline(input=in_dir)
        stdout, _ = cline()

        # read the command line output string
        align = AlignIO.read(StringIO(stdout), "fasta")
    return align


if __name__ == "__main__":
    SEQ_FILE = "data/spike_data_708.fasta"
    sequence = align_seq(SEQ_FILE, "clustalw", True)
