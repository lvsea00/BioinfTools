import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str) -> str:
    """
        Converts multiline fasta to oneline fasta.

    Args:
        input_fasta (str): Absolute path to multiline fasta file.
        output_fasta (str): Name of resulted oneline fasta file.

    Returns:
        str: Program work result information.
    """
    output_fasta_path = os.path.join(os.getcwd(), output_fasta)
    if not os.path.isfile(output_fasta_path):
        with open(input_fasta) as input_fasta, open(output_fasta, mode='a+') as output_fasta:
            seq = ''
            for line in input_fasta:
                line = line.strip()
                if line.startswith('>'):
                    if len(seq) != 0:
                        output_fasta.write(seq+'\n')
                        seq = ''
                    seq_id = line
                    output_fasta.write(seq_id+'\n')
                else:
                    seq += line
            output_fasta.write(seq+'\n')
        return "Convert multiline fasta to oneline"
    else:
        return "Check the name of the output file! Risk of overwriting!"
