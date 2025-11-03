import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str):
    """
    Converts multiline fasta to oneline fasta.

    Args:
        input_fasta (str): Absolute path to multiline fasta file.
        output_fasta (str): Name of resulted oneline fasta file.

    Returns:
        Program work result information.
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
        print("Convert multiline fasta to oneline")
    else:
        print("Check the name of the output file! Risk of overwriting!")


def parse_blast_output(input_file: str, output_file: str):
    """
    Parses blast output file to extract the names of best matched sequences.

    Args:
        input_file (str): Absolute path to blast results file.
        output_file (str): Name of the file with names of best matched sequence for each query.

    Returns:
        Program work result information.
    """
    output_file_path = os.path.join(os.getcwd(), output_file)
    if not os.path.isfile(output_file_path):
        with open(input_file) as input_file, open(output_file, mode='w') as output_file:
            best_matches = []
            for line in input_file:
                if line.startswith('Description'):
                    char_count = 0
                    for char in line:
                        char_count += 1
                        if char == 'N':
                            break
                    line = input_file.readline()
                    line = (line[:(char_count-1)]).strip()
                    best_matches.append(line)
            for name in sorted(best_matches):
                output_file.write(name+'\n')
        print("Extract the names of best matched sequences.")
    else:
        print("Check the name of the output file! Risk of overwriting!")
