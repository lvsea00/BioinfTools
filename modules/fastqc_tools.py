import os


def check_gc_count(seq: str, gc_lower: float, gc_upper: float) -> bool:
    """
    Checks if GC percentage of the read matches given conditions
    """
    gc_content = (seq.count('G') + seq.count('C')) * 100 / len(seq)
    return gc_lower <= gc_content <= gc_upper


def check_len(seq: str, len_min: int, len_max: int) -> bool:
    """
    Checks if length of the read matches given conditions
    """
    return len_min <= len(seq) <= len_max


def check_quality(seq: str, quality_threshold: int) -> bool:
    """
    Checks if quality of the read matches given conditions
    """
    quality_sum: int = 0
    for char in seq:
        quality_sum += ord(char) - 33
    mean_quality = quality_sum/len(seq)
    return mean_quality > quality_threshold


def fastq_to_dict(input_fastq: str) -> dict:
    """
    Converts fastq file to dictionary.

    Args:
        input_fastq (str): Absolute path to fastq file.

    Returns:
        dict: Dictionary with fastq sequences.
    """
    fastq_seqs = dict()
    with open(input_fastq) as fastq_file:
        for line in fastq_file:
            if line.startswith('@'):
                seq_id = line.strip()
                seq = fastq_file.readline().strip()
                fastq_file.readline().strip()
                quality = fastq_file.readline().strip()
            fastq_seqs[seq_id] = (seq, quality)
    return fastq_seqs


def save_filtered(filtered_seqs: dict, output_fastq: str):
    """
    Save filtered fastq sequences from dictionary to the file.

    Args:
        filtered_seqs (dict): Dictionary with filtered fastq sequences
        output_fastq (str): Name of fastq file with filtered sequences.

    Returns:
        Program work result information.
    """
    if 'filtered' not in os.listdir():
        os.mkdir('filtered')
    filtered_path = os.path.join(os.getcwd(), 'filtered')
    if not os.path.isfile(os.path.join(filtered_path, output_fastq)):
        with open(os.path.join(filtered_path, output_fastq), mode='w') as fastq_file:
            for name, (seq, quality) in filtered_seqs.items():
                fastq_file.write(name + '\n' + seq + '\n' + '+' + name[1:] + '\n' + quality + '\n')
        print("Sequences are filtered out")
    else:
        print("Check the name of the output file! Risk of overwriting!")
