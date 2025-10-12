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

    Arg:
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
