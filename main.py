from typing import Union
from modules.seq_validation import is_nucleic_acid, length_check
from modules.dna_rna_tools import transcribe, reverse, complement, reverse_complement
from modules.fastqc_tools import check_gc_count, check_len, check_quality, fastq_to_dict, save_filtered


def run_dna_rna_tools(*args) -> Union[list, str]:
    """Performs basic operations with DNA or RNA sequences"""
    *seqs, operation = args
    operations = {"is_nucleic_acid": is_nucleic_acid,
                  "transcribe": transcribe,
                  "reverse": reverse,
                  "complement": complement,
                  "reverse_complement": reverse_complement}
    result = []
    invalid_seq: str = "Invalid sequence in the input!"
    for seq in seqs:
        result.append(operations[operation](seq))
        if invalid_seq in result:
            break
    return length_check(result)


def filter_fastq(input_fastq: str,
                 output_fastq: str = 'output_fastq',
                 gc_bounds: Union[tuple, float] = (0, 100),
                 length_bounds: Union[tuple, int] = (0, 2**32),
                 quality_threshold: int = 0) -> str:
    """
    Filters reads based on QC.

    Args:
        input_fastq (str): Absolute path to fastq file with unfiltered sequences.
        output_fastq (str): Name of fastq file with filtered sequences.
        gc_bounds (tuple, float): GC percentage bounds.
        length_bounds (tuple, int): Read length bounds.
        quality_threshold (int): Min quality threshold.

    Returns:
        str: Program work result information.

    """
    fastq_seqs = fastq_to_dict(input_fastq)
    filtered_seqs = dict()
    gc_lower, gc_upper = gc_bounds if isinstance(gc_bounds, tuple) else (0, gc_bounds)
    len_min, len_max = length_bounds if isinstance(length_bounds, tuple) else (0, length_bounds)
    for name, (seq, quality) in fastq_seqs.items():
        if check_gc_count(seq, gc_lower, gc_upper) and check_len(seq, len_min, len_max) and check_quality(seq, quality_threshold):
            filtered_seqs[name] = (seq, quality)
    return save_filtered(filtered_seqs, output_fastq)
