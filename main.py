from typing import Union
from modules.seq_validation import is_nucleic_acid, length_check
from modules.dna_rna_tools import transcribe, reverse, complement, reverse_complement
from modules.fastqc_tools import count_gc, mean_quality


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


def filter_fastq(seqs: dict[str, tuple],
                 gc_bounds: Union[tuple, float] = (0, 100),
                 length_bounds: Union[tuple, int] = (0, 2**32),
                 quality_threshold: int = 0) -> dict:
    """
    Filters reads based on QC.
    Returns dict with reads that match all given conditions
    """
    qc_seqs = dict()
    for key, value in seqs.items():
        if type(gc_bounds) is tuple:
            gc_bool = gc_bounds[0] <= count_gc(value[0]) <= gc_bounds[1]
        else:
            gc_bool = count_gc(value[0]) < gc_bounds
        if type(length_bounds) is tuple:
            len_bool = length_bounds[0] <= len(value[0]) <= length_bounds[1]
        else:
            len_bool = len(value[0]) < length_bounds
        if gc_bool and len_bool:
            if mean_quality(value[1]) > quality_threshold:
                qc_seqs[key] = value
    return qc_seqs
