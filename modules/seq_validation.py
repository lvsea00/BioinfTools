from typing import Union


def length_check(sequences: list) -> Union[list, str]:
    """Checks the length of sequences list"""
    if len(sequences) == 1:
        return sequences[0]
    return sequences


def is_dna(seq: str) -> bool:
    """Checks if input seq is DNA"""
    dna = {"A", "T", "G", "C"}
    unique_chars = set(seq.upper())
    return unique_chars <= dna


def is_rna(seq: str) -> bool:
    """Checks if input seq is RNA"""
    rna = {"A", "U", "G", "C"}
    unique_chars = set(seq.upper())
    return unique_chars <= rna


def is_nucleic_acid(seq: str) -> bool:
    """Checks if input seq is nucleic acid"""
    return is_dna(seq) or is_rna(seq)