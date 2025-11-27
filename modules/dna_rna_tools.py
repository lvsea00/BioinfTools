from modules.seq_validation import is_nucleic_acid, is_rna, is_dna


def transcribe(seq: str) -> str:
    """Returns RNA sequence transcribed from DNA"""
    transcribed_nucl = {"A": "A", "T": "U", "G": "G", "C": "C",
                        "a": "a", "t": "u", "g": "g", "c": "c"}
    if is_dna(seq):
        transcribed_seq = "".join([transcribed_nucl[char] for char in seq])
    else:
        transcribed_seq = "Invalid sequence in the input!"
    return transcribed_seq


def reverse(seq: str) -> str:
    """Returns reverse sequence"""
    if is_nucleic_acid(seq):
        reverse_seq = seq[::-1]
    else:
        reverse_seq = "Invalid sequence in the input!"
    return reverse_seq


def complement(seq: str) -> str:
    """Returns complement sequence"""
    if is_dna(seq):
        complement_seq = ""
        complements = {"A": "T", "T": "A", "G": "C", "C": "G",
                       "a": "t", "t": "a", "g": "c", "c": "g"}
        for char in seq:
            complement_seq += complements[char]
    elif is_rna(seq):
        complement_seq = ""
        complements = {"A": "U", "U": "A", "G": "C", "C": "G",
                       "a": "u", "u": "a", "g": "c", "c": "g"}
        for char in seq:
            complement_seq += complements[char]
    else:
        complement_seq = "Invalid sequence in the input!"
    return complement_seq


def reverse_complement(seq: str) -> str:
    """Returns reverse complement sequence"""
    rev_compl_seq = reverse(complement(seq))
    return rev_compl_seq