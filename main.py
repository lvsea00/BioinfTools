from abc import ABC, abstractmethod

class BiologicalSequence(ABC):
    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)
    
    def __getitem__(self, index):
        return self.seq[index]

    def __str__(self):
        return f"{self.__class__.__name__}: {self.seq}"
    
    @abstractmethod
    def validate(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    complement_dict = None

    def __init__(self, seq):
        super().__init__(seq)
        if not self.validate():
            raise ValueError(f"Invalid sequence")
    
    def validate(self):
        for char in self.seq:
            if char not in self.alphabet:
                return False
        return True

    def reverse(self):
        """Returns reverse sequence"""
        return self.__class__(self.seq[::-1])

    def complement(self):
        """Returns complement sequence"""
        if self.complement_dict is None:
            raise NotImplementedError
        complement_seq = ''.join(self.complement_dict[char] for char in self.seq)
        return self.__class__(complement_seq)

    def reverse_complement(self):
        """Returns reverse complement sequence"""
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    alphabet = set('ATGC')
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    def transcribe(self):
        """Returns RNA sequence transcribed from DNA"""
        transcribed_seq = self.seq.replace('T', 'U')
        return RNASequence(transcribed_seq)


class RNASequence(NucleicAcidSequence):
    alphabet = set('AUGC')
    complement_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


class AminoAcidSequence(BiologicalSequence):
    alphabet = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'O', 'S', 'U', 'T', 'W', 'Y', 'V'}

    def __init__(self, seq):
        super().__init__(seq)
        if not self.validate():
            raise ValueError(f"Invalid sequence")
    
    def validate(self):
        for char in self.seq:
            if char not in self.alphabet:
                return False
        return True           
    
    def met_start(self):
        return self.seq[0].upper() == 'M'
        


import os
from typing import Union
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def filter_fastq(input_fastq: str,
                 output_fastq: str = 'output_fastq',
                 gc_bounds: Union[tuple, float] = (0, 100),
                 length_bounds: Union[tuple, int] = (0, 2**32),
                 quality_threshold: int = 0) -> str:
    """
    Filters reads based on QC.

    Args:
        input_fastq (str): Path to fastq file with unfiltered sequences.
        output_fastq (str): Name of fastq file with filtered sequences.
        gc_bounds (tuple, float): GC percentage bounds.
        length_bounds (tuple, int): Read length bounds.
        quality_threshold (int): Min quality threshold.

    Returns:
        str: Program work result information.

    """
    
    gc_lower, gc_upper = gc_bounds if isinstance(gc_bounds, tuple) else (0, gc_bounds)
    len_min, len_max = length_bounds if isinstance(length_bounds, tuple) else (0, length_bounds)

    output_dir = 'filtered'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, output_fastq)

    with open(input_fastq, 'r') as input_file, open(output_path, mode='w') as output_file:
        records = SeqIO.parse(input_file, "fastq")
        for record in records:
            seq = record.seq
            gc_content = gc_fraction(seq) * 100
            qualities = record.letter_annotations["phred_quality"]
            mean_quality = sum(qualities) / len(qualities)

            if (len_min <= len(seq) <= len_max) and (gc_lower <= gc_content <= gc_upper) and (mean_quality > quality_threshold):
                SeqIO.write(record, output_file, "fastq")

    print("Sequences are filtered out")
                
                      




