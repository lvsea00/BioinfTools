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
        


dna = DNASequence("ATGC")
rna = RNASequence("AUGC")
protein = AminoAcidSequence("MVLSPAD")

# Проверяем методы
print(dna)                          # DNASequence: ATGC
print(dna.complement())              # DNASequence: TACG
print(dna.reverse())                 # DNASequence: CGTA  
print(dna.reverse_complement())      # DNASequence: GCAT

print(dna.transcribe())              # RNASequence: AUGC

print(rna.complement())              # RNASequence: UACG

print(protein)                       # AminoAcidSequence: MVLSPAD
print(protein.met_start()) 



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





from typing import Union


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
        if check_gc_count(seq, gc_lower, gc_upper) and check_len(seq, len_min, len_max) and check_quality(quality, quality_threshold):
            filtered_seqs[name] = (seq, quality)
    return save_filtered(filtered_seqs, output_fastq)
