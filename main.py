import os
from typing import Union
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from abc import ABC, abstractmethod
import argparse
import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    style="{",
    format='{asctime}  {levelname}  {message}',
    handlers=[
        logging.FileHandler('fastq_filtrator.log', mode='w', encoding='utf-8'),
        logging.StreamHandler(sys.stderr)
    ]
)

class BiologicalSequence(ABC):
    def __init__(self, seq: str) -> None:
        self.seq = seq

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, index: int) -> str:
        return self.seq[index]

    def __str__(self) -> str:
        return f"{self.__class__.__name__}: {self.seq}"

    @abstractmethod
    def validate(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    complement_dict = None

    def __init__(self, seq: str):
        super().__init__(seq)
        if not self.validate():
            raise ValueError("Invalid sequence")

    def validate(self) -> bool:
        for char in self.seq:
            if char not in self.alphabet:
                return False
        return True

    def reverse(self) -> 'NucleicAcidSequence':
        """Returns reverse sequence"""
        return self.__class__(self.seq[::-1])

    def complement(self) -> 'NucleicAcidSequence':
        """Returns complement sequence"""
        if self.complement_dict is None:
            raise NotImplementedError
        complement_seq = ''.join(self.complement_dict[char] for char in self.seq)
        return self.__class__(complement_seq)

    def reverse_complement(self) -> 'NucleicAcidSequence':
        """Returns reverse complement sequence"""
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    alphabet = set('ATGC')
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    def transcribe(self) -> 'RNASequence':
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
            raise ValueError("Invalid sequence")

    def validate(self) -> bool:
        for char in self.seq:
            if char not in self.alphabet:
                return False
        return True

    def is_starts_from_met(self) -> bool:
        return self.seq[0].upper() == 'M'


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

    filtered_count = 0
    total_count = 0

    try:
        with open(input_fastq, 'r') as input_file, open(output_path, mode='w') as output_file:
            records = SeqIO.parse(input_file, "fastq")
            logging.info(f"Starting filtering: {input_fastq}")
            for record in records:
                total_count += 1
                seq = record.seq
                gc_content = gc_fraction(seq) * 100
                qualities = record.letter_annotations["phred_quality"]
                mean_quality = sum(qualities) / len(qualities)

                if (len_min <= len(seq) <= len_max) and (gc_lower <= gc_content <= gc_upper) and (mean_quality > quality_threshold):
                    SeqIO.write(record, output_file, "fastq")
                    filtered_count += 1

        logging.info(f"Filtered {filtered_count} of {total_count} sequences from {input_fastq}")
    
    except Exception as e:
        logging.error(f"Error to process {input_fastq}: {e}", exc_info=True)
        raise

def main():
    parser = argparse.ArgumentParser(
                    prog='Fastq Filtrator',
                    description='Filters reads based on QC',
                    epilog='')
    
    parser.add_argument('input_fastq', type=str, help='Path to input FASTQ file')

    parser.add_argument('-o', '--output_fastq', type=str, default='output_fastq',
                       help='Output fastq filename (default: output_fastq)')
    
    parser.add_argument('-gc', '--gc_bounds', type=float, nargs=2, default=[0, 100],
                   metavar=('MIN', 'MAX'), help='GC percentage bounds (default: 0-100)')
    
    parser.add_argument('-len', '--length_bounds', type=int, nargs=2, default=[0, 2**32],
                   metavar=('MIN', 'MAX'), help='Read length bounds (default: 0-2**32)')
    
    parser.add_argument('-q', '--quality_threshold', type=int, default=0,
                       help='Min quality threshold (default: 0)')
    
    args = parser.parse_args()
    gc_bounds = tuple(args.gc_bounds)
    length_bounds = tuple(args.length_bounds)

    filter_fastq(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=args.quality_threshold)


if __name__ == "__main__":
    main()