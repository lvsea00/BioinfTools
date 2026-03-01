# 🧬 BioinfTools 
__Bioinformatics tools for working with nucleic acid or protein sequences__

- Basic operations with nucleic acid sequences
- FastQC analysis
- Checking for sequence validity
- Reading bioinformatics files (fasta, fastq)

__Author: [Elizaveta Salnikova](https://github.com/lvsea00)__

## Main Features
### Basic operations with nucleic acid sequences and amino acid sequences:
- `validate` - checks the sequence validity
- `reverse` - returns reverse sequence
- `complement` - returns complement sequence
- `reverse_complement` - returns reverse complement sequence
- `transcribe` - returns RNA sequence transcribed from DNA
- `is_start_from_met` - checks if the amino acid sequence starts with M
  
### FastQC analysis
`filter_fastq` - filters reads based on QC
- Accepts 5 arguments as input:
  - *input_fastq* - path to fastq file with unfiltered sequences
  - *output_fastq* - name of fastq file with filtered sequences
  - *gc_bounds* - the GC percentage for filtering ((0, 100) by default)
  - *length_bounds* - the length interval for filtering ((0.2**32) by default)
  - *quality_threshold* - the threshold value of the average read quality (0 by default (phred33))
- Returns program work result information to notify user about performed operation
- Saves sequences that match all given conditions to the __output_fastq__ file in the directory __filtered__

## Examples
```
dna = DNASequence("ATGC")
rna = RNASequence("AUGC")
amino_acid = AminoAcidSequence("MVRPL")

print(dna)  # DNASequence: ATGC
print(dna.reverse())  # DNASequence: CGTA
print(dna.complement())  # DNASequence: TACG
print(dna.reverse_complement()) # DNASequence: GCAT
print(rna.complement()) # RNASequence: UACG
print(amino_acid.is_starts_from_met()) # True

filter_fastq(input_fastq, output_fastq, gc_bounds = 60, length_bounds = (89, 100), quality_threshold = 30)
# 'Sequences are filtered out' if everything worked correctly

```

