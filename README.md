# 🧬 BioinfTools 
Bioinformatics tools for working with nucleic acid sequences

__Author: [Elizaveta Salnikova](https://github.com/lvsea00)__

## Main Features:
### Checking for sequence validity:

*Module seq_validation*
 - `is_nucleic_acid` - checks if input seq is nucleic acid
 - `is_dna` - checks if input seq is DNA
 - `is_rna` - checks if input seq is RNA
 - `length_check` - checks the length of sequences list
 
### Basic operations with nucleic acid sequences:

`run_dna_rna_tools` - performs basic operations with DNA or RNA sequences
- Accept any number of positional arguments
- Case-sensitive
- Accepts DNA or RNA sequences and the name of the procedure to be performed
- Returns string if a single sequence is submitted or a list of rows if several sequences are submitted

*Module dna_rna_tools*
 - `transcribe` - returns RNA sequence transcribed from DNA
 - `reverse`- returns reverse sequence
 - `complement` - returns complement sequence
 - `reverse_complement` - returns reverse complement sequence
 
### FastQC analysis

`filter_fastq` - filters reads based on QC
- Accepts 4 arguments as input:
  - *seqs* - a dictionary consisting of fastq sequences (see example_fastq in examples)
  - *gc_bounds* - the GC percentage for filtering ((0, 100) by default)
  - *length_bounds* - the length interval for filtering ((0.2**32) by default)
  - *quality_threshold* - the threshold value of the average read quality (0 by default (phred33))
- Returns a dictionary with sequences that match all given conditions
 
 *Module fastqc_tools*
 - `count_gc` - counts GC percentage of the read
 - `mean_quality` - counts mean quality of the read

Functions in modules can be used individually or as part of __run_dna_rna_tools__ and __filter_fastq__

## Examples
```
run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False 
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']

filter_fastq(example_fastq, gc_bounds = 60, length_bounds = (89, 100), quality_threshold = 30)
# {'@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
  'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD')}`
```

Example file for __filter_fastq__ function
```
example_fastq = {
    # 'name' : ('sequence', 'quality')
    '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079803': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD')
    }
``` 
