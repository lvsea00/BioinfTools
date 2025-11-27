# 🧬 BioinfTools 
__Bioinformatics tools for working with nucleic acid or protein sequences__

- Basic operations with nucleic acid sequences
- FastQC analysis
- Checking for sequence validity
- Reading bioinformatics files (fasta, fastq)

__Author: [Elizaveta Salnikova](https://github.com/lvsea00)__

## Main Features
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
- Accepts 5 arguments as input:
  - *input_fastq* - absolute path to fastq file with unfiltered sequences
  - *output_fastq* - name of fastq file with filtered sequences
  - *gc_bounds* - the GC percentage for filtering ((0, 100) by default)
  - *length_bounds* - the length interval for filtering ((0.2**32) by default)
  - *quality_threshold* - the threshold value of the average read quality (0 by default (phred33))
- Returns program work result information to notify user about performed operation
- Saves sequences that match all given conditions to the __output_fastq__ file in the directory __filtered__
 
 *Module fastqc_tools*
 - `check_gc_count` - checks if GC percentage of the read matches given conditions
 - `check_len` - checks if length of the read matches given conditions
 - `check_quality` - checks if quality of the read matches given conditions
 - `fastq_to_dict` - converts fastq file to dictionary 
 - `save_filtered` - saves filtered fastq sequences from dictionary to the file

### Reading bioinformatics files

 `convert_multiline_fasta_to_oneline` - converts multiline fasta to oneline fasta
 - Accepts 2 arguments as input:
   - *input_fasta* - absolute path to multiline fasta file
   - *output_fasta* - name of resulted oneline fasta file
- Parses the input fasta file in which the sequence (DNA/RNA/protein/...) can be split into several lines and saves it to a new fasta file in which each sequence fits into one line
   
 `parse_blast_output` - parses blast output file to extract the names of best matched sequences
 - Accepts 2 arguments as input:
   - *input_file* - absolute path to blast results file
   - *output_file* - name of the file with names of best matched sequence for each query
 - For each query selects the first row from the *Description* column and saves the set of obtained proteins to a new file in one column sorted alphabetically

## Additional information
- Functions in modules can be used individually or as part of __run_dna_rna_tools__ and __filter_fastq__
- Functions are user-friendly - prevent data loss, notify user of file overwriting risk
- Functions return program work result information to notify user about performed operation
  

## Examples
```
run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False 
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']

filter_fastq(input_fastq, output_fastq, gc_bounds = 60, length_bounds = (89, 100), quality_threshold = 30)
# 'Sequences are filtered out' if everything worked correctly

convert_multiline_fasta_to_oneline (input_fasta, output_fasta)
# 'Convert multiline fasta to oneline' if everything worked correctly

parse_blast_output(input_file, output_file)
# 'Extract the names of best matched sequences.' if everything worked correctly

If the output file already exists, functions will notify user before performing the operation.
#'Check the name of the output file! Risk of overwriting!'
```

