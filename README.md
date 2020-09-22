# snpdiffrs

snpdiffrs finds differential single nucleotide polymorphisms between BAM
files (of the same species).

## Usage - N to N comparison

Create an `input.toml` containing, at minimum, the following

```toml
output_dir = 'tests/A_vs_B'
[samples]
    A = ['path_to/A.bam']
    B = ['path_to/B.bam']
```

and run using `snpdiffrs input.toml`.

See Modes for details and other run modes..

## Output
The output files are tab-seperated-value files,
with one row per SNP found.

|Column   |Description|
|---|---|
|chr| chromosome/BAM-reference   |
|pos| reference position (0-based)  |
|score| -1 * log<sub>e</sub>-likelihood according to internal model |
|A_A| Count of adenine in 1st sample |
|B_A| Count of adenine in 2nd sample |
|A_C| Count of cytosine in 1st sample|
|B_C| Count of cytosine in 2nd sample|
|A_G| Count of guanine in 1st sample |
|B_C| Count of guanine in 2nd sample |
|A_T| Count of thymine in 2nd sample |
|B_T| Count of thymine in 1st sample |
|haplotypeA| Most likely haplotype (diploid) of 1st sample  |
|haplotypeB| Most likely haplotype (diploid) of 2nd sample  |


## Options

### Top level options

| Option | Type | Description | 
|---|---|---|
|mode | enum | n\_to\_n, 1\_to\_n, quantify\_homopolymers. Default: n\_to\_n. See modes  | 
output\\_dir  | string | Write output files here. Directory will be created if missing. Files will be overwritten.
chromosomes| list | - restrict comparison to these chromosomes/BAM-references|
quality\_threshold | float | Base quality threshold (PHRED) applied before counting. Default 15.|
filter\_homo\_polymer\_threshold | int | Exclude reads that have a homopolymer of length >= filter\_homo\_polymer\_threshold. Default: No such filtering. |
min\_score | float | Minimum score (see output) to filter for. Default 50.0
ncores | int | How many CPU cores to use. Default: all of them.|
block\_size | int | Size of coverage region examined at once. Trade-off between disk seeks, RAM usage. Setting this too high will impair parallel computation though. Default 50\_000\_000|



## Modes

### n\_to\_n

Example:
```toml
mode = 'n_to_n'
output_dir = 'tests/A_vs_B'
[samples]
    A = ['path_to/A.bam']
    B = ['path_to/B.bam']
```


All samples will be compared vs each other, 
and you will receive one output file A_vs_B.tsv
in your output directory.

Each pair will be present only once, with the
sample names sorted in lexographical order.

### one\_to\_n

Example: 
```toml
mode = 'OneToN'
output_dir = 'tests/X_vs_all'
[query]
	X = ['path_to/X.bam', 'path_to/X2.bam'] # multiple BAMs per sample will be aggregated
[references]
	B = [...]
	C = [...]
	D = [...]
```

Compares one sample vs all others (e.g. to 
determine contamination).
Output files named like X_vs_B.tsv, X_vs_C.tsv, etc.


### quantify\_homopolymers

Example: 
```toml
mode = 'QuantifyHomopolymers'
output_dir = 'tests/quantification'
[samples]
	A = ['A.bam'] # multiple BAMs per sample will be aggregated
	B = ['B.bam'] # multiple BAMs per sample will be aggregated
	...
```

Provide a tab-seperated-value of reads-with-a-(maximum)-homopolymer-length
histogram per sample named like `A_homopolymer_histogram.tsv`.

Some sequencing runs apperantly create lot's of 'fake' homopolymers, 
which lead to massive false SNP calls. You can try to use the QuantifyHomopolymers
mode to gain insight, I've had better luck simply comparing a run with
a modest filter (e.g. filter homopolymers length >= 8) vs no homopolymer filtering,
and looking at the ratio of called snps for sample pairs.



