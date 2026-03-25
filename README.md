# Mandalorion #
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](http://perso.crans.org/besson/LICENSE.html)


v1.0 - I have a bad feeling... about this isoform

IG28 high accuracy full-length cDNA reads produced by ONT/R2C2 and/or PacBio/Kinnex in fasta/fastq/bam formats and identifies and quantifies high confidence isoforms.
You can mix and match R2C2/PacBio reads.

IG28 has not been tested for use with regular ONT cDNA reads. You can certainly try though!

## Short Example 

```bash
IG28 expects polyA trimmed reads in the 5->3 direction. If your reads aren't already polyA trimmed. You can remove them using the
remove_polyA.py script in utils/

Look at the reads in your hypothetical input file. 
If your reads look like this (some untrimmed adapter sequences, polyA tail still present)

CATGTGTATAAGAGACAGTAGTACTTGGTGGGG [RNA transcript 1]  AAAAAAAAAAAAAAAAAAAAAAGTGAACTGAGCGACTCTGCGT
CATGTGTATAAGAGACAGTAGTACTTGGTGGGG [RNA transcript 2]  AAAAAAAAAAAAAAAAAAAAAAGTGAACTGAGCGACTCTGCGT
CATGTGTATAAGAGACAGTAGTACTTGGTGGGG [RNA transcript 3]  AAAAAAAAAAAAAAAAAAAAAAGTGAACTGAGCGACTCTGCGT

|-------------33nt--------------|                                           |--------21nt-------|

remove untrimmed adapter sequences and polyA tails like this:

python3 IG28/utils/removePolyA.py -i input.5to3.fasta -o input.5to3.noPolyA.fasta -t 33,21

after trimming adapters and polyA sequences:

python3 IG28/IG28 -p ./ -g gencodeV49.gtf -G hg38.fasta -f input.5to3.noPolyA.fofn
```

IG28 will generate 

Isoform  models (.gtf) and quantification (.quant for read counts and .rpms for RPM values.

## Overview ##

IG28 contains several modules that are run sequentially by default (ADC) to define, align, and quantify high-quality isoforms.

## Download ##
Just navigate to the folder you want Mandalorion to live and git clone the repository

```bash
git clone https://github.com/christopher-vollmers/IG28.git
```

## Dependencies ##

To make IG28 executable without explicitly calling it using python3 you can change its permissions.

```bash
cd IG28
chmod +x IG28
```

You will need to install the following python3 libraries. If your python3 is <=3.6 you might be in trouble because of f-string formatting and 
fifo stuff that later version of python3 do that we might have inadvertendly relied on.  

numpy
pysam
mappy

You will also need to install 

minimap2
samtools
featureCounts (if using Module G for doing single cell analysis)
TD2 (if you want to add CDS coordinates to your gtf output)


in your PATH. 


## PolyA tail removal ##

IG28 requires reads that have no adapter sequences (trimmed), no polyA tails, and are in 5'->3' orientation. By default, C3POa (v3) does trim reads but doesn't do polyA removal and, depending on the adapters used, may not put reads in 5'->3' orientation..   
Therefore the removePolyA.py script in utils/ should be used to get both PacBio (depending on preprocessing) and R2C2 reads ready.
The script first removes a fixed number of bases (can't be zero because we aren't very good at programming) from the ends of the read and then removes the polyA tail on the 3' end of the read. Depending on the polyA position it will also reorient the read in 5'->3' direction.


## Detailed usage ##
You can call IG28 from the IG28 folder or from anywhere else so either

```bash
cd /path/to/IG28/
IG28 [OPTIONS]
```
or

```bash
python3 /path/to/IG28/IG28 [OPTIONS]
```

Running with default settings:
```bash
IG28 -p . -g gencodeV29.gtf -G hg38.fasta -f Consensus_reads.fofn
```

.fofn file structure is simply a text file with one line per input fasta/fastq/bam file that contains the path to the file and a tab separated sample name and a tab separated tissue name.
Reads from each unique sample will be combined for quantification. Reads from each tissue will be used to make separate isoform gtf/quant/rmps files -- in addition to a gtf/quant/rpms using all reads. 


You can mix and match fasta/fastq that are gzipped or not as well as bam files. 

```bash
/path/to/file1.fasta	sample1	tissue1
/path/to/file2.bam	sample2	tissue2
/path/to/file4.fastq.gz	sample1	tissue1
```

Here is a full list of options that you can use to modify IG28 behavior:

### Options: ###

```bash

  -p PATH, --path PATH  Directory to put output files into
  -u UPSTREAM_BUFFER, --upstream_buffer UPSTREAM_BUFFER
                        Defines upstream leniency window for polyA and TSS determination (default 10)
  -d DOWNSTREAM_BUFFER, --downstream_buffer DOWNSTREAM_BUFFER
                        Defines downstream leniency window for polyA and TSS determination (default 50)
  -g GENOME_ANNOTATION, --genome_annotation GENOME_ANNOTATION
                        Genome annotation file (gtf). Is used to identify individual annotated splice sites. If -W is set it is also used to white-list annotated polyA sites
  -G GENOME_SEQUENCE, --genome_sequence GENOME_SEQUENCE
                        Genome file (fasta)
  -r MINIMUM_RATIO, --minimum_ratio MINIMUM_RATIO
                        Proportion of reads that align to a locus required for an isoform (default 0.001)
  -i MINIMUM_INTERNAL_RATIO, --minimum_internal_ratio MINIMUM_INTERNAL_RATIO
  -R MINIMUM_READS, --minimum_reads MINIMUM_READS
                        Minimum number of reads to make an isoform (default 5)
  -f CONSENSUS_READS, --Consensus_reads CONSENSUS_READS
                        Fasta/fastq file with R2C2/PacBio consensus reads, can be entered as a single file path, a comma separated list of file paths, or a path to a file of filenames file (has to end on .fofn) that contains one file path
                        per line and a tab separated sample name
  -O OVERHANGS, --overhangs OVERHANGS
                        Defines bounds for unaligned bases on ends. Format: min5prime,max5prime,min3prime,max3prime (default 0,40,0,40)
  -t MINIMAP2_THREADS, --minimap2_threads MINIMAP2_THREADS
                        Number of threads to use when running minimap and consensus calling (default 8)
  -I MINIMUM_ISOFORM_LENGTH, --minimum_isoform_length MINIMUM_ISOFORM_LENGTH
                        Minimum length in nt of isoforms that will be considered (default 200)
  -n MINIMUM_FEATURE_COUNT, --minimum_feature_count MINIMUM_FEATURE_COUNT
                        Features (starts,ends, new splice sites) will be considered if they are present in this number of reads (default 2)
  -w SPLICE_SITE_WINDOW, --splice_site_window SPLICE_SITE_WINDOW
                        reads spliced within this number of nucleotides on each side of a splice site will be considered spliced at this site (default 1)
  -A ACUTOFF, --Acutoff ACUTOFF
                        Isoforms with A content of more than this cutoff in a 30nt window surrounding their polyA site will be discarded (default 0.5)
  -W WHITE_LIST_POLYA, --white_list_polyA WHITE_LIST_POLYA
                        If set, polyA sites that fall within +/-20nt of annotated transcript ends will not be filtered regardless of Acutoff set with -A. Annotated transcript ends will be taken from annotation file given with -g. Only
                        transcripts having one of the provided comma separated values in the line of their exon features will be used. E.g. setting -W to [SIRV,"basic"] will whitelist spike-in SIRV transcripts and transcripts considered
                        "basic", i.e. high confidence full length, in the gencode annotation. Setting -W to [exon] should include all transcripts in the gtf file. This is feature only checks whether the provided values are in the line, not
                        where they are. That means that setting -W to [chr1] will whitelist all transcripts on chromosome 1 but also transcripts with "chr1" in their name, so it's a bit dangerous
  -x META_DATA, --meta_data META_DATA
                        meta_data file with two columns: cell and celltype
  -m, --multi_exon_only
                        If used, IG28 will filter all single exon isoforms
  -j JUNCTIONS, --junctions JUNCTIONS
                        base context required (in lower-case for either strand) to call unannotated splice junctions (default gtag,gcag,atac,ctac,ctgc,gtat)
  -M MODULES, --Modules MODULES
                        Defines what modules of IG28 will be run. By default this includes: A - Alignment, D - defining isoforms, C - Adding CDS G - Make single cell gene expression matrix M - Modify meta data (mostly a single cell thing)
                        Each module needs the output of the modules run before it to function properly. Running individual modules can be useful if you want to for example filter with different parameters without rerunning the whole
                        pipeline. (default ADC)
  -v, --version         Prints IG28 version

```

## Outputs ##

When running Modules A and D back to back this is the output directory structure IG28 produces

output_dir
├── Combined
│   ├── Isoforms.filtered.clean.gtf
│   ├── Isoforms.filtered.clean.genes.quant
│   └── Isoforms.filtered.clean.genes.rpms
├── Tissue1
│   ├── Isoforms.filtered.clean.gtf
│   ├── Isoforms.filtered.clean.genes.quant
│   └── Isoforms.filtered.clean.genes.rpms
└── Tissue2
    ├── Isoforms.filtered.clean.gtf
    ├── Isoforms.filtered.clean.genes.quant
    └── Isoforms.filtered.clean.genes.rpms

For each tissue and all data combined, IG28 will generate the following output files in their respective output directories

Isoforms.filtered.clean.gtf: Isoform structures with small indels removed

Isoforms.filtered.clean.genes.quant: Number of reads from each sample associated with each isoform

Isoforms.filtered.clean.genes.rpms: Normalized expression values for each sample associated with each isoform

Isoforms.filtered.clean.genes: Isoform groups and gene association


#When running Module C as well an additional file#


Isoforms.filtered.clean.CDS.gtf 

will be added to each directory which adds CDS coordinates to each isoform by either lifting them from the supplied annotation or calling them using TD2

##Single Cell data ##

IG28 was designed with full-length cDNA single cell data in mind. 10X, Parse, etc... doesn't really matter. 

The workflow we envisioned for single cell data is as follows:

Your data should be in an (unaligned) bam file. It should be preprocessed to contain a CB (or equivalent) tag indicating what cell each read belongs to.
You sample sheet.fofn should be structured as follows

```bash
/path/to/sc-flcDNA1.bam	scRNA;flcDNA1;CB	None
/path/to/sc-flcDNA2.bam	scRNA;flcDNA2;CB	None
```

First column should contain the absulte path to your (unaligned) bam files. 
Second column should contain 'scRNA' to tell IG28 that this is single cell data and it has to parse this column. 
In the first line 'flcDNA1' would become the root name for each cell followed by what was in the bam file under the 'CB' tag
If the bam file has a different tag, just change it in the sample sheet.
 
#Run Module A to parse and align reads#

After you have a properly fomatted sample sheet (.fofn) file, run Module A (-M A) which parse and then uses minimap2 to align the reads.

#Run Module G to generate a genexcell expression matrix


After Module A finishes, run Module G (-M G) so IG28 calls featureCounts to quantify gene expression counts.
You will need to have featureCounts in your PATH for this module

```bash
IG28 -p ./ -g gencodeV49.gtf -G hg38.fasta -f scRNA_sampleSheet.fofn -M G
```

The matrix it produces can go straight in Seurat or what ever your scRNA analysis pipeline of choice is. 

We have a vignette for a public PBMC Kinnex data set here.

At this point, you can cluster your cells into celltypes in whatever way you please - go nuts. 

The only thing you will need your scRNA workflow to do is produce a metadata tab delimited file with the following format:

```bash
Cell_ID1	Celltype1
Cell_ID2	Celltype2
Cell_ID3	Celltype1
```
#Run Module M to modify IG28 metadata

After Module G finishes, and you ran your scRNA workflow of choice and generate a metadata file you can run Module M (-M M) to modify IG28 metadata.

```bash
IG28 -p ./ -g gencodeV49.gtf -G hg38.fasta -f scRNA_sampleSheet.fofn -M M -x metadata.tsv
```

In that way you prepare IG28 to create independent transcriptomes for each cell type in your data. 

#Run Module D and C as you would for bulk data

Once the metadata is modified, just run Module D as you would for bulk data. The output will look like this

output_dir
├── Combined
│   ├── Isoforms.filtered.clean.gtf
│   ├── Isoforms.filtered.clean.genes.quant
│   └── Isoforms.filtered.clean.genes.rpms
├── Celltype1
│   ├── Isoforms.filtered.clean.gtf
│   ├── Isoforms.filtered.clean.genes.quant
│   └── Isoforms.filtered.clean.genes.rpms
└── Celltype2
    ├── Isoforms.filtered.clean.gtf
    ├── Isoforms.filtered.clean.genes.quant
    └── Isoforms.filtered.clean.genes.rpms

The quant and rpms files will contain quantifications for every single cell.

## Utils ##

There is also a script to plot a genome browser shot as shown here. 

![test](https://github.com/christopher-vollmers/Mandalorion/assets/28308271/b44b972a-a9b9-4040-8f04-337647e59ada)
