# miRsnp

A Python package for detecting the effect of gene mutation (SNP/InDel) on microRNA binding target genes

## Features

+ Create mutant sequences according reference and mutations
+ Run `RNAhybrid` pipeline separately with mutant sequences and reference
+ Stats the difference on microRNA binding target genes between mutant sequences and reference

## Installation

+ Copy the package to a Linux system path, Python-3.6 required to be installed and other packages refer to `./requirements.txt`  
+ Software `RNAhybrid` need to be installed and add it to ENV  

## Run

+ refer to `./tests/run.sh`

    ```
    cd ./tests/

    ~/software/miniconda3/bin/python3 ../scripts/mir_snp.py \
        -i ../data/variants_input.txt \
        -f ../data/mirna.fa \
        -r ../data/Homo_sapiens.GRCh38.101.exon-example.fa \
        -g ../data/Homo_sapiens.GRCh38.101-example.gtf \
        -o ./output
    ```

## Input

+ Variants input file (`CSV` and `EXCEL` were supported too), format like this：

    ```
    rsid    chr pos ref alt gene
    rs1051042   12  112919432   G   C   OAS1
    rs1131476   12  112919404   G   A   OAS1
    ```

    for example: `./data/variants_input.txt`  
  
+ Query miRNA sequences(FASTA format):

    ```
    >1364 
    CGCAGGGAATGGATAATCTTGC
    ```

    for example: `./data/mirna.fa`  

+ Reference transcripts(FASTA format), complete file can be created with the command below:

    ```
    gffread -w Homo_sapiens.GRCh38.101.exon.fa -W -g \
    /path/to/Homo_sapiens.GRCh38.101.genome.fa \
    /path/to/Homo_sapiens.GRCh38.101.gtf
    ```

    for example: `./data/Homo_sapiens.GRCh38.101.exon-example.fa`  

+ GTF：  
    Downloaded from Ensembl Genomes  
    for example: `./data/Homo_sapiens.GRCh38.101-example.gtf`  

## Output

    ```
    ./tests/output/result.txt
    ```

## Contact

Liuke Yang (892297359@qq.com)
