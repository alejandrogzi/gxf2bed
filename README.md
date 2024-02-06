# gxf2bed
A high-performance bed-to-gtf converter written in Rust.

translates

```
chr27 gxf2bed gene 17266470 17285418 . + . gene_id "ENSG00000151743";

chr27 gxf2bed transcript 17266470 17281218 . + . gene_id "ENSG00000151743"; transcript_id "ENST00000541931.8";

chr27 gxf2bed exon 17266470 17266572 . + . gene_id "ENSG00000151743"; transcript_id "ENST00000541931.8"; exon_number "1"; exon_id "ENST00000541931.8.1";

...
```
into

```
chr27 17266469 17281218 ENST00000541931.8 1000 + 17266469 17281218 0,0,200 2 103,74, 0,14675,
```
in a few seconds.


Converts
- *Homo sapiens* GRCh38 GENCODE 44 (252,835 transcripts) in 3.34 seconds.
- *Mus musculus* GRCm39 GENCODE 44 (149,547 transcritps) in 2.11 seconds.
- *Canis lupus familiaris* ROS_Cfam_1.0 Ensembl 110 (55,335 transcripts) in 0.99 seconds. 
- *Gallus galus* bGalGal1 Ensembl 110 (72,689 transcripts) in 1.14 seconds.


## Usage
``` rust
Usage: gxf2bed[EXE] --input/-i <GTF/GFF> --output/-o <BED>
 
Arguments:
    --input/-i <GTF/GFF>: a .gtf/.gff file
    --output/-o <BED>: path to output .bed file

Options:
    --help: print help
    --version: print version
    --threads/-t: number of threads (default: max ncpus)
```

## Installation
to install gxf2bed on your system follow this steps:
1. get rust: `curl https://sh.rustup.rs -sSf | sh` on unix, or go [here](https://www.rust-lang.org/tools/install) for other options
2. run `cargo install gxf2bed` (make sure `~/.cargo/bin` is in your `$PATH` before running it)
4. use `gxf2bed` with the required arguments
5. enjoy!

## Build
to build gxf2bed from this repo, do:

1. get rust (as described above)
2. run `git clone https://github.com/alejandrogzi/gxf2bed.git && cd gxf2bed`
3. run `cargo run --release -- -i <GTF/GFF> -o <BED> 

## Container image
to build the development container image:
1. run `git clone https://github.com/alejandrogzi/gxf2bed.git && cd gxf2bed`
2. initialize docker with `start docker` or `systemctl start docker`
3. build the image `docker image build --tag gxf2bed .`
4. run `docker run --rm -v "[dir_where_your_gtf_is]:/dir" gxf2bed -i /dir/<GTF/GFF> -o /dir/<BED>`

## Conda
to use gxf2bed through Conda just:
1. `conda install gxf2bed -c bioconda` or `conda create -n gxf2bed -c bioconda gxf2bed`
