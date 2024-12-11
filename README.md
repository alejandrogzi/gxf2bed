<p align="center">
  <h1 align="center">
    gxf2bed
  </h1>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.1.0dev-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.2.3-green">
    </a>
    <a href="https://crates.io/crates/gxf2bed" target="_blank">
      <img alt="Crates.io Version" src="https://img.shields.io/crates/v/gxf2bed">
    </a>
    <a href="https://github.com/alejandrogzi/gxf2bed" target="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/alejandrogzi/gxf2bed?color=blue">
    </a>
    <a href="https://crates.io/crates/gxf2bed" target="_blank">
      <img alt="Crates.io Total Downloads" src="https://img.shields.io/crates/d/gxf2bed">
    </a>
    <a href="https://anaconda.org/bioconda/gxf2bed" target="_blank">
      <img alt="Conda Platform" src="https://img.shields.io/conda/pn/bioconda/gxf2bed">
    </a>
  </p>


  <p align="center">
    The fastest G{F,T}F-to-BED converter around the block!

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
    before your eyes blink!
  </p>

</p>


check out the metrics below:

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `target/release/gxf2bed -i ./test/Homo_sapiens.GRCh38.112.gtf -o ./test/output.bed` | 1.592 ± 0.007 | 1.585 | 1.603 | 1.00 |
| `target/release/gxf2bed -i ./test/Homo_sapiens.GRCh38.112.gtf -o ./test/output.bed.gz` | 2.575 ± 0.033 | 2.551 | 2.631 | 1.62 ± 0.02 |
| `target/release/gxf2bed -i ./test/Homo_sapiens.GRCh38.112.gtf.gz -o ./test/output.bed` | 2.984 ± 0.011 | 2.965 | 2.993 | 1.87 ± 0.01 |
| `target/release/gxf2bed -i ./test/Homo_sapiens.GRCh38.112.gtf.gz -o ./test/output.bed.gz` | 3.946 ± 0.009 | 3.938 | 3.958 | 2.48 ± 0.01 |
| `target/release/gxf2bed -i ./test/gencode.v46.chr_patch_hapl_scaff.annotation.gff3.gz -o ./test/output.bed` | 2.046 ± 0.003 | 2.042 | 2.050 | 1.28 ± 0.01 |
| `target/release/gxf2bed -i ./test/gencode.v46.chr_patch_hapl_scaff.annotation.gff3.gz -o ./test/output.bed.gz` | 2.050 ± 0.014 | 2.038 | 2.072 | 1.29 ± 0.01 |
| `target/release/gxf2bed -i ./test/gencode.v46.chr_patch_hapl_scaff.annotation.gff3 -o ./test/output.bed` | 1.861 ± 0.014 | 1.847 | 1.879 | 1.17 ± 0.01 |
| `target/release/gxf2bed -i ./test/gencode.v46.chr_patch_hapl_scaff.annotation.gff3 -o ./test/output.bed.gz` | 2.914 ± 0.062 | 2.866 | 3.019 | 1.83 ± 0.04 |

Converts
- *Homo sapiens* GRCh38 Ensembl 112 GTF (254,129 transcripts) in 1.592 seconds.
- *Mus musculus* GRCm39 GENCODE 44 (149,547 transcritps) in 0.977 seconds.
- *Canis lupus familiaris* ROS_Cfam_1.0 Ensembl 112 (54,144 transcripts) in 0.501 seconds.
- *Gallus galus* bGalGal1 Ensembl 112 (72,689 transcripts) in 0.529 seconds.

> What's new on v.0.2.3
>
> - gxf2bed now is even x2 faster and can convert the human annotation in 1.5 seconds!
> - including bench to support results in README and update benchamrk results + metrics.
> - implements free compressed conversion between formats: gtf/gtf.gz/gff/gff.gz to bed/bed.gz!


## Usage
``` rust
Usage: gxf2bed[EXE] --input/-i <GTF/GFF> --output/-o <BED> [--parent/-p <PARENT>] [--child/-c <CHILD>] [--feature/-f <FEATURE>]

Arguments:
    --input/-i <GTF/GFF>: a .gtf/.gff file
    --output/-o <BED>: path to output .bed file
    --parent/-p <PARENT>: parent node [default: "transcript"]
    --child/-c <CHILD>: child node [default: "exon"]
    --feature/-f <FEATURE>: feature to extract from the attribute line [default: "transcript_id"]

Options:
    --help: print help
    --version: print version
    --threads/-t: number of threads (default: max ncpus)
```

> [!TIP]
> The interpretation of the `--parent/-p`, `--child/-c` and `--feature/-f` arguments is as follows:
> - `--parent/-p`: the parent node is the name of the record in the second column of the .gtf that will work as rule to extract the child nodes.
> - `--child/-c`: the child node is the name of the record in the second column of the .gtf that the tool will extract and build coordinates from
> - `--feature/-f`: the feature is the name of the attribute that will be extracted from the attribute line of the .gtf file in order to match parent and child
>
> The most common case is to use `--parent/-p "transcript" --child/-c "exon" --feature/-f "transcript_id"` to extract exons from transcripts, but the tool
> gives you the flexibility to extract any feature from any parent-child relationship in the .gtf file, like 3' UTRs, 5' UTRs, CDS, etc. For the latter,
> you can use `--parent/-p "transcript" --child/-c "three_prime_UTR" --feature/-f "trancript_id"` to extract 3'UTRs from genes, for example.

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
3. run `cargo run --release -- -i <GTF/GFF> -o <BED>`

## Container image
to build the development container image:
1. run `git clone https://github.com/alejandrogzi/gxf2bed.git && cd gxf2bed`
2. initialize docker with `start docker` or `systemctl start docker`
3. build the image `docker image build --tag gxf2bed .`
4. run `docker run --rm -v "[dir_where_your_gtf_is]:/dir" gxf2bed -i /dir/<GTF/GFF> -o /dir/<BED>`

## Conda
to use gxf2bed through Conda just:
1. `conda install gxf2bed -c bioconda` or `conda create -n gxf2bed -c bioconda gxf2bed`


## Benchmark + FAQ

| File | Format | Tool | Language | Time [s] | Fold-change | Size | n_transcripts |
|:---|:---:|:---|:---:|---:|---:|---:|---:|
| `Homo_sapiens.GRCh38.112.gtf` | GTF | gxf2bed | Rust | 1.592 | 1.00 | 1.4 GB | 254,129 |
|  |  | BEDOPS | C | 13.06 | 8.20 |  |  |
|  |  | ea | Perl | 15.69 | 9.86 |  |  |
|  |  | Signal & Brown | Python | 7.58 | 4.77 |  |  |
|  |  | UCSC's utils | C++ | 10.97 | 6.89 |  |  |
| `gencode.v46.chr_patch_hapl_scaff.annotation.gff3` | GFF | gxf2bed | Rust | 1.861 | 1.00 | 1.6 GB | 278,220 |
|  |  | BEDOPS | C | 18.24 | 9.80 |  |  |
|  |  | ea | Perl | 22.28 | 11.98 |  |  |
|  |  | Signal & Brown | Python | 47.02 | 25.24 |  |  |
|  |  | UCSC's utils | C++ | 17.42 | 9.35 |  |  |


If you google "gtf to bed" or "gff to bed" you'll find some posts about recommended tools. Part of these tools are already deprecated or too difficult to run (even in Linux), due to poor maintenance. This project was conceived to provide an easy way to convert GTF/GFF files to BED structures and finish a set of high-performance converters between gene model formats in Rust (bed2gtf, bed2gff and now gxf2bed).

To test the efficiency of gxf2bed, I took 4 GTF/GFF-to-BED converters and ran them with the current GRCh38 GENCODE annotation (v.44), that has ~250,000 transcripts [the values showed here are the mean of 5 consecutive runs]. Briefly, the competitors are:
- UCSC's utils: `gtfToGenePred | genePredToBed` & `gffToGenePred | genePredToBed` (taken from [here](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/))
- Signal & Brown: `gtf2bed.py` & `gff32gtf.py && gtf2bed.py` (taken from [here](https://github.com/signalbash/how_are_we_stranded_here/tree/master/how_are_we_stranded_here))
- ea-utils: `gtf2bed.pl` & `gff2gtf.pl | gtf2bed.pl` (taken from [here](https://github.com/ExpressionAnalysis/ea-utils/tree/master/clipper))
- BEDOPS: `gtf2bed` &  `gtf2bed` (taken from [here](https://github.com/bedops/bedops/tree/master/applications/bed/conversion/src/wrappers))

Besides the easy way to make `gxf2bed` run and the GFF-GTF channel that allows centralize both in a single converter, this tool showed a significant difference in computation time against the other tools in each one of the two formats. Even using a combination of different number of threads, times were practically the same 1.592s +/- 0.007s and the differences were maintained. Is important to note that GFF3 files took a lot more time than its GTF counterpart. This could be due to the need of chaining two different programs 1) GFF-to-GTF and 2) GTF-to-BED.

Taken together, `gxf2bed` offers an easier and faster way to convert GTF/GFF files to BED files. This tool could save you at least x2-3 times the time you used to spend using other tools to convert GTF files and x5 times if you want to convert GFF3 files!
