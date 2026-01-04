<p align="center">
  <p align="center">
    <img width=200 align="center" src="./supp/logo.png" >
  </p>

  <span>
    <h1 align="center">
        gxf2bed
    </h1>
  </span>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.3.0-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.3.0-green">
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
  </p>

  <p align="center">
    <samp>
        <span> versatile BED-to-G{T,F}F converter in Rust</span>
        <br>
        <br>
        <a href="https://docs.rs/gxf2bed/0.3.0/gxf2bed/">docs</a> .
        <a href="https://github.com/alejandrogzi/gxf2bed?tab=readme-ov-file#Usage">usage</a> .
        <a href="https://github.com/alejandrogzi/gxf2bed?tab=readme-ov-file#Features">features</a> .
        <a href="https://github.com/alejandrogzi/gxf2bed/blob/master/assets/EXAMPLES.md">examples</a>
    </samp>
  </p>

</p>


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


> What's new on v.0.3
>
> - Converts GTF and GFF3 files to BED format
> - Supports multiple BED formats (BED3, BED4, BED5, BED6, BED9, BED12)
> - Handles compressed input files (gzip, zstd, bzip2)
> - Customizable parent/child feature and attribute mapping
> - Optional additional fields from GTF/GFF attributes
> - Memory-mapped I/O for uncompressed files
> - Uses [genepred](https://github.com/alejandrogzi/genepred) for parsing and writing


## Usage
 ```bash
 gxf2bed -i <INPUT> -o <OUTPUT> [OPTIONS]

 Required arguments:
   -i, --input <GXF>          Path to GTF/GFF file
   -o, --output <BED>         Path to output BED file

 Optional arguments:
   -T, --threads <THREADS>    Number of threads (default: CPU count)
   -F, --parent-feature <PARENT>     Parent feature
   -f, --child-features <CHILDS>    Child features (comma-separated)
   -A, --parent-attribute <FEATURE>  Feature to extract
   -a, --child-attribute <CHILD>     Child feature to extract
   -t, --type <BED_TYPE>      BED type format (3, 4, 5, 6, 9, 12) [default: 12]
   -d, --additional-fields <ADDITIONAL>  BED additional fields (comma-separated)
   -c, --chunks <CHUNKS>      Chunk size for parallel processing [default: 15000]
   -h, --help                 Print help
   -V, --version              Print version
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

> [!WARNING]
> The data in this benchmark corresponds to version 0.2.* and below. Versions >0.3 are 
> a bit slower than the previous versions due to algorithmic changes.
## Benchmark + FAQ
<p align="center">
    <img width=700 align="center" src="./supp/gxf2bed.jpg">
</p>

If you google "gtf to bed" or "gff to bed" you'll find some posts about recommended tools. Part of these tools are already deprecated or too difficult to run (even in Linux), due to poor maintenance. This project was conceived to provide an easy way to convert GTF/GFF files to BED structures and finish a set of high-performance converters between gene model formats in Rust (bed2gtf, bed2gff and now gxf2bed). 

To test the efficiency of gxf2bed, I took 4 GTF/GFF-to-BED converters and ran them with the current GRCh38 GENCODE annotation (v.44), that has ~250,000 transcripts [the values showed here are the mean of 5 consecutive runs]. Briefly, the competitors are:
- UCSC's utils: `gtfToGenePred | genePredToBed` & `gffToGenePred | genePredToBed` (taken from [here](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/))
- Signal & Brown: `gtf2bed.py` & `gff32gtf.py && gtf2bed.py` (taken from [here](https://github.com/signalbash/how_are_we_stranded_here/tree/master/how_are_we_stranded_here))
- ea-utils: `gtf2bed.pl` & `gff2gtf.pl | gtf2bed.pl` (taken from [here](https://github.com/ExpressionAnalysis/ea-utils/tree/master/clipper))
- BEDOPS: `gtf2bed` &  `gtf2bed` (taken from [here](https://github.com/bedops/bedops/tree/master/applications/bed/conversion/src/wrappers)) 

Besides the easy way to make `gxf2bed` run and the GFF-GTF channel that allows centralize both in a single converter, this tool showed a significant difference in computation time against the other tools in each one of the two formats. Even using a combination of different number of threads, times were practically the same 3.5s +/- 0.4s and the differences were maintained. Is important to note that GFF3 files took a lot more time than its GTF counterpart. This could be due to the need of chaining two different programs 1) GFF-to-GTF and 2) GTF-to-BED.

Taken together, `gxf2bed` offers an easier and faster way to convert GTF/GFF files to BED files. This tool could save you at least x2-3 times the time you used to spend using other tools to convert GTF files and x5 times if you want to convert GFF3 files!


