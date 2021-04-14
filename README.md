# nanopore-preprocessing-assembly

## preprocessing and assembly workflows for nanopore data

### NanoporeGuppyAssembly.wdl

The NanoporeGuppyAssembly.wdl workflow was developed for the preprocessing and assembly of Oxford Nanopore sequencing of SARS-CoV-2 to be run on the GCP Terra platform. It takes basecalled fastq files as input (e.g. from MinKNOW).

The columns for the input data table in Terra should be arranged as:
1. entity:sample_id (column of sample names/ids). If there is more than one data table in the Terra Workspace, add a number after the word sample (e.g. entity:sample2_id).
2. barcode (the ONT barocde for the sample, e.g. barcode01)
3. fastq_dir (the path to the google bucket directory containing the basecalled fastq files)

Also needed as input: covid_genome, save the path to the google bucket directory containing the SARS-CoV-2 reference genome as Workspace data and include in Terra input field for the workflow.

The NanoporeGuppyAssembly.wdl workflow will:
1. Demultiplex bascalled fastq files using guppy_barcoder
2. Merge the demultiplexed fastq files into a single merged fastq
3. Perform quality filtering on the fastq files using FiltLong
4. Align reads to the SARS-CoV-2 reference genome using Minimap2
5. Generate base and alignment quality and coverage statistics on the bam file using Samtools
6. Call variants on the bam file and generate a vcf using Medaka
7. Generate a draft assembly with Medaka
8. Scaffold assembly with pyScaf 
9. Rename consensus to CO-CDPHE-{sample_id}
10. Calculates percent coverage

External tools used in this workflow were from publicly available Docker images:
1. General utilies docker images: us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6, theiagen/utility:1.0, and mchether/py3-bio:v4
2. Guppy: Oxford Nanopore Technologies Guppy, https://community.nanoporetech.com/
  docker image: genomicpariscentre/guppy:4.4.1
4. FiltLong: https://github.com/rrwick/Filtlong
  docker image: staphb/filtlong:0.2.0-cv1
6. Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191
  docker image: matmu/nanopore
8. Samtools: Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, Twelve years of SAMtools and BCFtools, GigaScience (2021) 10(2) giab008 [33590861]
  docker image: staphb/samtools:1.10
10. Medaka: https://nanoporetech.github.io/medaka/installation.html
  docker image: staphb/artic-ncov2019-medaka:1.1.0
12. PyScaf: https://github.com/lpryszcz/pyScaf
  docker image: chrishah/pyscaf-docker
