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
2. Perform quality filtering using guppyplex
4. Align reads to the SARS-CoV-2 reference genome using Minimap2
5. Run artic minion --medaka for variat calling and to generate a consensus fasta
8. Scaffold assembly with pyScaf 
9. Rename consensus to CO-CDPHE-{sample_id}
10. Generate bam quality statstics using Samtools
11. Calculates percent coverage

External tools used in this workflow were from publicly available Docker images:

1. General utilies docker images: us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6, theiagen/utility:1.0, and mchether/py3-bio:v4
2. Guppy: Oxford Nanopore Technologies Guppy, https://community.nanoporetech.com/ docker image: genomicpariscentre/guppy:4.4.1
3. Samtools: Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, Twelve years of SAMtools and BCFtools, GigaScience (2021) 10(2) giab008 [33590861] docker image: staphb/samtools:1.10
4. Medaka: https://nanoporetech.github.io/medaka/installation.html docker image: staphb/artic-ncov2019-medaka:1.1.0
5. PyScaf: https://github.com/lpryszcz/pyScaf docker image: chrishah/pyscaf-docker
6. This workflow follows the SARS-CoV-2 nanopore whole genome assembly guidelines from Artic: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html using the docker image from theiagen/artic-ncov2019:1.1.3
