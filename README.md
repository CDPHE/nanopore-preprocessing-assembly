# nanopore-preprocessing-assembly

## preprocessing and assembly workflows for nanopore data

### NanoporeGuppyAssembly.wdl

The NanoporeGuppyAssembly.wdl workflow was developed for the preprocessing and assembly of Oxford Nanopore sequencing of SARS-CoV-2 ro be run on the GCP Terra platform. It takes basecalled fastq files as input (e.g. from MinKNOW).

The columns for the input data table in Terra should be arranged as:
1. entity:sample_id (column of sample names/ids). If there is more than one data table in the Terra Workspave, add a number after the word sample (e.g. entity:sample2_id).
2. barcode (the ONT barocde for the sample, e.g. barcode01)
3. fastq_dir (the path to the google bucket directory containing the basecalled fastq files)
4. out_dir (google bucket path where you would like the output files transferred)

Also needed as input: covid_genome, save the path to the google bucket directory containing the SARS-CoV-2 reference genome as Workspace data and include in Terra input field for the workflow.

The NanoporeGuppyAssembly.wdl workflow will:
1. Demultiplex bascalled fastq files using guppy_barcoder
2. Merge the demultiplexed fastq files into a single merged fastq
3. Perform quality filtering on the fastq files using FiltLong
4. Align reads to the SARS-CoV-2 reference genome using Minimap2
5. Generate base and alignment quality and coverage statistics on the bam file using Samtools
6. Call variants on the bam file and generate a vcf using iVar
7. Generate a consensus genome fasta using iVar consensus
8. Assign SARS-CoV-2 lineages using Pangolin
9. Assign clades with Nextclade
10. Tranfer outputs to your chosen google bucket
