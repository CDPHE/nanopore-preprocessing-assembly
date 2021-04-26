version 1.0

workflow NanoporeGuppyAssembly {

    input {
        String    gcs_fastq_dir
        String    sample_id
        String    barcode
        File    covid_genome
        File    preprocess_python_script
    }
    call ListFastqFiles {
        input:
            gcs_fastq_dir = gcs_fastq_dir
    }
    call Demultiplex {
        input:
            fastq_files = ListFastqFiles.fastq_files,
            barcode = barcode
    }
    call Read_Filtering {
        input:
            fastq_files = Demultiplex.guppy_demux_fastq,
            barcode = barcode,
            sample_id = sample_id
    }
    call Medaka {
        input:
            filtered_reads = Read_Filtering.guppyplex_fastq,
            sample_id = sample_id,
            barcode = barcode
    }
    call Bam_stats {
        input:
            bam = Medaka.trim_sorted_bam,
            sample_id = sample_id,
            barcode = barcode
    }
    call Scaffold {
        input:
            sample_id = sample_id,
            barcode = barcode,
            ref = covid_genome,
            fasta = Medaka.consensus
    }

    call rename_fasta {
        input:
            sample_id = sample_id,
            fasta = Scaffold.scaffold_consensus
    }

    call calc_percent_cvg {
        input:
            sample_id = sample_id,
            fasta = rename_fasta.renamed_consensus,
            preprocess_python_script = preprocess_python_script
    }

    output {
        File barcode_summary = Demultiplex.barcode_summary
        Array[File] guppy_demux_fastq = Demultiplex.guppy_demux_fastq
        File filtered_fastq = Read_Filtering.guppyplex_fastq
        File sorted_bam = Medaka.sorted_bam
        File trim_sorted_bam = Medaka.trim_sorted_bam
        File trim_sorted_bai = Medaka.trim_sorted_bai
        File flagstat_out = Bam_stats.flagstat_out
        File samstats_out = Bam_stats.stats_out
        File covhist_out = Bam_stats.covhist_out
        File cov_out = Bam_stats.cov_out
        File variants = Medaka.variants
        File consensus = Medaka.consensus
        File scaffold_consensus = Scaffold.scaffold_consensus
        File renamed_consensus = rename_fasta.renamed_consensus
        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv
        String artic_pipeline_version = Medaka.artic_pipeline_version
    }
}

task ListFastqFiles {
    input {
        String gcs_fastq_dir
    }

    String indir = sub(gcs_fastq_dir, "/$", "")

    command <<<
        gsutil ls ~{indir}/**.fastq > fastq_files.txt
    >>>

    output {
        Array[File] fastq_files = read_lines("fastq_files.txt")
    }

    runtime {
        cpu:    1
        memory:    "1 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
}

task Demultiplex {
    input {
        Array[File] fastq_files
        String barcode
    }

    Int disk_size = 3 * ceil(size(fastq_files, "GB"))

    command <<<
        set -e
        mkdir fastq_files
        ln -s ~{sep=' ' fastq_files} fastq_files
        ls -alF fastq_files
        guppy_barcoder --require_barcodes_both_ends --barcode_kits "EXP-NBD196" --fastq_out -i fastq_files -s demux_fastq
        ls -alF demux_fastq
    >>>

    output {
        Array[File] guppy_demux_fastq = glob("demux_fastq/${barcode}/*.fastq")
        File barcode_summary = "demux_fastq/barcoding_summary.txt"
    }

    runtime {
        cpu:    8
        memory:    "16 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    30
        preemptible:    0
        maxRetries:    0
        gpuType:    "nvidia-tesla-p100"
        gpuCount:    1
        nvidiaDriverVersion:    "418.152.00"
        zones:    ["us-east1-c"]
        cpuPlatform:    "Intel Haswell"
        docker:    "genomicpariscentre/guppy"
    }
}

task Read_Filtering {
    input {
        Array[File] fastq_files 
        String barcode
        String sample_id
    }

    command <<<
        set -e
        mkdir fastq_files
        ln -s ~{sep=' ' fastq_files} fastq_files
        ls -alF fastq_files
        
        artic guppyplex --min-length 400 --max-length 700 --directory fastq_files --output ~{sample_id}_~{barcode}.fastq

    >>>

    output {
        File guppyplex_fastq = "${sample_id}_${barcode}.fastq"
    }

    runtime {
        docker: "theiagen/artic-ncov2019:1.1.3"
        memory: "16 GB"
        cpu: 8
        disks: "local-disk 100 SSD"
        preemptible: 0
    }
}

task Medaka {
    input {
        String barcode
        String sample_id
        File filtered_reads
    }

    command <<<
    
        artic -v > VERSION
        artic minion --medaka --normalise 20000 --threads 8 --scheme-directory /artic-ncov2019/primer_schemes --read-file ~{filtered_reads} nCoV-2019/V3 ~{barcode}
        
        cp ~{barcode}.consensus.fasta ~{sample_id}_~{barcode}.consensus.fasta
        cp ~{barcode}.trimmed.rg.sorted.bam ~{sample_id}_~{barcode}.trimmed.rg.sorted.bam
        cp ~{barcode}.primertrimmed.rg.sorted.bam ~{sample_id}_~{barcode}.primertrimmed.rg.sorted.bam
        cp ~{barcode}.primertrimmed.rg.sorted.bam.bai ~{sample_id}_~{barcode}.primertrimmed.rg.sorted.bam.bai
        cp ~{barcode}.pass.vcf.gz ~{sample_id}_~{barcode}.pass.vcf.gz

    >>>

    output {
        File consensus = "${sample_id}_${barcode}.consensus.fasta"
        File sorted_bam = "${sample_id}_${barcode}.trimmed.rg.sorted.bam"
        File trim_sorted_bam = "${sample_id}_${barcode}.primertrimmed.rg.sorted.bam"
        File trim_sorted_bai = "${sample_id}_${barcode}.primertrimmed.rg.sorted.bam.bai"
        File variants = "${sample_id}_${barcode}.pass.vcf.gz"
        String artic_pipeline_version = read_string("VERSION")
    }

    runtime {
        docker: "theiagen/artic-ncov2019:1.1.3"
        memory: "16 GB"
        cpu: 8
        disks: "local-disk 100 SSD"
        preemptible: 0
    }
}

task Bam_stats {
    input {
        String sample_id
        String barcode
        File bam
    }

    Int disk_size = 3 * ceil(size(bam, "GB"))

    command {

        samtools flagstat ${bam} > ${sample_id}_${barcode}_flagstat.txt

        samtools stats ${bam} > ${sample_id}_${barcode}_stats.txt

        samtools coverage -m -o ${sample_id}_${barcode}_coverage_hist.txt ${bam}

        samtools coverage -o ${sample_id}_${barcode}_coverage.txt ${bam}

    }

    output {

        File flagstat_out  = "${sample_id}_${barcode}_flagstat.txt"
        File stats_out  = "${sample_id}_${barcode}_stats.txt"
        File covhist_out  = "${sample_id}_${barcode}_coverage_hist.txt"
        File cov_out  = "${sample_id}_${barcode}_coverage.txt"
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/samtools:1.10"
    }
}

task Scaffold {
    input {
        String sample_id
        String barcode
        File fasta
        File ref
    }

    Int disk_size = 3 * ceil(size(fasta, "GB"))

    command {

        pyScaf.py -f ${fasta} -o ${sample_id}_${barcode}_consensus_scaffold.fa -r ${ref}

    }

    output {
        File scaffold_consensus = "${sample_id}_${barcode}_consensus_scaffold.fa"
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "chrishah/pyscaf-docker"
    }
}

task rename_fasta {

    input {

        String sample_id
        File fasta
    }

    command {

        sed 's/>.*/>CO-CDPHE-~{sample_id}/' ~{fasta} > ~{sample_id}_consensus_renamed.fa

    }

    output {

        File renamed_consensus  = "${sample_id}_consensus_renamed.fa"

    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task calc_percent_cvg {

    input {

        File fasta
        String sample_id
        File preprocess_python_script
    }

    command {
        python ~{preprocess_python_script} \
          --sample_id ~{sample_id} \
          --fasta_file ~{fasta}
      }

    output {

      File percent_cvg_csv  = "${sample_id}_consensus_cvg_stats.csv"

    }

    runtime {

      docker: "mchether/py3-bio:v4"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }

}
