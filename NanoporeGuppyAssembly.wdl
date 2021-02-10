version 1.0

workflow NanoporeGuppyAssembly {

    input {
        String    gcs_fastq_dir
        String    sample_id
        String    barcode
        File    covid_genome
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
    call MergeFastq {
        input:
            fastq_files = Demultiplex.guppy_demux_fastq
    }
    call FiltLong {
        input:
            merged_fastq = MergeFastq.merged_fastq,
            sample_id = sample_id,
            barcode = barcode
    }
    call Minimap2 {
        input:
            reads = FiltLong.filtered_fastq,
            ref = covid_genome,
            sample_id = sample_id,
            barcode = barcode
    }
    call Bam_stats {
        input:
            bam = Minimap2.sorted_bam,
            sample_id = sample_id,
            barcode = barcode
    }
    call Variants {
        input:
            bam = Minimap2.sorted_bam,
            ref = covid_genome,
            sample_id = sample_id,
            barcode = barcode
    }
    call Consensus {
        input:
            ref = covid_genome,
            sample_id = sample_id,
            barcode = barcode,
            reads = FiltLong.filtered_fastq
    }
    call Pangolin {
        input :
            sample_id = sample_id,
            barcode = barcode,
            fasta = Consensus.consensus
    }
    output {
        File barcode_summary = Demultiplex.barcode_summary
        Array[File] guppy_demux_fastq = Demultiplex.guppy_demux_fastq
        File merged_fastq = MergeFastq.merged_fastq
        File filtered_fastq = FiltLong.filtered_fastq
        File sorted_bam = Minimap2.sorted_bam
        File flagstat_out = Bam_stats.flagstat_out
        File samstats_out = Bam_stats.stats_out
        File covhist_out = Bam_stats.covhist_out
        File cov_out = Bam_stats.cov_out
        File variants = Variants.vcf_final
        File consensus = Consensus.consensus
        File pangolin = Pangolin.lineage
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
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    30
        preemptible:    0
        maxRetries:    0
        gpuType:    "nvidia-tesla-p100"
        gpuCount:    1
        nvidiaDriverVersion:    "418.152.00"
        zones:    ["us-east1-c"]
        cpuPlatform:    "Intel Haswell"
        docker:    "genomicpariscentre/guppy:4.4.1"
    }
}

task MergeFastq {
    input {
        Array[File] fastq_files
    }

    Int disk_size = 3 * ceil(size(fastq_files, "GB"))

    command <<<

        cat ~{sep=" " fastq_files} > merged.fastq

    >>>

    output {
        File merged_fastq = "merged.fastq"
    }

    runtime {
        cpu:    1
        memory:    "4 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
}

task FiltLong {
    input {
        String sample_id
        String barcode
        File merged_fastq
    }

    Int disk_size = 3 * ceil(size(merged_fastq, "GB"))

    command <<<
    
        filtlong --min_length 400 ~{merged_fastq} | gzip > ~{sample_id}_~{barcode}_filtered_reads.fastq.gz
    
    >>>

    output {
        File filtered_fastq = "${sample_id}_${barcode}_filtered_reads.fastq.gz"
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/filtlong:0.2.0-cv1"
    }
}

task Minimap2 {
    input {
        String sample_id
        String barcode
        File reads
        File ref
    }

    Int disk_size = 3 * ceil(size(reads, "GB"))

    command <<<

        minimap2 -a -x map-ont ~{ref} ~{reads} > ~{sample_id}_~{barcode}_minimap2.sam
        
        samtools sort -O bam -o ~{sample_id}_~{barcode}_sorted.bam ~{sample_id}_~{barcode}_minimap2.sam
    
    >>>

    output {
        File sam = "${sample_id}_${barcode}_minimap2.sam"
        File sorted_bam = "${sample_id}_${barcode}_sorted.bam"
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "matmu/nanopore"
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

task Variants {
    input {
        String sample_id
        String barcode
        File bam
        File ref
    }

    Int disk_size = 3 * ceil(size(bam, "GB"))

    command <<<
    
        samtools index ~{bam}
        
        medaka_variant -i ~{bam} -f ~{ref}
        
        cp medaka_variant/round_1.vcf medaka_variant/~{sample_id}_~{barcode}.vcf
    
    >>>

    output {
        File vcf_final = "medaka_variant/${sample_id}_${barcode}.vcf"
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/artic-ncov2019-medaka:1.1.0"
    }
}

task Consensus {
    input {
        String sample_id
        String barcode
        File ref
        File reads
    }

    Int disk_size = 3 * ceil(size(reads, "GB"))

    command {
    
        medaka_consensus -i ${reads} -d ${ref}
        
        cp medaka/consensus.fasta medaka/${sample_id}_${barcode}_consensus.fa
    
    }

    output {
        File consensus = "medaka/${sample_id}_${barcode}_consensus.fa"
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/artic-ncov2019-medaka:1.1.0"
    }
}

task Pangolin {
    input {
        String sample_id
        String barcode
        File fasta
    }

    Int disk_size = 3 * ceil(size(fasta, "GB"))

    command {
    
        pangolin ${fasta} --outfile ${sample_id}_${barcode}_lineage_report.csv
    
    }

    output {
        File lineage = "${sample_id}_${barcode}_lineage_report.csv"
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/pangolin"
    }
}