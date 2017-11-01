shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")  # bash "safe mode"

samples = [x.rsplit()[0] for x in open("./RawData/accession_sra.tsv", "r")]

rule all:
    """ Map all """
    input:
        expand("ProcessedData/{sample}/{sample}_kbs_npq.bg", sample=samples),

rule download:
    """ download data from SRA """
    output:
        "ProcessedData/{sample}/{sample}_1.fastq.gz"
    shell:
        """
        ./code/download_from_sra.sh -f ./RawData/accession_sra.tsv -p 10 -q
        python3 ./code/move_sra_files.py
        """

rule trim:
    """ trim reads """
    input:
        read1="ProcessedData/{sample}/{sample}_1.fastq.gz",
        read2="ProcessedData/{sample}/{sample}_2.fastq.gz"
    output:
        read1="ProcessedData/{sample}/{sample}_1_trim.fq",
        read2="ProcessedData/{sample}/{sample}_2_trim.fq",
        se1="ProcessedData/{sample}/{sample}_1_se.fq",
        se2="ProcessedData/{sample}/{sample}_2_se.fq"
    threads: 10
    shell:
        """
        /usr/bin/java -jar /home/san/tstuart/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
          -threads {threads} \
          {input.read1} {input.read2} \
          {output.read1} {output.se1} {output.read2} {output.se2} \
          ILLUMINACLIP:/home/san/tstuart/tools/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule compress:
    """ gzip reads """
    input:
        "ProcessedData/{sample}/{sample}_1_trim.fq",
        "ProcessedData/{sample}/{sample}_2_trim.fq",
        "ProcessedData/{sample}/{sample}_1_se.fq",
        "ProcessedData/{sample}/{sample}_2_se.fq"
    output:
        "ProcessedData/{sample}/{sample}_1_trim.fq.gz",
        "ProcessedData/{sample}/{sample}_2_trim.fq.gz",
        "ProcessedData/{sample}/{sample}_1_se.fq.gz",
        "ProcessedData/{sample}/{sample}_2_se.fq.gz"
    threads: 20
    shell:
        """
        pigz -p {threads} {input}
        """

rule qc:
    """ run fastqc """
    input:
        "ProcessedData/{sample}/{sample}_1_trim.fq.gz",
        "ProcessedData/{sample}/{sample}_2_trim.fq.gz",
        "ProcessedData/{sample}/{sample}_1_se.fq.gz",
        "ProcessedData/{sample}/{sample}_2_se.fq.gz"
    shell:
        """
        fastqc -t 4 {input}
        """

rule map:
    """ align reads """
    input:
        read1="ProcessedData/{sample}/{sample}_1_trim.fq.gz",
        read2="ProcessedData/{sample}/{sample}_2_trim.fq.gz"
    output:
        kbs="ProcessedData/{sample}/{sample}_unsort_kbs.bam",
        tair10="ProcessedData/{sample}/{sample}_unsort_tair10.bam"
    threads: 20
    shell:
        """
        bowtie2 --local --dovetail -p {threads} --fr -q -R5 -N1 -x /home/san/tstuart/github/npq-qtl/RawData/kbs -X 200 \
           -1 {input.read1} -2 {input.read2} \
          | samtools view -b - > {output.kbs}

        bowtie2 --local --dovetail -p {threads} --fr -q -R5 -N1 -x /home/san/tstuart/github/npq-qtl/RawData/tair10 -X 200 \
           -1 {input.read1} -2 {input.read2} \
          | samtools view -b - > {output.tair10}
        """

rule sort:
    """ sort reads """
    input:
        kbs="ProcessedData/{sample}/{sample}_unsort_kbs.bam",
        tair10="ProcessedData/{sample}/{sample}_unsort_tair10.bam"
    output:
        kbs="ProcessedData/{sample}/{sample}_kbs.bam",
        tair10="ProcessedData/{sample}/{sample}_tair10.bam"
    threads: 20
    shell:
        """
        samtools sort -@ {threads} -O bam {input.kbs} > {output.kbs}
        samtools sort -@ {threads} -O bam {input.tair10} > {output.tair10}
        """

rule cleanup:
    """ remove unsorted bams if sorted bam present """
    input:
        "ProcessedData/{sample}/{sample}_kbs.bam",
        "ProcessedData/{sample}/{sample}_tair10.bam",
        kbs="ProcessedData/{sample}/{sample}_unsort_kbs.bam",
        tair10="ProcessedData/{sample}/{sample}_unsort_tair10.bam"
    shell:
        """
        rm {input.kbs} {input.tair10}
        """

rule index:
    """ index bam files """
    input:
        kbs="ProcessedData/{sample}/{sample}_kbs.bam",
        tair10="ProcessedData/{sample}/{sample}_tair10.bam"
    output:
        kbs="ProcessedData/{sample}/{sample}_kbs.bam.bai",
        tair10="ProcessedData/{sample}/{sample}_tair10.bam.bai"
    shell:
        """
        samtools index {input.kbs}
        samtools index {input.tair10}
        """

rule coverage:
    """ compute coverage over NPQ region """
    input:
        kbs="ProcessedData/{sample}/{sample}_kbs.bam",
        tair10="ProcessedData/{sample}/{sample}_tair10.bam"
    output:
        kbs="ProcessedData/{sample}/{sample}_kbs.bg",
        tair10="ProcessedData/{sample}/{sample}_tair10.bg"
    threads: 20
    shell:
        """
        bamCoverage -b {input.kbs} -o {output.kbs} -p {threads} --normalizeUsingRPKM -of bedgraph
        bamCoverage -b {input.tair10} -o {output.tair10} -p {threads} --normalizeUsingRPKM -of bedgraph
        """

rule region:
    """ extract the NPQ region """
    input:
        "ProcessedData/{sample}/{sample}_kbs.bam.bai",
        "ProcessedData/{sample}/{sample}_tair10.bam.bai",
        kbs_bg="ProcessedData/{sample}/{sample}_kbs.bg",
        tair10_bg="ProcessedData/{sample}/{sample}_tair10.bg",
        kbs_bam="ProcessedData/{sample}/{sample}_kbs.bam",
        tair10_bam="ProcessedData/{sample}/{sample}_tair10.bam"
    output:
        kbs_bam="ProcessedData/{sample}/{sample}_kbs_npq.bam",
        tair10_bam="ProcessedData/{sample}/{sample}_tair10_npq.bam",
        kbs_bg="ProcessedData/{sample}/{sample}_kbs_npq.bg",
        tair10_bg="ProcessedData/{sample}/{sample}_tair10_npq.bg"
    shell:
        """
        samtools view -b {input.kbs_bam} contig30:900000-950000 > {output.kbs_bam}
        samtools view -b {input.tair10_bam} 1:16851824-16891823 > {output.tair10_bam}
        samtools index {output.kbs_bam}
        samtools index {output.tair10_bam}

        bedtools intersect -a {input.kbs_bg} -b RawData/npq_kbs.bed > {output.kbs_bg}
        bedtools intersect -a {input.tair10_bg} -b RawData/npq_tair10.bed > {output.tair10_bg}
        """
