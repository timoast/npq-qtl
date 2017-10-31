shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")  # bash "safe mode"

samples = [x.rsplit()[0] for x in open("./RawData/accession_sra.tsv", "r")]

rule all:
    """ Map all """
    input:
        expand("ProcessedData/{sample}/{sample}_kbs_npq.bw", sample=samples),

rule download:
    """ download data from SRA """
    output:
        "ProcessedData/{sample}/{sample}_1.fastq.gz"
    shell:
        """
        ./code/download_from_sra.sh -f ./RawData/accession_sra.tsv -p 10 -q
        python3 ./code/rename_files.py
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
    threads = 10
    shell:
        """
        /usr/bin/java -jar /home/san/tstuart/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
          -threads {threads} \
          {input.read1} {input.read2} \
          {output.read1} {output.se1} {output.read1} {output.se2} \
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
        tai10="ProcessedData/{sample}/{sample}_unsort_tair10.bam"
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
        tai10="ProcessedData/{sample}/{sample}_unsort_tair10.bam"
    output:
        kbs="ProcessedData/{sample}/{sample}_kbs.bam",
        tai10="ProcessedData/{sample}/{sample}_tair10.bam"
    threads: 20
    shell:
        """
        samtools sort -@ {threads} -T temp -O bam {input.kbs} -o {output.kbs}
        samtools sort -@ {threads} -T temp -O bam {input.tair10} -o {output.tair10}
        """

rule region:
    """ extract the NPQ region """
    input:
        region_kbs="RawData/npq_kbs.bed",
        region_tair10="RawData/npq_tair10.bed",
        kbs="ProcessedData/{sample}/{sample}_kbs.bam",
        tai10="ProcessedData/{sample}/{sample}_tair10.bam"
    output:
        kbs="ProcessedData/{sample}/{sample}_kbs_npq.bam",
        tai10="ProcessedData/{sample}/{sample}_tair10_npq.bam"
    shell:
        """
        samtools view -b {input.kbs} {input.region_kbs} > {output.kbs}
        samtools view -b {input.tair10} {input.region_tair10} > {output.tair10}
        """

rule coverage:
    """ compute coverage over NPQ region """
    input:
        kbs="ProcessedData/{sample}/{sample}_kbs_npq.bam",
        tai10="ProcessedData/{sample}/{sample}_tair10_npq.bam"
    output:
        kbs="ProcessedData/{sample}/{sample}_kbs_npq.bw",
        tai10="ProcessedData/{sample}/{sample}_tair10_npq.bw"
    threads: 20
    shell:
        """
        bamCoverage -b {input.kbs} -o {output.kbs} -p {threads}
        bamCoverage -b {input.tair10} -o {output.tair10} -p {threads}
        """
