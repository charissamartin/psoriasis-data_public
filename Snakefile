#############################################
# Snakefile for analyzing psoriasis RNA-Seq #
# Charissa Martin                           #
# 2020-01-22                                #
#############################################

configfile: "snakemake_config.yaml"

########################
# Function definitions #
########################
def makeURL( SRRnumber, base_url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/' ):
    pattern1 = SRRnumber[:6]
    pattern2 = '00'+ SRRnumber[-1]
    link = base_url + pattern1 + '/' + pattern2 + '/' + SRRnumber
    return link

#######################################
# Read SRR identifiers from data file #
#######################################
with open( "info/full_sample_meta_data.txt", 'r' ) as file:
    SRR_identifer_list = [line.split( '\t' )[0] for line in file]
    SRR_identifer_list = SRR_identifer_list[1:]

###############################
# download rule for snakemake #
###############################
url_dict = {}
for SRR in SRR_identifer_list:
    url_dict[SRR] = makeURL(SRR)

##########################
# all rule for snakemake #
##########################
rule all:
    input:
        expand("output/raw_fastqc/{identifier}_{number}_fastqc.html", identifier = SRR_identifer_list, number = [ 1,2 ]),
        expand("output/trimmed_fastqc/{identifier}_{number}_clean_fastqc.html", identifier = SRR_identifer_list, number = [ 1,2 ]),
        expand("output/STAR/{identifier}/{star_output_files}", identifier = SRR_identifer_list, star_output_files = ['Aligned.out.bam', 'Log.final.out', 'Log.out', 'Log.progress.out', 'SJ.out.tab.gz', 'Unmapped.out.mate1.gz', 'Unmapped.out.mate2.gz' ] ),
        expand("output/featureCounts/{identifier}.counts.txt", identifier = SRR_identifer_list ),
        expand("output/salmon/{identifier}/{salmon_output_files}", identifier = SRR_identifer_list, salmon_output_files = [ 'quant.sf', 'quant.genes.sf', 'cmd_info.json', 'lib_format_counts.json', 'abundance.h5' ] )

##########################
# implicit download rule #
##########################
rule download:
    params:
        url = lambda wildcards: expand( "{url}/{identifier}_{number}.fastq.gz",
        url = url_dict[wildcards.identifier],
        identifier = wildcards.identifier,
        number = wildcards.number )
    output:
        "data/fastq_raw/{identifier}_{number}.fastq.gz"
    shell:
        "wget -O {output} {params.url}"

############################
# generic compression rule #
############################
rule pigz:
    input:
        "{file}"
    output:
        "{file}.gz"
    threads: 4
    shell:
        "pigz --best --processes {threads} {input}"

#############################
# fastqc rule for snakemake #
#############################
rule fastqc:
    input:
        "data/fastq_raw/{identifier}_{number}.fastq.gz"
    output:
        "output/raw_fastqc/{identifier}_{number}_fastqc.html",
        temp( "output/raw_fastqc/{identifier}_{number}_fastqc.zip" )
    params:
        fastqc_path = config[ "fastqc_path" ]
    threads: 1
    shell:
        "{params.fastqc_path} -q -f fastq -o output/raw_fastqc --threads {threads} {input}"

###############################
# cutadapt rule for snakemake #
###############################
rule cutadapt:
    input:
        r1 = "data/fastq_raw/{identifier}_1.fastq.gz",
        r2 = "data/fastq_raw/{identifier}_2.fastq.gz"
    output:
        r1 = "output/fastq_clean/{identifier}_1_clean.fastq.gz",
        r2 = "output/fastq_clean/{identifier}_2_clean.fastq.gz"
    log:
        "output/fastq_clean/{identifier}.log"
    shell:
        "cutadapt -m 15 -q 10 -O 5 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o {output.r1} -p {output.r2} {input.r1} {input.r2} 1>{log} 2>&1"

#######################################
# post trim fastqc rule for snakemake #
#######################################
rule trimmed_fastqc:
    input:
        "output/fastq_clean/{identifier}_{number}_clean.fastq.gz"
    output:
        "output/trimmed_fastqc/{identifier}_{number}_clean_fastqc.html",
        temp( "output/trimmed_fastqc/{identifier}_{number}_clean_fastqc.zip" )
    params:
        fastqc_path = config[ "fastqc_path" ]
    threads: 1
    shell:
        "{params.fastqc_path} -q -f fastq -o output/trimmed_fastqc --threads {threads} {input}"

###############################
# STAR rule for snakemake #
###############################
rule star:
    input:
        r1 = "output/fastq_clean/{identifier}_1_clean.fastq.gz",
        r2 = "output/fastq_clean/{identifier}_2_clean.fastq.gz"
    output:
        "output/STAR/{identifier}/Aligned.out.bam",
        "output/STAR/{identifier}/Log.final.out",
        "output/STAR/{identifier}/Log.out",
        "output/STAR/{identifier}/Log.progress.out",
        "output/STAR/{identifier}/SJ.out.tab",
        "output/STAR/{identifier}/Unmapped.out.mate1",
        "output/STAR/{identifier}/Unmapped.out.mate2"
    threads: 8
    params:
        output_dir = lambda wildcards: "output/STAR/%s/" % ( wildcards.identifier ),
        star_reference = config[ 'star_index_path' ]
    shell:
        "mkdir -p {params.output_dir} && "
        "STAR "
        "--runThreadN {threads} "
        "--genomeDir {params.star_reference} "
        "--readFilesCommand zcat "
        "--outSAMtype BAM Unsorted "
        "--outFileNamePrefix {params.output_dir} "
        "--outReadsUnmapped Fastx "
        "--readFilesIn {input}"

####################################
# featureCounts rule for snakemake #
####################################
rule featureCounts:
    input:
        "output/STAR/{identifier}/Aligned.out.bam"
    output:
        "output/featureCounts/{identifier}.counts.txt"
    log:
        "output/featureCounts/{identifier}.log"
    params:
        saf_file = config[ 'saf_format_genes' ],
        strand_string = config[ 'feature_count_strand']
    threads: 4
    shell:
        "featureCounts "
        "-a {params.saf_file} "
        "-p "
        "-F SAF "
        "-Q 5 "
        "-T {threads} "
        "-s {params.strand_string} "
        "-o {output} "
        "{input}"

#############################
# salmon rule for snakemake #
#############################
rule salmon:
    input:
        r1 = "output/fastq_clean/{identifier}_1_clean.fastq.gz",
        r2 = "output/fastq_clean/{identifier}_2_clean.fastq.gz"
    output:
        "output/salmon/{identifier}/quant.sf",
        "output/salmon/{identifier}/quant.genes.sf",
        "output/salmon/{identifier}/cmd_info.json",
        "output/salmon/{identifier}/lib_format_counts.json"
    log:
        "output/salmon/{identifier}/logs"
    params:
        output_dir = lambda wildcards: "output/salmon/%s" % ( wildcards.identifier ),
        salmon_index = config[ 'salmon_index_path' ],
        gene_map = config [ 'salmon_gene_transcript_map' ],
        num_bootstraps = "%s" % ( config[ 'salmon_bootstraps' ] )
    threads: 4
    shell:
        "salmon quant "
        "--libType ISF "
        "-1 {input.r1} "
        "-2 {input.r2} "
        "--index {params.salmon_index} "
        "--seqBias "
        "--gcBias "
        "--validateMappings "
        "--geneMap {params.gene_map} "
        "--threads {threads} "
        "--numBootstraps {params.num_bootstraps} "
        "--output {params.output_dir}"

############################
# wasabi of salmon results #
############################
rule wasabi:
    input:
        "output/salmon/{identifier}/quant.sf"
    output:
        "output/salmon/{identifier}/abundance.h5"
    log:
        "logs/wasabi_{identifier}.log"
    script:
        "snake_scripts/run_wasabi.R"
