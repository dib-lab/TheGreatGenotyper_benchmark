configfile: "config.json"

input_reference=config['data']['reference']
reference_version=config['data']['reference_version']
outname=config['parameters']['outname_reads']

prefix = "" if "37" in reference_version or "19" in reference_version else "chr"


##############################################################
##################    prepare read data     ##################
##############################################################

#### compute coverage of full data ####
rule compute_bam_coverage:
	input:
		full_cov_bam='{results}/{sample}/aligned/{sample}-full.bam',
		full_cov_bai='{results}/{sample}/aligned/{sample}-full.bam.bai'
	output:
		"{results}/{sample}/raw/{sample}-coverage.cov"
	conda:
		'../env/genotyping.yml'
	resources:
		mem_mb=10000,
		runtime_hrs=5,
		runtime_min=1,
		runtime=5*60
	log:
		"{results}/{sample}/raw/{sample}-coverage.log"
	shell:
		"bash ../scripts/compute-coverage.sh {input.full_cov_bam} {output} &> {log}"
		

#### downsample reads ####
# rule downsample_reads:
# 	input:
# 		reads1=lambda wildcards: config['data'][wildcards.sample]['reads'][0],
# 		reads2=lambda wildcards: config['data'][wildcards.sample]['reads'][1],
# 		coverage="{results}/{sample}/raw/{sample}-coverage.cov"
# 	output:
# 		sampled1=("{results}/{sample}/raw/{sample}-{fraction, [0-9.]+}_1.fastq"),
# 		sampled2=("{results}/{sample}/raw/{sample}-{fraction, [0-9.]+}_2.fastq")
# 	conda:
# 		'../env/genotyping.yml'
# 	resources:
# 		mem_mb=20000,
# 		runtime_hrs=5,
# 		runtime_min=1,
# 		time=5*60,
# 		partition="med2"
# 	log:
# 		"{results}/{sample}/raw/{sample}-{fraction, [0-9.]+}.downsample.log"
# 	shell:
# 		"bash ../scripts/downsample-fasta.sh {input.coverage} {wildcards.fraction} {input.reads1} {input.reads2} {output.sampled1} {output.sampled2} &> {log}"


# #### data for mapping free approaches ####

# def combine_reads_input(wildcards):
# 	if wildcards.fraction == "full":
# 		return config['data'][wildcards.sample]['reads']
# 	else:
# 		return expand("{results}/{sample}/raw/{sample}-{fraction}_{r}.fastq", results=wildcards.results, sample=wildcards.sample, fraction=wildcards.fraction, r=[1,2])

# # generate combined fastq file
# rule combine_reads:
# 	input:
# 		combine_reads_input
# 	output:
# 		"{results}/{sample}/raw/{sample}-{fraction, (full|[0-9.]+)}.fastq"
# 	shell:
# 		"cat {input} > {output}"


#### data for mapping based approaches ####

# index fasta
rule bwa_index:
	input:
		input_reference
	output:
		input_reference + ".ann"
	log:
		"{results}/reference-indexing.log".format(results=outname)
	conda:
		'../env/genotyping.yml'
	resources:
		mem_mb=5000,
		runtime=3*60
	shell:
		"(/usr/bin/time -v bwa index {input}) &> {log}"

# create fasta.fai file
rule samtools_faidx:
	input:
		input_reference
	output:
		input_reference + '.fai'
	conda:
		'../env/genotyping.yml'
	threads: 1
	resources:
		mem_mb=1024,
		runtime = lambda wildcards, attempt: 60  *4* attempt,
		partition = "med2"
	shell:
		"samtools faidx {input}"

# align illumina reads
def bwa_mem_input(wildcards):
	if wildcards.fraction == 'full':
		return [config['data'][wildcards.sample]['reads'][0], config['data'][wildcards.sample]['reads'][1]]
	else:
		return expand("{results}/{sample}/raw/{sample}-{fraction}_{r}.fastq", results=wildcards.results, sample=wildcards.sample, fraction=wildcards.fraction, r=[1,2])

rule bwa_mem:
	input:
		reads=bwa_mem_input,
		fasta=input_reference,
		index=input_reference + '.ann',
		fai=input_reference + '.fai'
	output:
		'{results}/{sample}/aligned/{sample}-{fraction, (full|[0-9.]+)}.bam'
	log:
		'{results}/{sample}/aligned/{sample}-{fraction, (full|[0-9.]+)}.log'
	threads: 24
	resources:
		mem_mb=60000,
		runtime_hrs=25,
		runtime_min=1,
		time=25*60,
		meduim=1,
		partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}")
	conda:
		'../env/genotyping.yml'
	shell:
		'(/usr/bin/time -v bwa mem -t {threads} -M {input.fasta} -R "@RG\\tID:{wildcards.sample}\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:{wildcards.sample}" {input.reads} | samtools view -bS | samtools sort -o {output} - ) &> {log}'

# index BAM file
rule samtools_index:
	input:
		"{filename}.bam"
	output:
		"{filename}.bam.bai"
	log:
		"{filename}-index.log"
	conda:
		'../env/genotyping.yml'
	threads: 1
	resources:
		mem_mb=1024,
		runtime = lambda wildcards, attempt: 60  *4* attempt,
		partition = "med2"
	shell:
		"(/usr/bin/time -v samtools index {input}) &> {log}"

# split BAM by chromosome
rule split_bam_by_chromosome:
	input:
		bam='{prefix}/{sample}-{fraction, [0-9.]+}.bam',
		bai='{prefix}/{sample}-{fraction, [0-9.]+}.bam.bai'
	output:
		'{prefix}/{sample}-{fraction}.chr{chrom, X|Y|[0-9]+}.bam'
	conda:
		'../env/genotyping.yml'
	threads: 1
	resources:
		mem_mb=1024,
		runtime = lambda wildcards, attempt: 60  *4* attempt,
		partition = "med2"
	shell:
		'samtools view -h {input.bam} {prefix}{wildcards.chrom} | samtools view -Sb -> {output}'
