configfile: 'config.json'
include: "prepare-reads.smk"

tempFolder= config['parameters']['temp']
### Parameters ###

# input data for genotyping
samples = []
for sample,data in config['data'].items():
    if isinstance(data, dict):
        samples.append(sample)
        samples.extend(data['sample'])



input_reference = config['data']['reference']
reference_version = config['data']['reference_version']
repeats_bed = config['data']['repeats']
complex_bed = config['data']['complex']

# programs
pangenie = config['programs']['pangenie']
#platypus = config['programs']['platypus']
# kmc = config['programs']['kmc']
# bayestyper = config['programs']['bayestyper']
# bayestyper_tools = config['programs']['bayestyper_tools']
paragraph = config['programs']['paragraph']
graphtyper = config['programs']['graphtyper']
GG = config['programs']['TheGreatGenotyper']
beagle= config['programs']['beagle']
beagleMap= config['programs']['beagleMap']
gatk = 'gatk'
picard = 'picard'

# parameters
chromosomes = config['parameters']['chromosomes']
# outname = config['parameters']['outname']
# outname_reads = config['parameters']['outname_reads']
#downsampling = [5, 10, 20, 30]
downsampling = [30]
GG_INDEX_Prefix=  config['parameters']['GG_INDEX']
#bayestyper_reference_canon = config['parameters']['bayestyper_reference_canon']
#bayestyper_reference_decoy = config['parameters']['bayestyper_reference_decoy']
#other_methods=['platypus', 'bayestyper', 'gatk', 'paragraph', 'graphtyper']
other_methods = ['gatk', 'paragraph', 'graphtyper', "TheGreatGenotyper-kmersOnly","TheGreatGenotyper-SecondPassHMM", "TheGreatGenotyper-HMM" ]

variants = ['snp', 'small-deletion', 'small-insertion', 'midsize-deletion', 'midsize-insertion', 'large-deletion',
            'large-insertion', 'indel', 'sv', 'small', 'midsize', 'large']

metric_to_script = {
    'precision-recall-all': '../scripts/plot-precision-recall.py',
    'precision-recall-typable': '../scripts/plot-precision-recall.py',
    'concordance': '../scripts/plot-concordances.py',
    'fscore': '../scripts/plot-fscores.py',
    'newFscore': '../scripts/plot_accuracy_new.py'
}

import random

def getHighPartition(sample):
    partitions=["high2","bmh"]
    p=hash(sample)%len(partitions)
    return partitions[p]

def getMeduimPartition(sample):
    partitions=["bmm"]
    p=hash(sample)%len(partitions)
    return partitions[p]

def getLowPartition(sample,attempt):
    if attempt>1:
       return getMeduimPartition(sample)
    partitions=["med2"]
    p=hash(sample)%len(partitions)
    return partitions[p]




localrules: extract_vcf,tabix,split_vcf_by_chromosome,paragraph_manifest,paragraph_preprocessing,paragraph_vcf,paragraph_postprocessing,normalize_callsets,merge_vcfs,bayestyper_make_samples_file,graphtyper_preprocess,prepare_regions,evaluation_region_graph,evaluation_region_external,rtg_format,plot_results,bgzip

##############################################################
##################    prepare input VCF    ##################
##############################################################


# uncompress vcf
rule extract_vcf:
    input: "{results}/{sample}/{sample_test}/truth/no{sample}_{mode}.vcf.gz"
    output:
        "{results}/{sample}/{sample_test}/{mode, biallelic|reference-panel}/{sample}-all.vcf"
    conda:
        "../env/genotyping.yml"
    shell:
        "gunzip -c {input} > {output}  "

# rule extract_gold:
#     input: lambda wildcards: config['data'][wildcards.sample]['truth']
#     output:
#         "{results}/{sample}/{sample_test}/truth/only{sample}-test.vcf"
#     conda:
#         "../env/genotyping.yml"
#     shell:
#         "cat {input} > {output}"


ruleorder: split_vcf_by_chromosome >  remove_untypable > paragraph_preprocessing > paragraph_genotyping >bgzip
# rule bgzip
rule bgzip:
    input:
        "{filename}.vcf"
    output:
        "{filename}.vcf.gz"
    conda:
        "../env/genotyping.yml"
    shell:
        "bgzip -c {input} > {output}"

# rule tabix
rule tabix:
    input:
        "{filename}.vcf.gz"
    output:
        "{filename}.vcf.gz.tbi"
    conda:
        "../env/genotyping.yml"
    shell:
        "tabix -p vcf {input}"

# def getTestVCF(wildcards):
#     vcf="%s/%s/%s/truth/no%s_biallelic.vcf.gz"%(wildcards.results,wildcards.sample,wildcards.sample_test,wildcards.sample)"
#     if wildcards.mode != "reference-panel":
#        vcf="%s/%s/%s/truth/no%s_biallelic.vcf.gz"%(wildcards.results,wildcards.sample,wildcards.sample_test,wildcards.sample)"
#     return [vcf,vcf+".tbi"]
    
# split VCF by chromosome
rule split_vcf_by_chromosome:
    input: 
        vcf="{results}/{sample}/{sample_test}/truth/no{sample}_{mode}.vcf.gz", 
        tbi="{results}/{sample}/{sample_test}/truth/no{sample}_{mode}.vcf.gz.tbi" 
    output:
        vcf="{results}/{sample}/{sample_test}/{mode}/{sample}-chr{chrom}.vcf",
        gz="{results}/{sample}/{sample_test}/{mode}/{sample}-chr{chrom}.vcf.gz"
    wildcard_constraints:
        chrom="X|Y|[0-9]+",
        mode="reference-panel|biallelic"
    params:
        prefix='' if ('37' in reference_version) or ('19' in reference_version) else 'chr'
    conda:
        "../env/genotyping.yml"
    shell:
        """
        bcftools view {input.vcf} -r {params.prefix}{wildcards.chrom} > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
    #        tabix -p vcf {output.gz}
        """


#########################################################
#################    Prepare running scripts     ########
#########################################################

rule getSampleNames:
    input:
        vcf=lambda wildcards: config["data"][wildcards.sample]["sample"][wildcards.sample_test],
    output:
        "{results}/{sample}/{sample_test}/truth/sample_names.txt"
    conda:
        "env.yaml"
    wildcard_constraints:
        sample="|".join(samples),
        sample_test="|".join(samples)
    resources:
        mem_mb=10000,
        runtime_hrs=1,
        runtime_min=1,
        runtime = lambda wildcards, attempt: 60  * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt)
    shell:
        "zgrep -P '^#[^#]' {input} | cut -f10- |tr -s $'\t' ',' >{output}"



# create a multisample VCF containing sample and gold
rule mergeGoldandCallset:
    input:
        vcf_1=lambda wildcards: config["data"][wildcards.sample]["truth"]+".gz",
        tbi_1=lambda wildcards: config["data"][wildcards.sample]["truth"]+".gz.tbi",
        vcf_2=lambda wildcards: config["data"][wildcards.sample]["sample"][wildcards.sample_test],
        tbi_2=lambda wildcards: config["data"][wildcards.sample]["sample"][wildcards.sample_test]+".tbi"
    output:
        "{results}/{sample}/{sample_test}/truth/merged-variants.vcf"
    conda:
        "env.yaml"
    wildcard_constraints:
        sample="|".join(samples),
        sample_test="|".join(samples)
    resources:
        mem_mb=10000,
        runtime_hrs=1,
        runtime_min=1,
        runtime = lambda wildcards, attempt: 60  * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt)
    shell:
        "bcftools merge -m none --missing-to-ref {input.vcf_1} {input.vcf_2}|grep -v '<DUP>'| python3 ../scripts/assign-variant-ids.py   >  {output}  "



def subsets_to_samples(wildcards):
    if wildcards.subset == "all":
       return wildcards.sample+",$samples"
    elif wildcards.subset == "no"+wildcards.sample:
       return "$samples"
    elif wildcards.subset == "only"+wildcards.sample:
       return wildcards.sample

rule merge_haplotypes:
    input:
        vcf = "{results}/{sample}/{sample_test}/truth/merged-variants.vcf.gz",
        tbi = "{results}/{sample}/{sample_test}/truth/merged-variants.vcf.gz.tbi",
        names = "{results}/{sample}/{sample_test}/truth/sample_names.txt",
        reference = config['data']['reference']
    output:
        biallelic= "{results}/{sample}/{sample_test}/truth/{subset}_biallelic.vcf",
        vcf = "{results}/{sample}/{sample_test}/truth/{subset}_reference-panel.vcf"
                # results/HG00731/NA12878/truth/noHG00731_reference-panel.vcf
    # wildcard_constraints:
    #     subset = "{sample_test}|{sample}"
    wildcard_constraints:
        sample="|".join(samples),
        sample_test="|".join(samples),
        #subset = "{sample_test}|{sample}|all"
    params:
        chrom = ','.join(["chr"+c for c in chromosomes]),
        samples = lambda wildcards : subsets_to_samples(wildcards)
    log:
        "{results}/{sample}/{sample_test}/truth/{subset}.log"
    conda:
        "env.yaml"
    resources:
        mem_mb=10000,
        runtime_hrs=1,
        runtime_min=1,
        runtime = lambda wildcards, attempt: 60  * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt)
    shell:
        """
        samples=$(cat {input.names} |tr -d $'\n')
        echo bcftools view -s {params.samples} {input.vcf}
        bcftools view -s {params.samples} {input.vcf} | bcftools view --min-ac 1 > {output.biallelic}
        python3 ../scripts/merge_vcfs.py merge -vcf {output.biallelic} -r {input.reference} -ploidy 2 -chromosomes {params.chrom} 2> {log} 1> {output.vcf}

        """


rule untypable_ids:
    input:
        "{results}/{sample}/{sample_test}/truth/all_biallelic.vcf"
    output:
        ids1="{results}/{sample}/{sample_test}/truth/untypable-ids/{sample}-untypable.tsv",
        #ids2="{results}/{sample}/{sample_test}/truth/untypable-ids/{sample_test}-untypable.tsv",
        summary= "{results}/{sample}/{sample_test}/truth/untypable-ids_summary.tsv"
    params:
        out= "{results}/{sample}/{sample_test}/truth/untypable-ids/"
    resources:
        mem_mb=10000,
        runtime_hrs=1,
        runtime_min=1,
        runtime = lambda wildcards, attempt: 60  * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt)
    conda:
        "env.yaml"
    shell:
        "cat {input} | python3 ../scripts/untypable-ids.py {params.out} > {output.summary}"



#########################################################
#################    run paragraph     ##################
#########################################################


# compute depth
rule paragraph_depth:
    input:
        bam='{prefix}/{sample}-{fraction}.bam',
        bai='{prefix}/{sample}-{fraction}.bam.bai',
        fasta=input_reference
    output:
        '{prefix}/depth_{sample}-{fraction}.json'
    log:
        '{prefix}/depth_{sample}-{fraction}.log'
    resources:
        mem_mb=6000,
        runtime_hrs=0,
        runtime_min=10,
        runtime = lambda wildcards, attempt: 10  * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),    
    threads:
        24
    conda:
        "env.yaml"
    shell:
        "/usr/bin/time -v idxdepth -b {input.bam} -r {input.fasta} -o {output} --threads {threads} &> {log}"

# write manifest file
rule paragraph_manifest:
    output:
        "{results}/{sample}/{sample_test}/paragraph/manifest/manifest_{sample}-{fraction}.txt"
    params:
        bam= lambda wildcards: config['data'][wildcards.sample]["mapping_prefix"] + wildcards.sample + '-' + wildcards.fraction+ '.bam',
        json= lambda wildcards: config['data'][wildcards.sample]["mapping_prefix"] + "depth_" + wildcards.sample + "-" + wildcards.fraction + ".json"
    run:
        with open(output[0],"w") as paragraph_manifest_file:
            paragraph_manifest_file.write("id\tpath\tidxdepth\n")
            paragraph_manifest_file.write("{sample}\t{bam}\t{json}\n".format(sample=wildcards.sample,bam=params.bam,json=params.json))  

# produce input vcf file
rule paragraph_preprocessing:
    input:
        vcf="{results}/{sample}/{sample_test}/reference-panel/{sample}-all.vcf",
        fasta=input_reference
    output:
        "{results}/{sample}/{sample_test}/paragraph/variants/variants_{sample}.vcf.gz"
    conda:
        "../env/genotyping.yml"
    shell:
        "bash ../scripts/paragraph-preprocess.sh {input.fasta} {input.vcf} {output} "

# run genotyping algorithm
rule paragraph_genotyping:
    input:
        vcf = "{results}/{sample}/{sample_test}/paragraph/variants/variants_{sample}.vcf.gz",
        tbi = "{results}/{sample}/{sample_test}/paragraph/variants/variants_{sample}.vcf.gz.tbi",
        manifest = "{results}/{sample}/{sample_test}/paragraph/manifest/manifest_{sample}-{fraction}.txt",
        fasta = input_reference,
        depth = lambda wildcards: config['data'][wildcards.sample]["mapping_prefix"] + "depth_" + wildcards.sample + "-" + wildcards.fraction + ".json" 
    output:
        temp("{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}/genotypes.vcf.gz"),
        temp("{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}/genotypes.json.gz"),
        temp("{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}/variants.json.gz"),
        temp("{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}/variants.vcf.gz")
    log:
        genotyping="{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}/genotypes.log"
    params:
        outname="{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}",
        max_depth=lambda wildcards: int(wildcards.fraction) * 20
    threads: 24
    wildcard_constraints:
        sample="|".join(samples),
        sample_test="|".join(samples)
    resources:
        mem_mb=30000,
        runtime_hrs=10,
        runtime=10 * 60,
        meduim=1,
        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}")
    conda:
        "env.yaml"
    shell:
        "/usr/bin/time -v multigrmpy.py -i {input.vcf} -M {params.max_depth} -m {input.manifest} -r {input.fasta} -o {params.outname} --threads {threads} --scratch-dir {params.outname}/tmp &> {log.genotyping} "

rule paragraph_vcf:
    input:
        "{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}/genotypes.vcf.gz"
    output:
        temp("{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}/genotypes.vcf")
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        sample="|".join(samples),
        sample_test="|".join(samples)
    shell:
        "gunzip -c {input} > {output}"

# postprocess paragraph output
rule paragraph_postprocessing:
    input:
        "{results}/{sample}/{sample_test}/paragraph/{sample}-{fraction}/genotypes.vcf"
    output:
        "{results}/{sample}/{sample_test}/paragraph/{sample}_{fraction}_paragraph_genotyping.vcf"
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        sample="|".join(samples),
        sample_test="|".join(samples)
    shell:
        "python3 ../scripts/paragraph-helper.py postprocess {input} > {output}"


#########################################################
##################      run GATK       ##################
#########################################################

# sequence dictionary
rule sequence_dictionary:
    input:
        "{filename}.fa"
    output:
        "{filename}.dict"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=10000,
        runtime_hrs=1,
        runtime=60
    shell:
        "{picard} CreateSequenceDictionary R={input} O={output}"

# mark duplicates
rule mark_duplicates:
    input:
        bam='{prefix}/{sample}-{fraction}.chr{chrom}.bam',
        bai='{prefix}/{sample}-{fraction}.chr{chrom}.bam.bai'
    output:
        bam= '{prefix}/marked/{sample}-{fraction}.chr{chrom}.bam',
        metrics= "{prefix}/marked/{sample}-{fraction}.chr{chrom}.marked.txt"
    log:
        "{prefix}/marked/{sample}-{fraction}.chr{chrom}.marked.log"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),    
    shell:
        "bash ../scripts/mark-duplicates.sh {input.bam} {output.bam} {output.metrics} {log}"

# prepare input vcf (fix ##config in header)
rule gatk_preprocessing:
    input:
        vcf="{results}/{sample}/{sample_test}/biallelic/{sample}-chr{chrom, X|Y|[0-9]+}.vcf",
        fai=config['data']['reference'] + '.fai'
    output:
        "{results}/{sample}/{sample_test}/gatk/variants/variants_{sample}.chr{chrom, X|Y|[0-9]+}.vcf.gz"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        "python3 ../scripts/gatk-helper.py preprocess {input.fai} {input.vcf} -l 50 | bgzip >   {output}"

#        shell("tabix -p vcf {output}")

# retype variants
rule haplotype_caller_retype:
    input:
        bam= lambda wildcards: getBAM(wildcards),
        bai= lambda wildcards: getBAM(wildcards) +'.bai',
        reference=config['data']['reference'],
        reference_dict=config['data']['reference'][:-3] + '.dict',
        vcf="{results}/{sample}/{sample_test}/gatk/variants/variants_{sample}.chr{chrom, X|Y|[0-9]+}.vcf.gz",
        tbi="{results}/{sample}/{sample_test}/gatk/variants/variants_{sample}.chr{chrom, X|Y|[0-9]+}.vcf.gz.tbi"
    output:
        "{results}/{sample}/{sample_test}/gatk/{sample}_{fraction}_gatk_genotyping-raw.chr{chrom, X|Y|[0-9]+}.vcf"
    log:
        "{results}/{sample}/{sample_test}/gatk/{sample}_{fraction}_gatk_genotyping.chr{chrom, X|Y|[0-9]+}.log"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=10000,
        runtime_hrs=4,
        runtime_min=59,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt)
    shell:
        "/usr/bin/time -v {gatk} HaplotypeCaller --reference {input.reference}  --input {input.bam} --output {output} --intervals chr{wildcards.chrom} --minimum-mapping-quality 20 --genotyping-mode GENOTYPE_GIVEN_ALLELES --alleles {input.vcf} &> {log}"

# discover variants
rule haplotype_caller_discover:
    input:
        bam= lambda wildcards: getBAM(wildcards),
        bai= lambda wildcards: getBAM(wildcards) +'.bai',
        reference=config['data']['reference'],
        reference_dict=config['data']['reference'][:-3] + '.dict'
    output:
        "{results}/{sample}/{sample_test}/gatk/{sample}_{fraction}_gatk_discovery-raw.chr{chrom, X|Y|[0-9]+}.vcf"
    log:
        "{results}/{sample}/{sample_test}/gatk/{sample}_{fraction}_gatk_discovery.chr{chrom, X|Y|[0-9]+}.log"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=10000,
        runtime_hrs=4,
        runtime_min=59,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        "/usr/bin/time -v {gatk} HaplotypeCaller --reference {input.reference} --input {input.bam} --output {output} --intervals chr{wildcards.chrom} --minimum-mapping-quality 20 --genotyping-mode DISCOVERY &> {log}"


#########################################################
##################    run platypus     ##################
#########################################################


# run Platypus
rule platypus_retype:
    input:
        vcf="{results}/{sample}/{sample_test}/biallelic/{sample}-chr{chrom, X|Y|[0-9]+}.vcf.gz",
        tbi="{results}/{sample}/{sample_test}/biallelic/{sample}-chr{chrom, X|Y|[0-9]+}.vcf.gz.tbi",
        fasta=input_reference,
        fai=input_reference + ".fai",
        bam= lambda wildcards: getBAM(wildcards),
        bai= lambda wildcards: getBAM(wildcards) +'.bai'
    output:
        "{results}/{sample}/{sample_test}/platypus/{sample}_{fraction}_platypus_genotyping-raw.chr{chrom, X|Y|[0-9]+}.vcf"
    log:
        "{results}/{sample}/{sample_test}/platypus/{sample}_{fraction}_platypus_chr{chrom, X|Y|[0-9]+}.log"
    conda:
        "env-Platypus.yaml"
    resources:
        mem_mb=3000,
        runtime_hrs=4,
        runtime_min=59,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        "bash ../scripts/run-platypus.sh {platypus} {input.bam} {input.fasta} {output} {input.vcf} {log}"


# run Platypus in discovery mode
rule platypus_discovery:
    input:
        fasta=input_reference,
        fai=input_reference + ".fai",
        bam= lambda wildcards: getBAM(wildcards),
        bai= lambda wildcards: getBAM(wildcards) +'.bai'
    output:
        "{results}/{sample}/{sample_test}/platypus/{sample}_{fraction}_platypus_discovery-raw.chr{chrom, X|Y|[0-9]+}.vcf"
    log:
        "{results}/{sample}/{sample_test}/platypus/{sample}_{fraction}_platypus_chr{chrom, X|Y|[0-9]+}.log"
    resources:
        mem_mb=3000,
        runtime_hrs=1,
        runtime_min=59,
        runtime = lambda wildcards, attempt: 60  * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    conda:
        "env-Platypus.yaml"
    shell:
        "bash ../scripts/run-platypus-discover.sh {platypus} {input.bam} {input.fasta} {output} {log}"


# normalize platypus/gatk/graphtyper VCFs
rule normalize_callsets:
    input:
        vcf="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}-raw.chr{chrom}.vcf",
        ref=input_reference
    output:
        "{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}.chr{chrom}.vcf"
    wildcard_constraints:
        method="platypus|gatk",
        run_mode="genotyping|discovery",
        chrom="X|Y|[0-9]+"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        """
    bgzip -c {input.vcf} > {input.vcf}.gz
    tabix -p vcf {input.vcf}.gz
        bcftools norm -f {input.ref} {input.vcf}.gz -m -any > {output}

        """

# merge chromosome-wise VCFs
rule merge_vcfs:
    input:
        expand("{{results}}/{{sample}}/{{sample_test}}/{{method}}/{{sample}}_{{fraction}}_{{method}}_{{run_mode}}.chr{chrom}.vcf",chrom=chromosomes)
    output:
        "{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}.vcf"
    wildcard_constraints:
        method="platypus|gatk|graphtyper|graphtyper",
        run_mode="genotyping|discovery|graphtyper"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    run:
        shell("grep '^#' {input[0]} > {output}")
        shell("grep -v '^#' --no-filename {input} >> {output}")


##############################################################################
##################    run BayesTyper genotyping pipeline    ##################
##############################################################################


# # run kmc to count kmers
# rule run_kmc:
#     input:
#         outname_reads + "/{sample}/raw/{sample}-{fraction, [0-9.]+}.fastq",
#     output:
#         suf="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_suf",
#         pre="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_pre",
#     #        tmp=temp("{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/tmp-{sample}/")
#     log:
#         "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}_kmc.log"
#     params:
#         out_prefix="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}",
#         tmp="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/tmp-{sample}/"
#     threads: 40
#     conda:
#         "env.yaml"
#     resources:
#         mem_mb=15000,
#         runtime_hrs=1,
#         runtime_min=59,
#         runtime=60,
#         meduim=1,
#         partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}")
#     shell:
#         "/usr/bin/time -v {kmc} -k55 -t{threads} -ci1 {input} {params.out_prefix} {params.tmp} > {log} 2>&1"

# create bloomfilter
rule create_bloomfilter:
    input:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_pre",
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_suf"
    output:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.bloomData",
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.bloomMeta"
    log:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}_bloom.log"
    params:
        out_prefix="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}"
    conda:
        "env.yaml"
    resources:
        mem_mb=20000,
        runtime_hrs=4,
        runtime_min=59,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        "/usr/bin/time -v {bayestyper_tools} makeBloom -k {params.out_prefix} > {log} 2>&1"

# create samples file
rule bayestyper_make_samples_file:
    output:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/{sample}.tsv"
    params:
        prefix="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}",
        sex=lambda wildcards: config['data'][wildcards.sample]['sex']
    run:
        with open(output[0],"w") as bayestyper_samples_file:
            bayestyper_samples_file.write("{sample}\t{sex}\t{prefix}\n".format(sample=wildcards.sample,sex=params.sex,prefix=params.prefix))

# bayestyper cluster
checkpoint bayestyper_cluster:
    input:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.bloomData",
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.bloomMeta",
        variants="{results}/{sample}/{sample_test}/reference-panel/{sample}-all.vcf",
        samples="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/{sample}.tsv"
    output:
        dir=directory("{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/clusters/")
    params:
        out_prefix="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/clusters/bayestyper"
    log:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/clusters/clusters-log.log"
    resources:
        mem_mb=60000,
        runtime_hrs=4,
        runtime_min=59,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    threads: 24
    conda:
        "env.yaml"
    shell:
        """
        mkdir -p {output.dir}
        /usr/bin/time -v {bayestyper} cluster -v {input.variants} -s {input.samples} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} \
                    -p {threads} -o {params.out_prefix} > {log} 2>&1
        """

# bayestyper genotype
rule run_bayestyper_genotype:
    input:
        kmc_pre="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_pre",
        kmc_suf="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_suf",
        samples="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/{sample}.tsv",
        unit="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/clusters/bayestyper_unit_{unit_id}/variant_clusters.bin"
    output:
        genotypes="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
        kmer_coverage_file="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper_genomic_parameters.txt"
    log:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper.log"
    params:
        cluster_data_dir="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/clusters/bayestyper_cluster_data",
        out_prefix="{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper"
    threads: 24
    conda:
        "env.yaml"
    resources:
        mem_mb=50000,
        runtime_hrs=1,
        runtime_min=59,
        runtime = lambda wildcards, attempt: 60  * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        """        
        /usr/bin/time -v {bayestyper} genotype -v {input.unit} -s {input.samples} -c {params.cluster_data_dir} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} \
        -p {threads} -z -o {params.out_prefix} > {log} 2>&1
        # fix the vcf ...
        gunzip -c {output.genotypes} | bgzip > {output.genotypes}-tmp
        mv {output.genotypes}-tmp  {output.genotypes}
        tabix -p vcf {output.genotypes}
        """


# combine vcfs
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output[0]
    result = expand("{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
        results=wildcards.results,
        sample=wildcards.sample,
        fraction=wildcards.fraction,
        unit_id=glob_wildcards(os.path.join(checkpoint_output,"bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
    return sorted(result)


rule bcftools_concat_units:
    input:
        aggregate_input
    output:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}_bayestyper_genotyping.vcf"
    log:
        "{results}/{sample}/{sample_test}/bayestyper/{sample}_{fraction}_bayestyper_genotyping.log"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        "/usr/bin/time -v bcftools concat -a -o {output} {input} &> {log}"


########################################################
##################    run PanGenie    ##################
########################################################

# run pangenie
rule pangenie:
    input:
        reads= lambda wildcards: config['data'][wildcards.sample]["reads_prefix"] + wildcards.sample + "-" + wildcards.fraction + ".fastq",
        fasta=input_reference,
        vcf="{results}/{sample}/{sample_test}/reference-panel/{sample}-all.vcf",
    output:
        "{results}/{sample}/{sample_test}/pangenie/{sample}_{fraction}_pangenie_genotyping.vcf"
    log:
        "{results}/{sample}/{sample_test}/pangenie/{sample}_{fraction}_pangenie_log.log"
    threads: 24
    params:
        out_prefix="{results}/{sample}/{sample_test}/pangenie/{sample}_{fraction}_pangenie"
    resources:
        mem_mb=100000,
        runtime_hrs=4,
        runtime_min=59,
        runtime=24 * 60,
        meduim=1,
        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}")    
    shell:
        """
        (/usr/bin/time -v {pangenie}  -i {input.reads}    -v <(sed -e 's/\([0-9]\)\/\([0-9]\)/\1|\2/g' {input.vcf}) -r {input.fasta} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g) &> {log}   
    #
     """



########################################################
############    run The Great Genotyper    #############
########################################################

#GG_INDEX_Prefix="/group/ctbrowngrp2/mshokrof/metagraph_2/smooth_1_noLog_noClean/"
#GG_INDEX_Prefix="/group/ctbrowngrp2/mshokrof/metagraph_3/smooth_1_noLog_noClean/"
GG_log= ""
# # run The Great Genotyper
# rule TheGreatGenotyper:
#     input:
#         graph      = lambda wildcards: config['data'][wildcards.sample]['metagraph_index'] + "graph.dbg",
#         annotation = lambda wildcards: config['data'][wildcards.sample]['metagraph_index'] + "annotation.relaxed.row_diff_int_brwt.annodbg",
#         desc       = lambda wildcards: config['data'][wildcards.sample]['metagraph_index'] + "graph.desc.tsv",    
#         ref        = input_reference,
#         vcf        = "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.vcf",
#     output:
#         expand("{results}/{sample}/{{sample_test}}/TheGreatGenotyper/{sample}_{fraction}_TheGreatGenotyper_genotyping.vcf",results= "{results}",sample= "{sample}",fraction=downsampling) 
#     log:
#         "{results}/{sample}/{sample_test}/TheGreatGenotyper/{sample}.log"
#     threads: 32
#     retries: 1
#     params:
#         out_prefix="{results}/{sample}/{sample_test}/TheGreatGenotyper/"
#     resources:
#         mem_mb=100000,
#         cores=24,
#         nodes = 1,
#         meduim=1,
#         time = 60 * 48,
#         partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}"),
#         #tmp= lambda wildcards: "%sgenotyper_%s_%s/"%(tempFolder,f"{wildcards.sample}",f"{wildcards.graph}")
#     shell:
#         """
#     /home/mshokrof/TheGreatGenotyper/build2/pangenie/src/TheGreatGenotyper {GG_log} -g -a {input.annotation} -i {input.graph} -f {input.desc} -j {threads} -t {threads} -r {input.ref} -v <(sed -e 's/\([0-9]\)\/\([0-9]\)/\1|\2/g' {input.vcf} ) -o {params.out_prefix}  &> {log}
#     cd {params.out_prefix}
#      ls _*genotyping.vcf  | sed -e 's/^_//' | parallel --gnu 'mv _{{}} {{}}'
#     ####
#     """

rule unphase_vcf:
    input:
        "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.vcf"
    output:
        "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.unphased.vcf"
    log:
        "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.unphased.log"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=5000,
        runtime=2 * 60,
        meduim=1,
        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}")
    shell:
        """
    bgzip -c {input} > tmp.$$.vcf.gz
    tabix -p vcf tmp.$$.vcf.gz
    python ../scripts/unphase.py tmp.$$.vcf.gz {output}
    rm  tmp.$$.vcf.gz    
    """




    
    
        
rule TheGreatGenotyperkmersOnly:
    input:
        graphFolders      = GG_INDEX_Prefix+"chunk.{chunk}",    
        ref        = input_reference,
        vcf        = "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.unphased.vcf",
    output:
        merged_vcf="{results}/{sample}/{sample_test}/TheGreatGenotyper-kmersOnly/population.chunk-{chunk}.vcf.gz"
    log:
        "{results}/{sample}/{sample_test}/TheGreatGenotyper-kmersOnly/population.{chunk}.log"
    threads: 32
    retries: 0
    resources:
        mem_mb= 200 *1024,
        cores=24,
        nodes = 1,
        meduim=1,
        runtime = 60 * 24,
        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}"),
        tmpdir= lambda wildcards: "%sTheGreatGenotyper-kmersOnly_%s_%s/"%(tempFolder,f"{wildcards.sample}",f"{wildcards.chunk}")
    shell:
        r"""
            mkdir -p {resources.tmpdir}
            {GG} -f  -a   -g  -i {input.graphFolders}  -j {threads} -t {threads} -r {input.ref}  -y  {resources.tmpdir}emissions -v {input.vcf} -o -  2> {log} | bgzip > {resources.tmpdir}/tmp.vcf.gz 
            mv {resources.tmpdir}/tmp.vcf.gz {output.merged_vcf}
            tabix -p vcf {output.merged_vcf}
            rm -rf {resources.tmpdir}
             ###
        """
rule rephase_vcf:
    input:
        vcf= "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.unphased.vcf",
        panel= "{results}/{sample}/{sample_test}/TheGreatGenotyper-kmersOnly/population.vcf.gz"
    output:
        "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.rephased.vcf"
    log:
        "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.rephased.log"
    params:
        prefix="{results}/{sample}/{sample_test}/reference-panel/{sample}-all.rephased"
    threads: 64
    resources:
        mem_mb=80*1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 10,
        tmpdir= lambda wildcards: "%srephase_%s/"%(tempFolder,f"{wildcards.sample}"),
        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}"),
    conda:
        "../env/genotyping.yml"
    shell:
        """
            mkdir -p {resources.tmpdir}
            bcftools view --force-samples  -s ^NA12878 {input.panel}  2> {log}|bgzip > {resources.tmpdir}panel.vcf.gz
            tabix -p vcf {resources.tmpdir}panel.vcf.gz
            java -Xmx70G -jar {beagle} gt={resources.tmpdir}panel.vcf.gz out={resources.tmpdir}tmp2.$$ nthreads={threads}  map={beagleMap}  >> {log}
            tabix -p vcf {resources.tmpdir}tmp2.$$.vcf.gz
            bcftools view -Oz -o {resources.tmpdir}tmp.input.$$.vcf.gz {input.vcf}
            java -Xmx70G -jar {beagle} gt={resources.tmpdir}tmp.input.$$.vcf.gz  ref={resources.tmpdir}tmp2.$$.vcf.gz out={params.prefix} nthreads={threads}  map={beagleMap}  >> {log}
            gzip -d {params.prefix}.vcf.gz
            rm -rf {resources.tmpdir}    
        """






rule TheGreatGenotyperHMM:
    input:
        graphFolders      =GG_INDEX_Prefix+"chunk.{chunk}",    
        ref        = input_reference,
        vcf        = "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.vcf",
    output:
        merged_vcf="{results}/{sample}/{sample_test}/TheGreatGenotyper-HMM/population.chunk-{chunk}.vcf.gz"
    log:
        "{results}/{sample}/{sample_test}/TheGreatGenotyper-HMM/population.{chunk}.log"
    threads: 32
    retries: 0
    resources:
        mem_mb= 400*1024,
        cores=24,
        nodes = 1,
        meduim=1,
        runtime = 60 * 48,
        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}"),
        tmpdir= lambda wildcards: "%sTheGreatGenotyper-HMM_%s_%s/"%(tempFolder,f"{wildcards.sample}",f"{wildcards.chunk}")
    shell:
        r"""
            mkdir -p {resources.tmpdir}
            {GG}  -f  -g  -i {input.graphFolders}  -j {threads} -t {threads} -r {input.ref}  -y   {resources.tmpdir}emissions.TheGreatGenotyperHMM -v {input.vcf} -o -  2> {log} | bgzip > {resources.tmpdir}/tmp.vcf.gz
            mv {resources.tmpdir}/tmp.vcf.gz {output.merged_vcf}
            tabix -p vcf {output.merged_vcf}
            rm -rf {resources.tmpdir}
            ####
        """

rule TheGreatGenotyperSecondPassHMM:
    input:
        graphFolders = GG_INDEX_Prefix+"chunk.{chunk}",    
        ref        = input_reference,
        vcf        = "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.rephased.vcf",
    output:
        merged_vcf="{results}/{sample}/{sample_test}/TheGreatGenotyper-SecondPassHMM/population.chunk-{chunk}.vcf.gz"
    log:
        "{results}/{sample}/{sample_test}/TheGreatGenotyper-SecondPassHMM/population.{chunk}.log"
    threads: 32
    retries: 0
    resources:
        mem_mb=400*1024,
        cores=24,
        nodes = 1,
        meduim=1,
        runtime = 60 * 48,
        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}"),
        tmpdir= lambda wildcards: "%sTheGreatGenotyper-SecondPass_%s_%s/"%(tempFolder,f"{wildcards.sample}",f"{wildcards.chunk}")
    shell:
        r"""
            mkdir -p {resources.tmpdir}
            {GG}   -f  -g  -i {input.graphFolders}  -j {threads} -t {threads} -r {input.ref}  -y  {resources.tmpdir}emissions.TheGreatGenotyperSecondPassHMM -v {input.vcf} -o -  2> {log} | bgzip > {resources.tmpdir}/tmp.vcf.gz
            mv {resources.tmpdir}/tmp.vcf.gz {output.merged_vcf}
            tabix -p vcf {output.merged_vcf}
            rm -rf {resources.tmpdir}
            ##
        """

ruleorder: remove_untypable > convert_genotyping_to_biallelic >  FilterAndBeagle > TheGreatGenotyperSecondPassHMM > TheGreatGenotyperHMM >  TheGreatGenotyperkmersOnly > bgzip


rule merge_population_chunks:
    input:
        population = expand("{{results}}/{{sample}}/{{sample_test}}/{{tool}}/population.chunk-{chunk}.vcf.gz", chunk= range(1,11)),
    output:
        "{results}/{sample}/{sample_test}/{tool}/population.vcf.gz" 
    log:
        "{results}/{sample}/{sample_test}/{tool}/population.log"
    threads: 16
    retries: 0
    resources:
        mem_mb=20*1024,
        cores=16,
        nodes = 1,
        meduim=1,
        runtime = 60 * 4,
        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}"),
        tmp= lambda wildcards: "%smerge%s_%s/"%(tempFolder,f"{wildcards.sample}",f"{wildcards.tool}")
    shell:
        r"""
            mkdir -p {resources.tmp}
            bcftools merge {input} | parallel -j {threads} --pipe --block 10m -k /home/mshokrof/TheGreatGenotyper/build6/PopMergeVCF/PopMergeVCF 2> {log} | bgzip > {resources.tmp}tmp.vcf.gz 
            mv {resources.tmp}tmp.vcf.gz {output}
            rm -rf {resources.tmp}
        #
          """


rule PreProcessPanel:
    input:
        population = "{results}/{sample}/{sample_test}/{tool}/population.vcf.gz" 
    output:
        "{results}/{sample}/{sample_test}/{tool}/population.preprocessed.vcf.gz"  
    log:
        "{results}/{sample}/{sample_test}/{tool}/population.preprocessed.log"
    threads: 64
    retries: 0
    params:
        out_prefix="{results}/{sample}/{sample_test}/{tool}/population.preprocessed"
    resources:
        mem_mb=80*1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 10,
        #partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}"),
        partition =  "bmm",
        tmp= lambda wildcards: "%sPreProcessPanel%s_%s/"%(tempFolder,f"{wildcards.sample}",f"{wildcards.tool}")
    shell:
        r"""
            mkdir -p {resources.tmp}
            bcftools view --header-only {input.population} | zgrep -P "^#CHROM"  |cut -f10- | tr -s $'\t' $'\n'| grep -P ".*_TheGreatGenotyper$"  > {resources.tmp}tmp.$$
            bcftools view -S ^{resources.tmp}tmp.$$ {input.population} 2> {log} |bgzip -c > {resources.tmp}tmp.$$.vcf.gz
            tabix -p vcf {resources.tmp}tmp.$$.vcf.gz 
           # /home/mshokrof/miniconda3/envs/cattle_sv/bin/plink2 --vcf {resources.tmp}tmp.$$.vcf.gz --het --missing --vcf-half-call missing --out {resources.tmp}plink 2>> {log}
           # tail -n +2 {resources.tmp}plink.het |awk '$5 < 0.45' |cut -f1 > {resources.tmp}samples.het0.2.txt
           # bcftools view -S {resources.tmp}samples.het0.2.txt {resources.tmp}tmp.$$.vcf.gz 2>> {log} |bgzip -c > {resources.tmp}tmp.2.$$.vcf.gz
           # tabix -p vcf {resources.tmp}tmp.2.$$.vcf.gz
           # /home/mshokrof/miniconda3/envs/cattle_sv/bin/plink2 --vcf {resources.tmp}tmp.2.$$.vcf.gz  --hwe 1e-50 'midp' --mind 0.2  --recode vcf --vcf-half-call missing  --output-chr 'chrM' --out {resources.tmp}tmp.3.$$ 2>> {log}
            java -Xmx70G -jar {beagle} gt={resources.tmp}tmp.$$.vcf.gz    out={params.out_prefix} nthreads={threads}  map={beagleMap} &> {log}
            rm -rf {resources.tmp}
          """
ruleorder: pangenie > FilterAndBeagle
rule FilterAndBeagle:
    input:
        population = "{results}/{sample}/{sample_test}/{tool}/population.vcf.gz" ,
        reference_panel = "{results}/{sample}/{sample_test}/{tool}/population.preprocessed.vcf.gz"
    output:
        "{results}/{sample}/{sample_test}/{tool}/{sample}_{fraction}_{tool}_genotyping.vcf" 
    log:
        "{results}/{sample}/{sample_test}/{tool}/{sample}_{fraction}_{tool}.log"
    threads: 64
    retries: 0
    params:
        out_prefix="{results}/{sample}/{sample_test}/{tool}/{sample}_{fraction}_{tool}_genotyping",
        sample_name= "{sample}_{fraction}_TheGreatGenotyper"
    wildcard_constraints:
        tool="TheGreatGenotyper-SecondPassHMM|TheGreatGenotyper-kmersOnly|TheGreatGenotyper-HMM"
    resources:
        mem_mb=80*1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 10,
#        partition =  lambda wildcards: getMeduimPartition(f"{wildcards.sample}"),
        partition =  "bmm",
        tmp= lambda wildcards: "%sbeagle%s_%s_%s/"%(tempFolder,f"{wildcards.sample}",f"{wildcards.tool}",f"{wildcards.fraction}")
    shell:
        r"""
            mkdir -p {resources.tmp}
            bcftools view -s {wildcards.sample}_{wildcards.fraction}_TheGreatGenotyper {input.population} 2> {log} |bgzip -c > {resources.tmp}tmp.$$.vcf.gz
            tabix -p vcf {resources.tmp}tmp.$$.vcf.gz 
            java -Xmx70G -jar {beagle} gt={resources.tmp}tmp.$$.vcf.gz ref={input.reference_panel}    out={resources.tmp}tmp4.$$ nthreads={threads}  map={beagleMap} 2>&1 >> {log}
            tabix -p vcf {resources.tmp}tmp4.$$.vcf.gz
            python ../scripts/unphase.py {resources.tmp}tmp4.$$.vcf.gz {output} 2>> {log}
            rm -rf {resources.tmp}
          """


def getBAM(wildcards):
    return config['data'][wildcards.sample]["mapping_prefix"] + wildcards.sample + '-' +\
 wildcards.fraction+ ".chr" +wildcards.chrom+ '.bam'


#########################################################
#################    run GraphTyper    ##################
#########################################################

# input VCF for graphtyper (biallelic, untypables removed since representation is changed)
rule graphtyper_preprocess:
    input:
        vcf="{results}/{sample}/{sample_test}/biallelic/{sample}-chr{chrom}.vcf.gz",
        tbi="{results}/{sample}/{sample_test}/biallelic/{sample}-chr{chrom}.vcf.gz.tbi"
    output:
        "{results}/{sample}/{sample_test}/graphtyper/variants/variants_chr{chrom}_{variant}.vcf.gz"
    wildcard_constraints:
        chrom="X|Y|[0-9]+",
        variant="indel|sv"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        """
        zcat {input.vcf} | python3 ../scripts/extract-varianttype.py  {wildcards.variant} | bgzip > {output}
#        tabix -p vcf {output}

        #
        """

# genotype SNPs/indels/SVs (/usr/bin/time -v {graphtyper} genotype {input.fasta} --vcf={input.vcf_small} --sam={input.bam} --region=chr{wildcards.chrom} --no_decompose --verbose --output={params.dir_small} --threads={threads}) &> {log.small} 
rule graphtyper_genotype:
    input:
        #    vcf_small="{results}/{sample}/{sample_test}/graphtyper/variants/variants_chr{chrom}_indel.vcf.gz",
        #    tbi_small="{results}/{sample}/{sample_test}/graphtyper/variants/variants_chr{chrom}_indel.vcf.gz.tbi",
        vcf_sv="{results}/{sample}/{sample_test}/graphtyper/variants/variants_chr{chrom}_sv.vcf.gz",
        tbi_sv="{results}/{sample}/{sample_test}/graphtyper/variants/variants_chr{chrom}_sv.vcf.gz.tbi",
        bam= lambda wildcards: getBAM(wildcards),
        bai= lambda wildcards: getBAM(wildcards) +'.bai',
        fasta=input_reference
    output:
        "{results}/{sample}/{sample_test}/graphtyper/{sample}_{fraction}_graphtyper_genotyping.chr{chrom, X|Y|[0-9]+}.vcf"
    params:
        #        dir_small="{results}/{sample}/{sample_test}/graphtyper/{sample}_{fraction}_{chrom}_small",
        dir_sv="{results}/{sample}/{sample_test}/graphtyper/{sample}_{fraction}_{chrom}_sv"
    log:
        #        small="{results}/{sample}/{sample_test}/graphtyper/{sample}_{fraction}_graphtyper_discovery.chr{chrom}.small.log",
        sv="{results}/{sample}/{sample_test}/graphtyper/{sample}_{fraction}_graphtyper_discovery.chr{chrom}.sv.log"
    conda:
        "env.yaml"
    threads:
        24
    resources:
        mem_mb=30000,
        runtime_hrs=4,
        runtime_min=59,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        """
        (/usr/bin/time -v {graphtyper} genotype_sv {input.fasta} {input.vcf_sv} --sam={input.bam} --region=chr{wildcards.chrom} --output={params.dir_sv} --threads={threads}) &> {log.sv}
        bcftools concat -a {params.dir_sv}/chr{wildcards.chrom}/*.vcf.gz | python3 ../scripts/graphtyper-postprocess.py {input.vcf_sv} > {output}
        #
        """


########################################################
##################      Evaluation    ##################
########################################################

def get_bed(wildcards):
    if wildcards.regions == "repeats-complex":
        return "{results}/bed/complex-rep.bed"
    elif wildcards.regions == "nonrep-complex":
        return "{results}/bed/complex-nonrep.bed"
    elif wildcards.regions == "repeats-simple":
        return "{results}/bed/simple-rep.bed"
    elif wildcards.regions == "nonrep-simple":
        return "{results}/bed/simple-nonrep.bed"
    elif wildcards.regions == "external":
        return config['data'][wildcards.sample]['external_bed']
    else:
        assert (False)


# prepare all regions
rule prepare_regions:
    input:
        fai=input_reference + '.fai',
        repeats=repeats_bed,
        complex=complex_bed
    output:
        tmp1=temp("{results}/bed/tmp1.txt"),
        tmp2=temp("{results}/bed/tmp2.txt"),
        tmp3=temp("{results}/bed/tmp4.txt"),
        complex_nonrep="{results}/bed/complex-nonrep.bed",
        complex_rep="{results}/bed/complex-rep.bed",
        simple_nonrep="{results}/bed/simple-nonrep.bed",
        simple_rep="{results}/bed/simple-rep.bed"
    conda:
        "../env/genotyping.yml"
    shell:
        "bash ../scripts/prepare-beds.sh {input.repeats} {input.complex} {input.fai} {output.tmp1} {output.tmp2} {output.tmp3} {output.complex_rep} {output.simple_rep} {output.complex_nonrep} {output.simple_nonrep}"


# prepare evaluation region (graph)
rule evaluation_region_graph:
    input:
        regions=get_bed,
        callable=lambda wildcards: config['data'][wildcards.sample]['bed']
    output:
        "{results}/{sample}/{sample_test}/bed/graph/{sample}_{regions}.bed"
    wildcard_constraints:
        regions="repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        "bedtools intersect -a {input.regions} -b {input.callable} > {output}"


# prepare evaluation region (external)
rule evaluation_region_external:
    input:
        regions=get_bed,
        callable="{results}/{sample}/{sample_test}/bed/graph/{sample}_external.bed"
    output:
        "{results}/{sample}/{sample_test}/bed/external/{sample}_{regions}.bed"
    wildcard_constraints:
        regions="repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        "bedtools intersect -a {input.regions} -b {input.callable} > {output}"

# annotate genotyped VCFs and turn them to bi-allelic
rule convert_genotyping_to_biallelic:
    input:
        callset="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_genotyping.vcf",
        graph=lambda wildcards: "{results}/{sample}/{sample_test}/biallelic/{sample}-all.vcf" if wildcards.method in ["platypus",
                                                                                                        "gatk",
                                                                                                        "graphtyper"] else "{results}/{sample}/{sample_test}/reference-panel/{sample}-all.vcf",
        biallelic= "{results}/{sample}/{sample_test}/truth/no{sample}_biallelic.vcf"
    output:
        "{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_genotyping-biallelic.vcf"
    wildcard_constraints:
        method= "|".join(other_methods+["pangenie"])
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=20000,
        runtime_hrs=2,
        runtime=2 * 60,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    shell:
        """
        cat {input.callset} | python3 ../scripts/annotate.py {input.graph} | python3 ../scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' >  {output}

        """


# # normalize and annotate discovery sets and assign variant IDs
# rule convert_discovery_to_biallelic:
#     input:
#         callset="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_discovery.vcf",
#         biallelic=config['data']['full_callset'],
#         reference=input_reference
#     output:
#         "{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_discovery-biallelic.vcf.gz"
#     conda:
#         "../env/genotyping.yml"
#     resources:
#         mem_mb=20000,
#         runtime_hrs=2,
#         runtime=2 * 60
#     shell:
#         """
#         cat {input.callset} | python3 ../scripts/annotate.py <(gzip -dc {input.biallelic}) | bgzip > {output}
# #        tabix -p vcf {output}
#         """

# # prepare external ground truth sets
# rule prepare_ground_truth:
#     input:
#         callset=lambda wildcards: config['data'][wildcards.sample]['external'],
#         truth=config['data']['full_callset'],
#         reference=input_reference
#     output:
#         vcf="{results}/{sample}/{sample_test}/external-truth/{sample}-truth.vcf.gz"
#     conda:
#         "../env/genotyping.yml"
#     resources:
#         mem_mb=20000,
#         runtime_hrs=2,
#         runtime=2 * 60
#     shell:
#         """
#         bcftools norm -f {input.reference} -m -any {input.callset} | bcftools sort | python3 ../scripts/annotate.py <(gzip -dc {input.truth}) | bgzip > {output.vcf}
#         tabix -p vcf {output}
#         """
print(samples)
# determine untypable IDs
rule remove_untypable:
    input:
        vcf="{results}/{sample}/{sample_test}/{prefix}{sample}{other}.vcf",
        ids=lambda wildcards: "{results}/{sample}/{sample_test}/truth/untypable-ids/{sample}-untypable.tsv" if wildcards.mode == "graph" else
        config['data'][wildcards.sample]['untypable_external']
    output:
        "{results}/{sample}/{sample_test}/{prefix}{sample}{other}-typable-{vartype}-{mode}.vcf"
    wildcard_constraints:
        sample="|".join(samples),
        sample_test="|".join(samples),
        vartype="|".join(variants),
        mode="external|graph"
    resources:
        mem_mb=20000,
        runtime_hrs=1,
        runtime=60
    shell:
        """
        cat {input.vcf} | python3 ../scripts/skip-untypable.py {input.ids} | python3 ../scripts/extract-varianttype.py  {wildcards.vartype} >{output}


        """

rule prepare_all:
    input:
        vcf="{path}{sample}{other}.vcf.gz",
        tbi="{path}{sample}{other}.vcf.gz.tbi"
    output:
        "{path}{sample}{other}-all-{vartype}-{mode}.vcf"
    wildcard_constraints:
        sample="|".join(samples),
        vartype="|".join(variants),
        mode="external|graph"
    resources:
        mem_mb=20000,
        runtime_hrs=1,
        runtime=60
    shell:
        """
        zcat {input.vcf} | python3 ../scripts/extract-varianttype.py  {wildcards.vartype}  > {output}


        """


########################## compute precision/recall ##############################

rule rtg_format:
    input:
        input_reference
    output:
        directory("{results}/evaluation/SDF")
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"SDF", attempt),
    shell:
        'rtg format -o {output} {input}'


# precision-recall
rule vcfeval:
    input:
    #results/HG00731/meduim-NA12878-ASSM/evaluation/precision-recall-typable/graph/TheGreatGenotyper-HMM-genotyping-biallelic/coverage-5_repeats-complex_small-deletion/qual_0/summary.txt
    #results/HG00731/meduim-NA12878-ASSM/evaluation/HG00731_5_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz
        callset="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz",
        callset_tbi="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz.tbi",
        baseline=lambda wildcards: "{results}/{sample}/{sample_test}/external-truth/{sample}-truth-{variantset}-{vartype}-{mode}.vcf.gz" if wildcards.mode == "external" else
        "{results}/{sample}/{sample_test}/truth/only{sample}_biallelic-{variantset}-{vartype}-{mode}.vcf.gz",
        baseline_tbi=lambda wildcards: "{results}/{sample}/{sample_test}/external-truth/{sample}-truth-{variantset}-{vartype}-{mode}.vcf.gz.tbi" if wildcards.mode == "external" else
        "{results}/{sample}/{sample_test}/truth/only{sample}_biallelic-{variantset}-{vartype}-{mode}.vcf.gz.tbi",
        #        regions=lambda wildcards: "{results}/{sample}/{sample_test}/bed/{sample}_{regions}.bed" if wildcards.mode == "graph" else config['data'][wildcards.sample]['external_bed'],
        regions="{results}/{sample}/{sample_test}/bed/{mode}/{sample}_{regions}.bed",
        sdf="{results}/evaluation/SDF"
    output:
        summary="{results}/{sample}/{sample_test}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_0/summary.txt"
      #  results/HG00731/meduim-NA12878-ASSM/evaluation/precision-recall-typable/graph/pangenie-genotyping-biallelic/coverage-5_repeats-complex_small-deletion/qual_0/summary.txt
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        sample="|".join(samples),
        mode="external|graph",
        run_mode="genotyping-biallelic|discovery-biallelic",
       # fraction="|".join([str(f) for f in downsampling]),
        regions="repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external",
        vartype="|".join(variants),
        variantset="typable|all"
    params:
        tmp="{results}/{sample}/{sample_test}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_0_tmp",
        outname="{results}/{sample}/{sample_test}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_0",
        which=lambda wildcards: "--all-records" if wildcards.run_mode == 'genotyping-biallelic' else ""
    resources:
        mem_mb=10000,
        runtime_hrs=0,
        runtime_min=20,
        runtime=20
    shell:
        """
        rm -rf {params.tmp}
        rtg vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 > {output.summary}.tmp
        mv {params.tmp}/* {params.outname}/
        mv {output.summary}.tmp {output.summary}
        rm -r {params.tmp}
        """


# precision-recall quality cutoff
rule vcfeval_cutoff:
    input:
        callset="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz",
        callset_tbi="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz.tbi",
        baseline=lambda wildcards: "{results}/{sample}/{sample_test}/external-truth/{sample}-truth-{variantset}-{vartype}-{mode}.vcf.gz" if wildcards.mode == "external" else
        "{results}/{sample}/{sample_test}/truth/only{sample}_biallelic-{variantset}-{vartype}-{mode}.vcf.gz",
        baseline_tbi=lambda wildcards: "{results}/{sample}/{sample_test}/external-truth/{sample}-truth-{variantset}-{vartype}-{mode}.vcf.gz.tbi" if wildcards.mode == "external" else
        "{results}/{sample}/{sample_test}/truth/only{sample}_biallelic-{variantset}-{vartype}-{mode}.vcf.gz.tbi",
        #        regions=lambda wildcards: "{results}/{sample}/{sample_test}/bed/{sample}_{regions}.bed" if wildcards.mode == "graph" else config['data'][wildcards.sample]['external_bed'],
        regions="{results}/{sample}/{sample_test}/bed/{mode}/{sample}_{regions}.bed",
        sdf="{results}/evaluation/SDF"
    output:
        tmp_vcf=temp("{results}/{sample}/{sample_test}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}_coverage-{fraction}_{regions}_{vartype}_qual_{qual}.vcf.gz"),
        tmp_tbi=temp("{results}/{sample}/{sample_test}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}_coverage-{fraction}_{regions}_{vartype}_qual_{qual}.vcf.gz.tbi"),
        summary="{results}/{sample}/{sample_test}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt"
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        sample="|".join(samples),
        mode="external|graph",
        run_mode="genotyping-biallelic|discovery-biallelic",
       # fraction="|".join([str(f) for f in downsampling]),
        regions="repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external",
        vartype="|".join(variants),
        qual="200",
        variantset="typable|all"
    params:
        tmp="{results}/{sample}/{sample_test}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}_tmp",
        outname="{results}/{sample}/{sample_test}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}",
        which=lambda wildcards: "--all-records" if wildcards.run_mode == 'genotyping-biallelic' else ""
    resources:
        mem_mb=10000,
        runtime_hrs=0,
        runtime_min=20,
        runtime=20
    shell:
        """
        bcftools view -i 'FMT/GQ>={wildcards.qual}' -O z {input.callset} > {output.tmp_vcf}
        tabix -p vcf {output.tmp_vcf}
        rtg vcfeval -b {input.baseline} -c {output.tmp_vcf} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 > {output.summary}.tmp
        mv {params.tmp}/* {params.outname}/
        mv {output.summary}.tmp {output.summary}
        rm -r {params.tmp}
        """


################################## compute genotype concordance ########################################


# determine the variants that went into re-typing per category
rule collected_typed_variants:
    input:
        graph=lambda wildcards: "{results}/{sample}/{sample_test}/truth/no{sample}_biallelic.vcf.gz",
        #        regions = lambda wildcards: "{results}/{sample}/{sample_test}/bed/{sample}_{regions}.bed" if wildcards.mode == "graph" else config['data'][wildcards.sample]['external_bed']
        regions="{results}/{sample}/{sample_test}/bed/{mode}/{sample}_{regions}.bed"
    output:
        "{results}/{sample}/{sample_test}/genotyped-ids/{mode}_{regions}_{vartype}.tsv"
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        sample="|".join(samples),
        mode="external|graph",
        regions="repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external",
        vartype="|".join(variants)
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt)
    shell:
        "zcat {input.graph} | python3 ../scripts/extract-varianttype.py  {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 ../scripts/get_ids.py > {output}    "


# compute concordances
rule genotype_concordances:
    input:
        callset="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}-typable-{vartype}-{mode}.vcf.gz",
        callset_tbi="{results}/{sample}/{sample_test}/{method}/{sample}_{fraction}_{method}_{run_mode}-typable-{vartype}-{mode}.vcf.gz.tbi",
        baseline=lambda wildcards: "{results}/{sample}/{sample_test}/external-truth/{sample}-truth-typable-{vartype}-{mode}.vcf.gz" if wildcards.mode == "external" else
        "{results}/{sample}/{sample_test}/truth/only{sample}_biallelic-typable-{vartype}-{mode}.vcf.gz",
        baseline_tbi=lambda wildcards: "{results}/{sample}/{sample_test}/external-truth/{sample}-truth-typable-{vartype}-{mode}.vcf.gz.tbi" if wildcards.mode == "external" else
        "{results}/{sample}/{sample_test}/truth/only{sample}_biallelic-typable-{vartype}-{mode}.vcf.gz.tbi",
        #        regions=lambda wildcards: "{results}/{sample}/{sample_test}/bed/{sample}_{regions}.bed" if wildcards.mode == "graph" else config['data'][wildcards.sample]['external_bed'],
        regions="{results}/{sample}/{sample_test}/bed/{mode}/{sample}_{regions}.bed",
        typed_ids="{results}/{sample}/{sample_test}/genotyped-ids/{mode}_{regions}_{vartype}.tsv"
    output:
        tmp_vcf1=temp("{results}/{sample}/{sample_test}/evaluation/concordance/{mode}/{method}-{run_mode}_coverage-{fraction}_{regions}_{vartype}_qual_{qual}_base.vcf"),
        tmp_vcf2=temp("{results}/{sample}/{sample_test}/evaluation/concordance/{mode}/{method}-{run_mode}_coverage-{fraction}_{regions}_{vartype}_qual_{qual}_call.vcf"),
        summary="{results}/{sample}/{sample_test}/evaluation/concordance/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt"
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        sample="|".join(samples),
        mode="external|graph",
        run_mode="genotyping-biallelic|discovery-biallelic",
    #    fraction="|".join([str(f) for f in downsampling]),
        regions="repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external",
        vartype="|".join(variants),
        qual="0|200"
    params:
        which=lambda wildcards: "" if wildcards.run_mode == 'genotyping-biallelic' else "--only_pass"
    log:
        "{results}/{sample}/{sample_test}/evaluation/concordance/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.log"
    resources:
        mem_mb=20000,
        runtime_hrs=0,
        runtime_min=20,
        runtime=20
    shell:
        """
        bedtools intersect -header -a {input.baseline} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf1}
        bedtools intersect -header -a {input.callset} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf2}
        python3 ../scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual {wildcards.qual} {params.which} 2> {log} 1> {output.summary}
        """


########################################################
##################       Plotting     ##################
########################################################

def plot_input(wildcards):
    output = []
    variant_list = ['indel', 'sv', 'large-insertion', 'large-deletion', 'large'] if (
                                                                                                wildcards.run_mode == "discovery-biallelic") or (
                                                                                                wildcards.mode != "graph") else variants
    all_regions = [wildcards.regions] if 'external' in wildcards.regions else [wildcards.regions + '-simple',
                                                                               wildcards.regions + '-complex']
    metric = 'precision-recall-typable' if wildcards.metric == 'fscore' else wildcards.metric
    for reg in all_regions:
        if wildcards.run_mode == "genotyping-biallelic":
            for fraction in downsampling:
                for var in variant_list:
                    for method in ['pangenie'] + other_methods:
                        output.append("{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{method}-genotyping-biallelic/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt".format(
                            results=wildcards.results,
                            sample_test=wildcards.sample_test,
                            vartype=var,
                            sample=wildcards.sample,
                            mode=wildcards.mode,
                            method=method,
                            fraction=fraction,
                            regions=reg,
                            qual="0",
                            metric=metric
                        )
                        )
            #         output.append("{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{method}-genotyping-biallelic/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt".format(
            #             results=wildcards.results,
            #             sample_test=wildcards.sample_test,
            #             vartype=var,
            #             sample=wildcards.sample,
            #             mode=wildcards.mode,
            #             fraction=fraction,
            #             regions=reg,
            #             qual="200",
            #             metric=metric,
            # method="pangenie"
            #         )
            #)
            # output.append("{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{method}-genotyping-biallelic/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt".format(
                    #     results=wildcards.results,
                    #     sample_test=wildcards.sample_test,
                    #     vartype=var,
                    #     sample=wildcards.sample,
                    #     mode=wildcards.mode,
                    #     fraction=fraction,
                    #     regions=reg,
                    #     qual="200",
                    #     metric=metric,
            #     method="TheGreatGenotyper"
                    # )

            #        )
        elif wildcards.run_mode == "discovery-biallelic":
            for fraction in downsampling:
                for method in ['gatk', 'platypus']:
                    for var in variant_list:
                        output.append("{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{method}-discovery-biallelic/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt".format(
                            results=wildcards.results,
                            sample_test=wildcards.sample_test,
                            sample=wildcards.sample,
                            vartype=var,
                            mode=wildcards.mode,
                            method=method,
                            fraction=fraction,
                            regions=reg,
                            qual="0",
                            metric=metric
                        )
                        )
        else:
            assert (False)
    #print( output)
    return output

# plot precision/recall or genotype concordance over all coverages
rule plot_results:
    input:
        plot_input
    output:
        "{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{run_mode}_{regions}_{metric}.pdf"
    log:
        "{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{run_mode}_{regions}_{metric}.log"
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        run_mode="genotyping-biallelic|discovery-biallelic",
        regions="repeats|nonrep|external",
        metric="precision-recall-typable|precision-recall-all|concordance|fscore"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    params:
        folder="{results}/{sample}/{sample_test}",
        script=lambda wildcards: metric_to_script[wildcards.metric],
        coverages=lambda wildcards: [max(downsampling)] if wildcards.metric == "fscore" else downsampling,
        run_mode=lambda wildcards: wildcards.run_mode if 'genotyping-biallelic' == wildcards.run_mode else [
            wildcards.run_mode, 'genotyping-biallelic'],
        varset=lambda wildcards: "-variantset all" if wildcards.metric == "precision-recall-all" else ""
    shell:
        "python3  {params.script}   {wildcards.mode} {wildcards.regions} -run_mode {params.run_mode} -folder {params.folder} -coverages {params.coverages} -outfile {output} -sample {wildcards.sample} {params.varset} &> {log}"



# plot precision/recall or genotype concordance over all coverages
rule plot_results_fscore:
    input:
        plot_input
    output:
        figure="{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{run_mode}_{regions}_{metric}_{coverage}.pdf",
        summary="{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{run_mode}_{regions}_{metric}_{coverage}.csv"
    log:
        "{results}/{sample}/{sample_test}/evaluation/{metric}/{mode}/{run_mode}_{regions}_{metric}_{coverage}.log"
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        run_mode="genotyping-biallelic|discovery-biallelic",
        regions="repeats|nonrep|external",
        metric="precision-recall-typable|precision-recall-all|concordance|fscore"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    params:
        folder="{results}/{sample}/{sample_test}",
        script=lambda wildcards: metric_to_script[wildcards.metric],
        run_mode=lambda wildcards: wildcards.run_mode if 'genotyping-biallelic' == wildcards.run_mode else [
            wildcards.run_mode, 'genotyping-biallelic'],
        varset=lambda wildcards: "-variantset all" if wildcards.metric == "precision-recall-all" else ""
    shell:
        "python3 {params.script} {wildcards.mode} {wildcards.regions} -run_mode {params.run_mode} -folder {params.folder} -coverages {wildcards.coverage} -outfile {output.figure} -sample {wildcards.sample} {params.varset} 2> {log} > {output.summary}"




rule plot_newFscore:
    input:
        summary_repeats="{results}/{sample}/{sample_test}/evaluation/fscore/{mode}/{run_mode}_repeats_fscore_{coverage}.csv",
        summary_nonrepeats="{results}/{sample}/{sample_test}/evaluation/fscore/{mode}/{run_mode}_nonrep_fscore_{coverage}.csv"
    output:
        figure= "{results}/{sample}/{sample_test}/evaluation/newFscore/{mode}/{run_mode}_fscore_{coverage}.png"    
    log:
        "{results}/{sample}/{sample_test}/evaluation/newFscore/{mode}/{run_mode}_fscore_{coverage}.log"
    conda:
        "../env/genotyping.yml"
    wildcard_constraints:
        run_mode="genotyping-biallelic|discovery-biallelic",
        regions="repeats|nonrep|external"
    resources:
        mem_mb=4000,
        runtime_hrs=2,
        runtime_min=120,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
    params:
        script=lambda wildcards: metric_to_script["newFscore"],
    shell:
        "python3 {params.script} {input.summary_repeats} {input.summary_nonrepeats} {output.figure} &> {log}"
