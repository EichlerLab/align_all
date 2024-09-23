import pandas as pd
import os

configfile: 'map_contigs.yaml'

MANIFEST = config['MANIFEST']
REF_DICT = config['REF']
ALN_PARAMS = config.get('ALN_PARAMS', '')
NBATCHES = config.get('NBATCHES', 15)

command_dict = {}
command_dict['CLR'] = "pbmm2 align --preset SUBREAD -j"
command_dict['CCS'] = "pbmm2 align --preset CCS -j"
command_dict['PacBio_HiFi'] = "pbmm2 align --preset HiFi -j"
command_dict['REVIO'] = "pbmm2 align --preset CCS -j"
command_dict['ONT'] = 'minimap2 -ax map-ont -I 8G -t'
command_dict['ONT_UL'] = 'minimap2 -ax map-ont -I 8G -t'
command_dict['Illumina'] = 'bwa mem -Y -K 100000000 -t'

manifest_df = pd.read_csv(MANIFEST, sep='\s+', header=0, dtype=str)
manifest_df = manifest_df.set_index(['SAMPLE', 'TYPE'], drop=False)
manifest_df.sort_index(inplace=True)

def find_read(wildcards):
    read_df = pd.read_csv(manifest_df.at[(wildcards.sample, wildcards.aln), 'FOFN'], header=None, sep='\t')
    return read_df.at[int(wildcards.read), 0].split(" ")

def find_read_batch(wildcards):
    read_df = pd.read_csv(manifest_df.at[(wildcards.sample, wildcards.aln), 'FOFN'], header=None, sep='\t', names=['FILE'])
    file = read_df.iloc[int(wildcards.read)]['FILE']
    return f'{file}.fai'

def combine_reads(wildcards):
    read_df = pd.read_csv(manifest_df.at[(wildcards.sample, wildcards.aln),'FOFN'], header=None, sep='\t')
    if wildcards.aln in ['ONT', 'ONT_UL', 'REVIO']:
        return expand(gather.split('tmp/{{ref}}/{{aln}}/{{sample}}.{{read}}_{scatteritem}.sorted.bam'), sample=wildcards.sample, read=read_df.index, ref=wildcards.ref, aln=wildcards.aln)
    else:
        return expand('{ref}/{aln}/{sample}.{read}.sorted.bam', sample=wildcards.sample, read=read_df.index, ref=wildcards.ref, aln=wildcards.aln)

def find_map(wildcards):
    if wildcards.aln == 'Illumina':
        return rules.index_ref.output.ref
    else:
        return REF_DICT[wildcards.ref]

scattergather:
    split=NBATCHES

def find_command(wildcards):
    return command_dict[manifest_df.at[(wildcards.sample, wildcards.aln), 'TYPE']]

def find_ref(wildcards):
    return REF_DICT[wildcards.ref]

def find_aln_params(wildcards):
    if wildcards.aln in ['CCS', 'PacBio_HiFi', 'CLR', 'REVIO']:
        return ALN_PARAMS+" "+f"--sample {wildcards.sample} --rg '@RG\\tID:{wildcards.read}'"
    elif wildcards.aln in ['ONT', 'ONT_UL']:
        if wildcards.aln == 'ONT':
            library="STD"
        else:
            library="UL"
        return ALN_PARAMS+" "+f"-R '@RG\\tID:{wildcards.read}\\tSM:{wildcards.sample}\\tPL:ONT\\tLB:{library}'"
    elif wildcards.aln in ['Illumina']:
        return ALN_PARAMS+" "+f"-R '@RG\\tID:{wildcards.read}\\tSM:{wildcards.sample}\\tPL:ILLUMINA'"
    else:
        return ALN_PARAMS

wildcard_constraints:
    sample='|'.join(manifest_df['SAMPLE'].unique()),
    aln='|'.join(command_dict),
    read='\d+',
    ref='|'.join(REF_DICT)


localrules: all, index_ref

rule all:
    input:
        expand(expand('{{ref}}/{aln}/{sample}.all.sorted.bam', zip, sample=manifest_df.index.get_level_values('SAMPLE'), aln=manifest_df.index.get_level_values('TYPE')), ref=REF_DICT)

rule index_ref:
    input:
        ref = find_ref
    output:
        ref = 'ref/{ref}.fa',
        index = 'ref/{ref}.fa.fai',
        bwa_index = 'ref/{ref}.fa.amb'
    resources:
        mem = 12,
        smem = 4,
    threads: 1
    singularity:
        "docker://eichlerlab/align-basics:0.2",
    shell:
        '''
        ln -s $( readlink -f {input.ref} ) {output.ref}
        samtools faidx {output.ref}
        bwa index {output.ref}
        '''

rule map_reads:
    input:
        ref = find_map,
        read = find_read
    output:
        bam = temp('{ref}/{aln}/{sample}.{read}.bam')
    resources:
        mem = 12,
        smem = 4,
        hrs = 96
    threads: 8
    params:
        command = find_command,
        aln_params = find_aln_params
    singularity:
        "docker://eichlerlab/align-basics:0.2",
    shell:
        '''
        {params.command} {threads} {params.aln_params} {input.ref} {input.read} | samtools view -b - > {output.bam}
        '''


rule get_batch_ids:
    input:
        fai=find_read_batch
    output:
        batches=temp(
            scatter.split(
                "tmp/splitBatchID/{{sample}}/{{aln}}/{{read}}_{scatteritem}.txt"
            )
        ),
    resources:
        mem=4,
        hrs=5,
    threads: 1
    script:
        "smk_scripts/split_faidx.py"


rule map_split:
    input:
        fastq=find_read,
        batch_file="tmp/splitBatchID/{sample}/{aln}/{read}_{scatteritem}.txt",
        ref=find_map,
    output:
        bam=temp(
            'tmp/{ref}/{aln}/{sample}.{read}_{scatteritem}.sorted.bam'
        )
    resources:
        mem=lambda wildcards, attempt: min(4**attempt,128),
        hrs=96,
    params:
        command = find_command,
        aln_params = find_aln_params
    threads: 4
    singularity:
        "docker://eichlerlab/align-basics:0.2",
    shell:
        """
        samtools fqidx {input.fastq} -r {input.batch_file} > {resources.tmpdir}/{wildcards.scatteritem}.fastq
        {params.command} {threads} {params.aln_params} {input.ref} {resources.tmpdir}/{wildcards.scatteritem}.fastq | samtools view -b - | sambamba sort -t {threads} -o {output.bam} -m 30G /dev/stdin
        rm {resources.tmpdir}/{wildcards.scatteritem}.fastq
        """


rule sort_indiv:
    input:
        bam = rules.map_reads.output.bam
    output:
        sort_bam = temp('{ref}/{aln}/{sample}.{read}.sorted.bam'),
        index = temp('{ref}/{aln}/{sample}.{read}.sorted.bam.bai')
    resources:
        mem = 18,
        smem = 4,
        hrs = 8
    threads: 12
    singularity:
        "docker://eichlerlab/align-basics:0.2",
    shell:
        '''
        sambamba sort -t {threads} -o {output.sort_bam} -m 125G {input.bam}
        '''

rule merge_maps:
    input:
        reads = combine_reads
    output:
        merged = '{ref}/{aln}/{sample}.all.sorted.bam',
        index = '{ref}/{aln}/{sample}.all.sorted.bam.bai'
    resources:
        mem = 12,
        smem = 4,
        hrs = 12
    threads: 8
    singularity:
        "docker://eichlerlab/align-basics:0.2",
    shell:
        '''
        if [[ $( echo {input.reads} | wc -w ) == 1 ]]; then
            cp -rl {input.reads} {output.merged}; samtools index {output.merged}
        else
            samtools merge -@ {threads} -o {output.merged} {input.reads}; samtools index {output.merged}
        fi
        '''


rule mark_duplicates:
    input:
        bam = rules.merge_maps.output.merged,
        index = rules.merge_maps.output.index
    output:
        bam = '{ref}/{aln}/{sample}.all.sorted.md.bam',
        index = '{ref}/{aln}/{sample}.all.sorted.md.bam.bai'
    resources:
        mem = 4,
        smem = 4,
        hrs = 48
    threads: 8
    singularity:
        "docker://eichlerlab/align-basics:0.2",
    shell:
        '''
        sambamba markdup -t {threads} --tmpdir={resources.tmpdir} {input.bam} {output.bam}
        samtools index {output.bam}
        '''

rule cram_conv:
    input:
        bam = rules.mark_duplicates.output.bam,
        ref = find_ref
    output:
        cram = '{ref}/{aln}/cram/{sample}.final.cram',
        crai = '{ref}/{aln}/cram/{sample}.final.cram.crai'
    resources:
        mem = 4,
        smem = 4,
        hrs = 48
    singularity:
        "docker://eichlerlab/align-basics:0.2",
    threads: 1
    shell:
        '''
        samtools view -C -T {input.ref} {input.bam} > {output.cram}
        samtools index {output.cram} 
        '''
        

