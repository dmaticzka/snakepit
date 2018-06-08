configfile: 'config.yaml'

include: 'rules/get_input_globs.smk'
include: 'rules/helpers.smk'
include: 'rules/targetdist.smk'

rule all:
    message: 'Help Text Here'

def make_peakachu_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/peakachu/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_peakachu.bed', fn), chdir))
    return(chsuff)

def make_peakachu_bam_conversion(wildcards):
    signal = get_input_bed("input/{}/signal".format(wildcards.id))
    control = get_input_bed("input/{}/control".format(wildcards.id))
    bed = signal + control
    bam = list(map(lambda fn: re.sub(r'.bed$', '.bam', fn), bed))
    bai = list(map(lambda fn: re.sub(r'.bed$', '.bai', fn), bed))
    return(bam + bai)

rule peakachu:
    input:
        make_peakachu_input()

rule peakachu_convert_to_bam:
    shell:
        'echo converting bed to bam'

rule peakachu_impl:
    input:
        # make_peakachu_bam_conversion,
        make_peakachu_bam_conversion,
        'input/{id}'
    output:
        'output/peakachu/{id}_peakachu.bed'
    log:
        'log/peakachu/{id}_peakachu.log'
    params:
        genome = config['genome']
    conda:
        'envs/peakachu.yaml'
    script:
        'scripts/peakachu.py'

rule bed_to_bam:
    input:
        bed = '{id}.bed',
        limits = lambda wildcards: "{}.limits".format(config["genome"])
    output:
        '{id}.bam'
    params:
        genome = config['genome']
    shell:
        'bedtools bedtobam -i {input.bed} -g {input.limits}'
### example bed to bam peakachu

# source /home/maticzkd/opt/miniconda3/bin/activate peakachu_0.1.0
# SLOP=10; GENOME=hg19;
# parallel --eta -j 12 "bedtools slop -i {} -b 10 -g ~/genomes/$GENOME.genome | bedtools bedtobam -i - -g ~/genomes/$GENOME.genome | samtools sort > {= s/.bed.gz// =}_slop$SLOP.bam; samtools index {= s/.bed.gz// =}_slop$SLOP.bam" ::: *.bed.gz
