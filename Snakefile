configfile: 'config.yaml'

include:
    'rules/get_input_globs.smk',
    'rules/helpers.smk'

rule all:
    message: 'Help Text Here'

def make_peakachu_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/peakachu/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_peakachu.bed', fn), chdir))
    return(chsuff)

def make_targetdist_input():
    fns = get_input_bed()
    chdir = list(map(lambda fn: fn.replace('input/','output/targetdist/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'.bed$', '.csv', fn), chdir))
    return(chsuff)

rule peakachu:
    input:
        make_peakachu_input()

rule peakachu_impl:
    input:
        'input/{id}'
    output:
        'output/peakachu/{id}_peakachu.bed'
    log:
        'log/peakachu/{id}_peakachu.log'
    params:
        genome = config['genome']
    conda:
        'envs/peakachu.yaml'
    shell:
        'echo {input} {output}'

rule targetdist:
    input:
        make_targetdist_input()

rule targetdist_impl:
    input:
        'input/{id}.bed'
    output:
        'output/targetdist/{id}.csv'
    log:
        'log/targetdist/{id}.log'
    params:
        genome = config['genome'],
        outprefix = 'output/targetdist/{id}'
    conda:
        'envs/targetdist.yaml'
    shell:
        '~/co/targetdist/targetdist_{params.genome}.sh {input} {params.outprefix} 2> {log}'

### example bed to bam peakachu

# source /home/maticzkd/opt/miniconda3/bin/activate peakachu_0.1.0
# SLOP=10; GENOME=hg19;
# parallel --eta -j 12 "bedtools slop -i {} -b 10 -g ~/genomes/$GENOME.genome | bedtools bedtobam -i - -g ~/genomes/$GENOME.genome | samtools sort > {= s/.bed.gz// =}_slop$SLOP.bam; samtools index {= s/.bed.gz// =}_slop$SLOP.bam" ::: *.bed.gz
