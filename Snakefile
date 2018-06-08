configfile: 'config.yaml'

include: 'rules/get_input_globs.smk'
include: 'rules/helpers.smk'
include: 'rules/converters.smk'
include: 'rules/targetdist.smk'
include: 'rules/peakachu.smk'

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
    bai = list(map(lambda fn: re.sub(r'.bam$', '.bam.bai', fn), bam))
    return(bam + bai)

def make_peakachu_size_factors(wildcards):
    id = wildcards.id
    signal = "1 " * len(glob.glob('input/{}/signal/*.bam'.format(id)))
    control = "0 " * len(glob.glob('input/{}/control/*.bam'.format(id)))
    return(signal + control)

rule peakachu:
    input:
        make_peakachu_input()

rule peakachu_impl:
    input:
        bam = make_peakachu_bam_conversion,
        dir = 'input/{id}'
    output:
        bed = 'output/peakachu/{id}_peakachu.bed'
    log:
        'log/peakachu/{id}_peakachu.log'
    params:
        genome = config['genome'],
        size_factors = make_peakachu_size_factors

    shadow: "shallow"
    threads: 1
    conda:
        'envs/peakachu.yaml'
    shell:
        'peakachu adaptive '
        '--exp_libs {input.dir}/signal/*.bam '
        '--ctr_libs {input.dir}/control/*.bam '
        '--pairwise_replicates '
        '--max_proc {threads} '
        '-m 0 -n manual --size_factors {params.size_factors} '
        '--output_folder {input.dir} 2>&1 > {log}; '
        'cat {input.dir}/peak_annotations/*.gff | '
        'bedtools sort -i - | '
        'gff2bed > {output.bed}'
