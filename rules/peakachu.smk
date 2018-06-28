localrules: bed_slop10, peakachu_bed, peakachu_initial_peaks


def make_peakachu_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/peakachu/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_peakachu.bed', fn), chdir))
    initpeaks = list(map(lambda fn: re.sub(r'$', '_peakachu_initial_peaks.csv.bed', fn), chdir))
    return(chsuff + initpeaks)


def make_peakachu_bed_input(wildcards):
    signal = get_input_bedngz("input/{}/signal".format(wildcards.id))
    control = get_input_bedngz("input/{}/control".format(wildcards.id))
    bed = list(map(lambda fn: re.sub(r'.gz$', '', fn), signal + control))
    return(bed)


def make_peakachu_bam_conversion(wildcards):
    signal = get_input_bedngz("input/{}/signal".format(wildcards.id))
    control = get_input_bedngz("input/{}/control".format(wildcards.id))
    chdir = list(map(lambda fn: fn.replace('input/', 'output/peakachu/'), signal + control))
    bed = list(map(lambda fn: re.sub(r'.gz$', '', fn), chdir))
    bam = list(map(lambda fn: re.sub(r'.bed$', '_slop10.bam', fn), bed))
    bai = list(map(lambda fn: re.sub(r'.bam$', '.bam.bai', fn), bam))
    return(bam + bai)


def make_peakachu_size_factors(wildcards):
    id = wildcards.id
    signal = "1 " * len(glob.glob('input/{}/signal/*.bed'.format(id)))
    control = "0.75 " * len(glob.glob('input/{}/control/*.bed'.format(id)))
    return(signal + control)


rule peakachu:
    input:
        make_peakachu_input()


rule bed_slop10:
    input:
        bed = 'input/{dir}/{lib}/{id}.bed',
        limits = lambda wildcards: "{}.limits".format(config["genome"])
    output:
        'output/peakachu/{dir}/{lib}/{id}_slop10.bed'
    conda:
        '../envs/bedtools.yaml'
    shell:
        'bedtools slop -i {input.bed} -g {input.limits} -b 10 > {output}'


rule peakachu_impl:
    input:
        bed = temporary(make_peakachu_bed_input),
        bam = temporary(make_peakachu_bam_conversion),
        dir = 'input/{id}'
    output:
        dir = 'output/peakachu/{id}'
    log:
        'log/peakachu/{id}_peakachu.log'
    params:
        genome = config['genome'],
        size_factors = make_peakachu_size_factors
    shadow: "shallow"
    threads: 6
    resources:
        vmem = lambda wildcards, attempt: int((2**(attempt+2))/6)
    conda:
        '../envs/peakachu.yaml'
    shell:
        'peakachu adaptive '
        '--exp_libs {output.dir}/signal/*.bam '
        '--ctr_libs {output.dir}/control/*.bam '
        '--pairwise_replicates '
        '--max_proc {threads} '
        '-m 0 -n manual --size_factors {params.size_factors} '
        '--output_folder {output.dir} 2>&1 > {log}; '


rule peakachu_bed:
    input:
        peakachudir = 'output/peakachu/{id}',
    output:
        peaks_bed = 'output/peakachu/{id}_peakachu.bed',
    shell:
        'GFF=({input.peakachudir}/peak_annotations/*.gff); '
        'if [[ -f ${{GFF[0]}} ]]; '
        'then '
        '  cat {input.peakachudir}/peak_annotations/*.gff | '
        '  bedtools sort -i - | '
        '  gff2bed | '
        '  cut -f 1-6 | '
        '  awk \'BEGIN{{OFS=\"\\t\"}}{{$5=255; print}}\' > {output.peaks_bed}; '
        'else '
        '  touch {output.peaks_bed}; '
        'fi; '


rule peakachu_initial_peaks:
    input:
        peakachudir = 'output/peakachu/{id}',
    output:
        initial_peaks = 'output/peakachu/{id}_peakachu_initial_peaks.csv',
        initial_peaks_bed = 'output/peakachu/{id}_peakachu_initial_peaks.csv.bed',
    conda:
        '../envs/misc_scripts.yaml'
    shell:
        'if [[ -f {input.peakachudir}/initial_peaks.csv ]]; '
        'then '
        '  cp {input.peakachudir}/initial_peaks.csv {output.initial_peaks}; '
        'else '
        '  touch {output.initial_peaks}; '
        'fi;'
        'peakachu_initial_peaks_to_bed.R {output.initial_peaks}; '
