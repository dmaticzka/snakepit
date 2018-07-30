
def make_bed_correlation_heatmap_input():
    # this expects all input bed files within a directory under folder input
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/bed_correlation_heatmap/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_correlation_heatmap.pdf', fn), chdir))
    chsuff_nonumbers = list(map(lambda fn: re.sub(r'$', '_correlation_heatmap_nonumbers.pdf', fn), chdir))
    return(chsuff + chsuff_nonumbers)

rule bed_correlation_heatmap:
    input:
        make_bed_correlation_heatmap_input()


rule bed_correlation_heatmap_npz:
    input:
        dir = 'input/{dir}',
        limits = get_genome_limits
    output:
        npz = 'output/bed_correlation_heatmap/{dir}_correlation_heatmap.npz'
    params:
        bed = collect_input_bedngz,
        window_size = config['bed_correlation_heatmap_window_size']
    conda:
        '../envs/bed_correlation_heatmap.yaml'
    shell:
        'LABELS=`basename -s .gz {params.bed} | sed "s:.bed$::g"`;'
        'bedtools makewindows -g {input.limits} -w {params.window_size} | '
        'multiBedSummary.py '
        '--regions - '
        '--outFileName {output.npz} '
        '--bedfiles {params.bed} '
        '--labels $LABELS; '


rule bed_correlation_heatmap_impl:
    input:
        npz = 'output/bed_correlation_heatmap/{dir}_correlation_heatmap.npz',
    output:
        pdf = 'output/bed_correlation_heatmap/{dir}_correlation_heatmap.pdf',
    conda:
        '../envs/bed_correlation_heatmap.yaml'
    shell:
        'plotCorrelation '
        '-in {input.npz} '
        '-o  {output.pdf} '
        '--outFileCorMatrix {output.csv} '
        '--corMethod spearman '
        '--whatToPlot heatmap '
        '--skipZeros '
        '--plotNumbers'


rule bed_correlation_heatmap_nonumbers_impl:
    input:
        npz = 'output/bed_correlation_heatmap/{dir}_correlation_heatmap.npz',
    output:
        pdf = 'output/bed_correlation_heatmap/{dir}_correlation_heatmap_nonumbers.pdf',
    conda:
        '../envs/bed_correlation_heatmap.yaml'
    shell:
        'plotCorrelation '
        '-in {input.npz} '
        '-o  {output.pdf} '
        '--corMethod spearman '
        '--whatToPlot heatmap '
        '--skipZeros '
