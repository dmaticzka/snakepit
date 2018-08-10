def make_bed_pca_input():
    # this expects all input bed files within a directory under folder input
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/bed_pca/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_pca.pdf', fn), chdir))
    return(chsuff)


rule bed_correlation_heatmap:
    input:
        make_bed_pca_input()


rule bed_pca_npz:
    input:
        dir = 'input/{dir}',
        limits = get_genome_limits
    output:
        npz = 'output/bed_pca/{dir}.npz'
    params:
        bed = collect_input_bedngz,
        window_size = config['bed_pca_window_size']
    conda:
        '../envs/bed_pca.yaml'
    shell:
        'LABELS=`basename -s .gz {params.bed} | sed "s:.bed$::g"`;'
        'bedtools makewindows -g {input.limits} -w {params.window_size} | '
        'multiBedSummary.py '
        '--regions - '
        '--outFileName {output.npz} '
        '--bedfiles {params.bed} '
        '--labels $LABELS; '


rule bed_pca_impl:
    input:
        npz = 'output/bed_pca/{dir}_correlation_heatmap.npz',
    output:
        pdf = 'output/bed_pca/{dir}_correlation_heatmap.pdf',
        csv = 'output/bed_pca/{dir}_correlation_heatmap.csv',
    conda:
        '../envs/bed_pca.yaml'
    shell:
        'plotCorrelation '
        '-in {input.npz} '
        '-o  {output.pdf} '
        '--outFileCorMatrix {output.csv} '
        '--corMethod spearman '
        '--whatToPlot heatmap '
        '--skipZeros '
        '--plotNumbers'
