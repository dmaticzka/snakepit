
def make_bed_correlation_heatmap_input():
    # this expects all input bed files within a directory under folder input
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/bed_correlation_heatmap/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_correlation_heatmap.pdf', fn), chdir))
    return(chsuff)

rule bed_correlation_heatmap:
    input:
        make_bed_correlation_heatmap_input()

rule bed_correlation_heatmap_impl:
    input:
        dir = 'input/{dir}',
        limits = get_genome_limits
    output:
        npz = 'output/bed_correlation_heatmap/{dir}_correlation_heatmap.npz',
        pdf = 'output/bed_correlation_heatmap/{dir}_correlation_heatmap.pdf'
    params:
        bed = collect_input_bedngz,
        labels = collect_input_bedngz
    conda:
        '../envs/bed_correlation_heatmap.yaml'
    shell:
        'bedtools makewindows -g {input.limits} -w 200000 | '
        'multiBedSummary.py '
        '--regions - '
        '--outFileName {output.npz} '
        '--bedfiles {params.bed} '
        '--labels {params.bed}'
#
# --labels
#
# plotCorrelation \
# -in counts_wins200.npz \
# -o  counts_wins200.pdf \
# --outFileCorMatrix  counts_wins200.csv \
# --corMethod spearman \
# --whatToPlot heatmap \
# --skipZeros \
# --plotNumbers
