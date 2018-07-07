rule gunzip:
    input:
        '{id}.gz'
    output:
        temporary('{id,^(.gz$)}')
    shell:
        'zcat {input} > {output}'


rule bed_to_bam:
    input:
        bed = '{id}.bed',
        limits = lambda wildcards: "{}.limits".format(config["genome"])
    output:
        temporary('{id}.bam')
    params:
        genome = config['genome']
    conda:
        '../envs/bedtobam.yaml'
    shell:
        'bedtools bedtobam -i {input.bed} -g {input.limits} | '
        'samtools sort > {output}'


rule index_bam:
    input:
        '{id}.bam'
    output:
        '{id}.bam.bai'
    conda:
         '../envs/samtools.yaml'
    shell:
        'samtools index {input}'


# only keep bam alignments that have their id in an accompanying bed
rule filter_bam_by_bed:
    input:
        bam = 'input/{id}.bam',
        bed = 'input/{id}.bed',
    output:
        bam = 'output/filter_bam_by_bed/{id}.bam',
    resources:
        vmem = lambda wildcards, attempt: int(2*(attempt+3))
    conda:
        '../envs/samtools.yaml'
    shell:
        '( '
        'samtools view -H {input.bam}; '
        'samtools view {input.bam} | '
        'fgrep -w -f <(cut -f 4 {input.bed}) '
        ') | '
        'samtools view -b - | '
        'samtools sort > {output.bam}; '
