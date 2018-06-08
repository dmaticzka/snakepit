rule bed_to_bam:
    input:
        bed = '{id}.bed',
        limits = lambda wildcards: "{}.limits".format(config["genome"])
    output:
        '{id}.bam'
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
