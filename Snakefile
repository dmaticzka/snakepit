# rule all:
    # input: 'targetdist/targetdist_hg19.csv'

rule targetdist:
    input:
        bed = 'input/{id}_{genome}.bed'
    params:
        genome = '{genome}',
        outprefix = 'targetdist/{id}_{genome}'
    output:
        stats = 'targetdist/{id}_{genome}.csv'
    conda:
        'envs/targetdist.yaml'
    shell:
        '~/co/targetdist/targetdist_{params.genome}.sh {input.bed} {params.outprefix}'
