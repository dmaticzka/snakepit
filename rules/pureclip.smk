localrules: pureclip, combine_bed_to_bam


def make_pureclip_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/pureclip/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_pureclip_sites.bed', fn), chdir))
    return(chsuff)


rule pureclip:
    input:
        make_pureclip_input()


rule pureclip_impl:
    input:
        sig_bam = 'output/pureclip/{id}/signal.bam',
        sig_bai = 'output/pureclip/{id}/signal.bam.bai',
        ctl_bam = 'output/pureclip/{id}/control.bam',
        ctl_bai = 'output/pureclip/{id}/control.bam.bai',
    output:
        sites = 'output/pureclip/{id}_pureclip_sites.bed',
        regions = 'output/pureclip/{id}_pureclip_regions.bed',
    params:
        genome = '~/genomes/{}.fa'.format(config['genome']),
    conda:
        '../envs/pureclip.yaml'
    threads: 8
    shell:
        'pureclip '
        '-i {input.sig_bam} -bai {input.sig_bai} '
        '-ibam {input.ctl_bam} -ibai {input.ctl_bai} '
        '-g {params.genome} '
        '-iv "chr1;chr2;chr3;" '
        '-nt {threads} '
        '-o {output.sites} '
        '-or {output.regions}; '


rule combine_bed_to_bam:
    input:
        dir = 'input/{id}/{sigtype}',
        limits = lambda wildcards: "{}.limits".format(config["genome"]),
    output:
        combined_bam = 'output/pureclip/{id}/{sigtype}.bam',
    conda:
        '../envs/bedtobam.yaml'
    shell:
        'cat {input.dir}/*.bed | '
        'bedtools bedtobam -i - -g {input.limits} | '
        'samtools sort > {output.combined_bam}; '
