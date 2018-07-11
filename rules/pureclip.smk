localrules: pureclip


def make_pureclip_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/pureclip/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_pureclip_sites.bed', fn), chdir))
    return(chsuff)


def make_pureclip_onlysignal_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/pureclip_onlysignal/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_pureclip_sites.bed', fn), chdir))
    return(chsuff)


def make_pureclip_bam_fmt_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/pureclip_bam_fmt/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_pureclip_sites.bed', fn), chdir))
    return(chsuff)


def make_pureclip_onlysignal_bam_fmt_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','output/pureclip_onlysignal_bam_fmt/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_pureclip_sites.bed', fn), chdir))
    return(chsuff)


rule pureclip:
    input:
        make_pureclip_input()


rule pureclip_onlysignal:
    input:
        make_pureclip_onlysignal_input()


rule pureclip_bam_fmt:
    input:
        make_pureclip_bam_fmt_input()


rule pureclip_onlysignal_bam_fmt:
    input:
        make_pureclip_onlysignal_bam_fmt_input()


rule pureclip_impl:
    input:
        sig_bam = 'output/pureclip/{id}/signal.bam',
        sig_bai = 'output/pureclip/{id}/signal.bam.bai',
        ctl_bam = 'output/pureclip/{id}/control.bam',
        ctl_bai = 'output/pureclip/{id}/control.bam.bai',
    output:
        sites = 'output/pureclip/{id}_pureclip_sites.bed',
        regions = 'output/pureclip/{id}_pureclip_regions.bed',
    log:
        'log/pureclip/{id}_pureclip.log',
    params:
        genome = '~/genomes/{}.fa'.format(config['genome']),
    conda:
        '../envs/pureclip.yaml',
    # use 20, 40, 80, 160 GB of RAM
    resources:
        vmem = lambda wildcards, attempt: math.ceil(10*(2**(attempt)))
    threads: 1
    shell:
        'pureclip '
        '-i {input.sig_bam} -bai {input.sig_bai} '
        '-ibam {input.ctl_bam} -ibai {input.ctl_bai} '
        '-g {params.genome} '
        '-iv "chr1;chr2;chr3;" '
        '-nt {threads} -nta {threads} '
        '-o {output.sites} '
        '-ld '
        '-or {output.regions} 2>&1 > {log}; '


rule pureclip_onlysignal_impl:
    input:
        sig_bam = 'output/pureclip_onlysignal/{id}/signal.bam',
        sig_bai = 'output/pureclip_onlysignal/{id}/signal.bam.bai',
    output:
        sites = 'output/pureclip_onlysignal/{id}_pureclip_sites.bed',
        regions = 'output/pureclip_onlysignal/{id}_pureclip_regions.bed',
    log:
        'log/pureclip_onlysignal_bam_fmt/{id}_pureclip.log',
    params:
        genome = '~/genomes/{}.fa'.format(config['genome']),
    conda:
        '../envs/pureclip.yaml',
    # use 20, 40, 80, 160 GB of RAM
    resources:
        vmem = lambda wildcards, attempt: math.ceil(10*(2**(attempt)))
    threads: 1
    shell:
        'pureclip '
        '-i {input.sig_bam} -bai {input.sig_bai} '
        '-g {params.genome} '
        '-iv "chr1;chr2;chr3;" '
        '-ld '
        '-nt {threads} -nta {threads} '
        '-o {output.sites} '
        '-or {output.regions} 2>&1 > {log}; '


rule pureclip_bam_fmt_impl:
    input:
        sig_bam = 'output/pureclip_bam_fmt/{id}/signal.bam',
        sig_bai = 'output/pureclip_bam_fmt/{id}/signal.bam.bai',
        ctl_bam = 'output/pureclip_bam_fmt/{id}/control.bam',
        ctl_bai = 'output/pureclip_bam_fmt/{id}/control.bam.bai',
    output:
        sites = 'output/pureclip_bam_fmt/{id}_pureclip_sites.bed',
        regions = 'output/pureclip_bam_fmt/{id}_pureclip_regions.bed',
    log:
        'log/pureclip_bam_fmt/{id}_pureclip.log',
    params:
        genome = '~/genomes/{}.fa'.format(config['genome']),
    conda:
        '../envs/pureclip.yaml',
    # use 24, 48, 96, ... GB of RAM
    resources:
        vmem = lambda wildcards, attempt: math.ceil((2**(attempt)))
    threads: 12
    shell:
        'pureclip '
        '-i {input.sig_bam} -bai {input.sig_bai} '
        '-ibam {input.ctl_bam} -ibai {input.ctl_bai} '
        '-g {params.genome} '
        '-iv "chr15" '
        '-bw 25 '
        '-nt {threads} -nta {threads} '
        '-o {output.sites} '
        '-or {output.regions} 2>&1 > {log}; '


rule pureclip_onlysignal_bam_fmt_impl:
    input:
        sig_bam = 'output/pureclip_onlysignal_bam_fmt/{id}/signal.bam',
        sig_bai = 'output/pureclip_onlysignal_bam_fmt/{id}/signal.bam.bai',
    output:
        sites = 'output/pureclip_onlysignal_bam_fmt/{id}_pureclip_sites.bed',
        regions = 'output/pureclip_onlysignal_bam_fmt/{id}_pureclip_regions.bed',
    log:
        'log/pureclip_onlysignal_bam_fmt/{id}_pureclip.log',
    params:
        genome = '~/genomes/{}.fa'.format(config['genome']),
    conda:
        '../envs/pureclip.yaml',
    # use 20, 40, 80, 160 GB of RAM
    resources:
        vmem = lambda wildcards, attempt: math.ceil(10*(2**(attempt)))
    threads: 1
    shell:
        'pureclip '
        '-i {input.sig_bam} -bai {input.sig_bai} '
        '-g {params.genome} '
        '-iv "chr1;chr2;chr3;" '
        '-ld '
        '-nt {threads} -nta {threads} '
        '-o {output.sites} '
        '-or {output.regions} 2>&1 > {log}; '


rule pureclip_combine_bed_to_bam:
    input:
        dir = 'input/{id}/{sigtype}',
        limits = lambda wildcards: "{}.limits".format(config["genome"]),
    output:
        combined_bam = 'output/pureclip/{id}/{sigtype}.bam',
    conda:
        '../envs/bedtobam.yaml'
    shell:
        'cat {input.dir}/*.bed | '
        'sort -k1,1 -k2,2n | '
        'bedtools slop -b 10 -g {input.limits} -i - | '
        'bedtools bedtobam -i - -g {input.limits} | '
        'samtools sort > {output.combined_bam}; '


rule pureclip_onlysignal_combine_bed_to_bam:
    input:
        dir = 'input/{id}/signal',
        limits = lambda wildcards: "{}.limits".format(config["genome"]),
    output:
        combined_bam = 'output/pureclip_onlysignal/{id}/signal.bam',
    conda:
        '../envs/bedtobam.yaml'
    shell:
        'cat {input.dir}/*.bed | '
        'bedtools slop -b 10 -g {input.limits} -i - | '
        'bedtools bedtobam -i - -g {input.limits} | '
        'samtools sort > {output.combined_bam}; '


rule pureclip_onlysignal_combine_bam_filter_fmt:
    input:
        dir = 'input/{id}/signal',
        limits = lambda wildcards: "{}.limits".format(config["genome"]),
    output:
        combined_bam = 'output/pureclip_onlysignal_bam_fmt/{id}/signal.bam',
    params:
        merged_bam = 'output/pureclip_onlysignal_bam_fmt/{id}/merged.bam',
    conda:
        '../envs/bedtobam.yaml'
    shell:
        'samtools merge -f -u {params.merged_bam} {input.dir}/*.bam; '
        'samtools index {params.merged_bam}; '
        'samtools view -hb -f 66 -o {output.combined_bam} {params.merged_bam}; '


rule pureclip_combine_bam_filter_fmt:
    input:
        dir = 'input/{id}/{sourcedir}',
        limits = lambda wildcards: "{}.limits".format(config["genome"]),
    output:
        combined_bam = 'output/pureclip_bam_fmt/{id}/{sourcedir}.bam',
    params:
        merged_bam = 'output/pureclip_bam_fmt/{id}/{sourcedir}_merged.bam',
    conda:
        '../envs/bedtobam.yaml'
    shell:
        'samtools merge -f -u {params.merged_bam} {input.dir}/*.bam; '
        'samtools index {params.merged_bam}; '
        'samtools view -hb -f 66 -o {output.combined_bam} {params.merged_bam}; '
