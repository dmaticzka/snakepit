def get_genome_limits(wildcards):
    return('output/{}.limits'.format(config['genome']))

def strip_bedngz(fns):
    rmgz = list(map(lambda fn: re.sub(r'$', '.gz', fn), fns))
    rmbed = list(map(lambda fn: re.sub(r'$', '.bed', fn), rmgz))
    return(rmbed)

rule genome_impl:
    output:
        limits = '{genome}.limits'
    conda:
        '../envs/mysql.yaml'
    shell:
        'mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e '
        '"select chrom, size from hg19.chromInfo"  > {output.limits}'

rule genome_nohead_impl:
    input:
        limits = '{genome}.limits'
    output:
        limits_nohead = '{genome}.limits_nohead'
    shell:
        'grep -v "^chrom" {input.limits} > {output.limits_nohead}'
