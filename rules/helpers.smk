rule genome_impl:
    output:
        limits = '{genome}.limits'
    # conda:
        # '../envs/mysql.yaml'
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
