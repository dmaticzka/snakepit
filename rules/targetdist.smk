localrules: get_targetdist

def make_targetdist_input():
    fns = get_input_bed('input/')
    chdir = list(map(lambda fn: fn.replace('input/','output/targetdist/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'.bed$', '.csv', fn), chdir))
    return(chsuff)

rule targetdist:
    input:
        make_targetdist_input()


rule get_targetdist:
    output:
        # candidate for temporary()
        'scripts/targetdist-1.0'
    shell:
        'wget https://github.com/dmaticzka/targetdist/archive/v1.0.tar.gz -O scripts/targetdist-1.0.tar.gz; '
        'tar xvf scripts/targetdist-1.0.tar.gz --directory=scripts/; '
        'wget http://www.bioinf.uni-freiburg.de/~maticzkd/targetdist_annotation_v0_1.tar.bz2 -O scripts/targetdist_annotation_v0_1.tar.bz2; '
        'tar xvf scripts/targetdist_annotation_v0_1.tar.bz2 --directory=scripts/targetdist-1.0/ '


rule targetdist_impl:
    input:
        annotation = 'scripts/targetdist-1.0',
        bed = 'input/{id}.bed'
    output:
        'output/targetdist/{id}.csv'
    log:
        'log/targetdist/{id}.log'
    params:
        genome = config['genome'],
        outprefix = 'output/targetdist/{id}',
    conda:
        '../envs/targetdist.yaml'
    resources:
        vmem = lambda wildcards, attempt: 2**attempt
    shell:
        'scripts/targetdist-1.0/targetdist_{params.genome}.sh {input.bed} {params.outprefix} 2> {log}'
