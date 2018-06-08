def make_targetdist_input():
    fns = get_input_bed()
    chdir = list(map(lambda fn: fn.replace('input/','output/targetdist/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'.bed$', '.csv', fn), chdir))
    return(chsuff)

rule targetdist:
    input:
        make_targetdist_input()

rule targetdist_impl:
    input:
        'input/{id}.bed'
    output:
        'output/targetdist/{id}.csv'
    log:
        'log/targetdist/{id}.log'
    params:
        genome = config['genome'],
        outprefix = 'output/targetdist/{id}'
    conda:
        '../envs/targetdist.yaml'
    shell:
        '~/co/targetdist/targetdist_{params.genome}.sh {input} {params.outprefix} 2> {log}'
