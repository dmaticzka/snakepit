import glob

rule all:
    message: 'Help Text Here'

def get_input_bed():
    return(glob.glob('input/*.bed'))

def make_targetdist_input():
    fns = get_input_bed()
    chdir = list(map(lambda fn: fn.replace('input/','targetdist/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'.bed$', '.csv', fn), chdir))
    return(chsuff)

rule targetdist_all:
    input:
        make_targetdist_input()

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
