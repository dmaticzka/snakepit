import glob
import os

rule all:
    message: 'Help Text Here'

def get_input_all():
    return(glob.glob('input/*'))

def get_input_dirs():
    return(filter(os.path.isdir, glob.glob('input/*')))

def get_input_bed():
    return(glob.glob('input/*.bed'))

def make_peakachu_input():
    # this expects directories under input/ that contain bed files
    # in directories named signal and control
    fns = get_input_dirs()
    chdir = list(map(lambda fn: fn.replace('input/','peakachu/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'$', '_peakachu.bed', fn), chdir))
    return(chsuff)

def make_targetdist_input():
    fns = get_input_bed()
    chdir = list(map(lambda fn: fn.replace('input/','targetdist/'), fns))
    chsuff = list(map(lambda fn: re.sub(r'.bed$', '.csv', fn), chdir))
    return(chsuff)

rule peakachu:
    input:
        make_peakachu_input()

rule targetdist:
    input:
        make_targetdist_input()

rule targetdist_impl:
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
