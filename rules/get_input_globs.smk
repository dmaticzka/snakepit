import glob
import os

## get functions work independently of the snakemake wildcards-style arguments

def get_input_all():
    return(glob.glob('input/*'))

def get_input_dirs():
    return(filter(os.path.isdir, glob.glob('input/*')))

def get_input_bed(base = '.'):
    return(glob.glob('{}/*.bed'.format(base)))

def get_input_bedngz(base = '.'):
    bed = glob.glob('{}/*.bed'.format(base))
    bedgz = glob.glob('{}/*.bed.gz'.format(base))
    return(bed + bedgz)

## collect functions use the snakemake wildcards-style arguments

def collect_input_bedngz(wildcards):
    bed = glob.glob('input/{}/*.bed'.format(wildcards.dir))
    bedgz = glob.glob('input/{}/*.bed.gz'.format(wildcards.dir))
    return(bed + bedgz)
