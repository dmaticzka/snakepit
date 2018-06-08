import glob
import os

def get_input_all():
    return(glob.glob('input/*'))

def get_input_dirs():
    return(filter(os.path.isdir, glob.glob('input/*')))

def get_input_bed(base = '.'):
    return(glob.glob('{}/*.bed'.format(base)))
