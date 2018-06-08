configfile: 'config.yaml'

include: 'rules/get_input_globs.smk'
include: 'rules/helpers.smk'
include: 'rules/converters.smk'
include: 'rules/targetdist.smk'
include: 'rules/peakachu.smk'

rule all:
    message: 'Help Text Here'
