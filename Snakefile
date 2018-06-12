configfile: 'config.yaml'

# generic helper functions
include: 'rules/get_input_globs.smk'
include: 'rules/helpers.smk'

# tools
include: 'rules/bed_correlation_heatmap.smk'
include: 'rules/converters.smk'
include: 'rules/targetdist.smk'
include: 'rules/peakachu.smk'

rule all:
    message: 'Help Text Here'
