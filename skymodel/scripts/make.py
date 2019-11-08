import os

# clean outout directory
outdir = '../output'
os.system('rm {}/*.xml'.format(outdir))
os.system('rm {}/*.fits'.format(outdir))
os.system('rm {}/*.txt'.format(outdir))

# run all preliminary scripts to generate model components

# templates
os.chdir('../known-sources/templates')
# os.system('python make_ic443_model.py')
# os.system('python make_halo_model.py')

# assemble final model
os.chdir('../../scripts')
os.system('python make_model.py')

# make more diagnostic plots
os.chdir('../scripts')
os.system('python show_cutoffs.py')
