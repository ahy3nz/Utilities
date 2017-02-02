import os
from optparse import OptionParser
"""
CD to each directory in the input file and submit sbatch.sbatch script
"""

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "somebilayer", dest = "filename")
parser.add_option("--ST", action="store_true", dest ="STrun")
parser.add_option("--MD", action="store_true", dest ="MDrun")
(options, args) = parser.parse_args()

inputfile = open(options.filename,'r')
filenames = inputfile.readlines()
for i, val in enumerate(filenames):
    if "#" not in val and val.rstrip():
        val = val.strip()
        print('*********************************')
        print('{:^20s}'.format(val))
        print('*********************************')
	os.chdir('/global/cscratch1/sd/ahy3nz/Trajectories/{}'.format(val))
	if options.STrun:
		os.system('sbatch {}STsbatch.sbatch'.format(val))
	elif options.MDrun:
		os.system('sbatch {}MDsbatch.sbatch'.format(val))
	else:
		pass



