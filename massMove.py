import os
from optparse import OptionParser
"""
Move things to rahman, from rahman, or to filey based on folder name in the input file
"""

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "somebilayer", dest = "filename")
parser.add_option("--toRAHMAN", action="store_true", dest = "toRahman")
parser.add_option("--toFILEY", action="store_true", dest="toFiley")
parser.add_option("--fromRAHMAN", action="store_true", dest = "fromRahman")
parser.add_option("--toACCRE", action="store_true", dest = "toACCRE")

(options, args) = parser.parse_args()

inputfile = open(options.filename,'r')
filenames = inputfile.readlines()
for i, val in enumerate(filenames):
    if "#" not in val and val.rstrip():
        val = val.strip()
        print('*********************************')
        print('{:^20s}'.format(val))
        print('*********************************')
        if options.toRahman:
	    	os.system('scp -r {}/{} $RAHMAN:Trajectories/'.format(os.getcwd(), val))
            os.system('scp -r {}/{} $RAHMAN:Trajectories/'.format(os.getcwd(), inputfile))
	    elif options.toFiley:
	    	os.system('scp -r {}/{} $FILEY:Trajectories/'.format(os.getcwd,val))
            os.system('scp -r {}/{} $FILEY:Trajectories/'.format(os.getcwd, inputfile))
        elif options.toACCRE:
            os.system('scp -r {}/{} $ACCRE:Trajectories'.format(os.getcwd, val))
            os.system('scp -r {}/{} $ACCRE:Trajectories'.format(os.getcwd, inputfile))
	    elif options.fromRahman:
	    	os.system('scp -r $RAHMAN:Trajectories/{} {}'.format(val.os.getcwd))
	    	#os.system('scp -r $RAHMAN:Programs/setup/Bilayer/{} /global/cscratch1/sd/ahy3nz/Trajectories'.format(val))
	    else:
	    	pass



