import os
import subprocess
from optparse import OptionParser
"""
Move things to rahman, from rahman, or to filey based on folder name in the input file
Compress into tar bal and move the tarball
"""

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "somebilayer", dest = "filename")
parser.add_option("--toRAHMAN", action="store_true", dest = "toRahman")
parser.add_option("--toFILEY", action="store_true", dest="toFiley")
parser.add_option("--fromRAHMAN", action="store_true", dest = "fromRahman")
parser.add_option("--toACCRE", action="store_true", dest = "toACCRE")
parser.add_option("--toTITAN", action="store_true", dest = "toTITAN")
parser.add_option("--toEDISON", action="store_true", dest = "toEDISON")

(options, args) = parser.parse_args()

with open(options.filename) as f:
    filenames = list(f)
#options.filename = '90_WVTR_4-29b'
#filenames = ['90_WVTR_4-29b']
# Compress everything before sending, folders and filename
to_compress = ' '.join([line.rstrip() for line in filenames])
tar_ball = options.filename[:-4]+'.tar'
os.system('tar cvf {} {} {}'.format(tar_ball, to_compress, options.filename))
print(to_compress)
if options.toRahman:
    os.system('cp -r {}/{} ~/Trajectories/'.format(os.getcwd(), tar_ball))
elif options.toFiley:
    os.system('scp -r {}/{} $FILEY:Trajectories/'.format(os.getcwd(), tar_ball))
elif options.toACCRE:
    os.system('scp -r {}/{} $ACCRE:Trajectories'.format(os.getcwd(), tar_ball))
elif options.toTITAN:
    os.system('scp -r {}/{} $TITAN:/lustre/atlas/scratch/ahy3nz/mat149/Trajectories'.format(os.getcwd(), tar_ball))
elif options.toEDISON:
    p=subprocess.Popen('scp -r {}/{} $EDISON:/scratch2/scratchdirs/ahy3nz/Trajectories'.format(os.getcwd(), tar_ball), shell= True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()
    print(p.stdout.read())
elif options.fromRahman:
    os.system('scp -r $RAHMAN:Trajectories/{} {}'.format(val,tar_ball))
	#os.system('scp -r $RAHMAN:Programs/setup/Bilayer/{} /global/cscratch1/sd/ahy3nz/Trajectories'.format(val))
else:
    pass



