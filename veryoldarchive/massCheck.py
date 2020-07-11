import numpy as np
import os
import sys
from optparse import OptionParser

"""
Given a list of filenames, check which simulations have finished
based on the respective grofile printed
"""
parser = OptionParser()
parser.add_option("-f", action = "store", type = "string", dest = "filename")
parser.add_option("--ST", action="store_true", dest = "STrun")
parser.add_option("--MD", action="store_true", dest = "MDrun")
(options, args) = parser.parse_args()

if not options.STrun and not options.MDrun:
    sys.exit("Specify ST or MD file extension")


inputfile = open(options.filename,'r')
filenames = inputfile.readlines()
success_list = []
fail_list = []

for i, val in enumerate(filenames):
    if "#" not in val and val.rstrip():
        val = val.strip()
        print('*********************************')
        print('{:^20s}'.format(val))
        print('*********************************')
        os.chdir('/raid6/homes/ahy3nz/Trajectories/{}'.format(val))
        if options.STrun:
            if os.path.isfile('/raid6/homes/ahy3nz/Trajectories/{}/ST_{}.gro'.format(val, val)):
                success_list.append(val)
            else:
                fail_list.append(val)
        elif options.MDrun:
            if os.path.isfile('/raid6/homes/ahy3nz/Trajectories/{}/md_{}.gro'.format(val, val)):
                success_list.append(val)
            else:
                fail_list.append(val)
        else:
            pass
inputfile.close()
os.chdir('/raid6/homes/ahy3nz/Trajectories')
successfile = open(str(options.filename[:-4] + '_done.dat'),'w')
for i, val in enumerate(success_list):
    successfile.write(val+"\n")
successfile.close()

failfile = open(str(options.filename[:-4] + '_ongoing.dat'),'w')
for i, val in enumerate(fail_list):
    failfile.write(val+"\n")
failfile.close()
