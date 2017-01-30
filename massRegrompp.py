import numpy as np
import os
import sys
from optparse import OptionParser
"""
    Given a inputfile of all filenames which have all completed ST runs
    Regrompp the simulated tempering output for md simulation
    """

parser = OptionParser()
parser.add_option("-f", action = "store", type = "string", dest = "filename")
(options, args) = parser.parse_args()

inputfile = open(options.filename,'r')
filenames = inputfile.readlines()
failed_list = []
for i, val in enumerate(filenames):
    if "#" not in val and val.rstrip():
        val = val.strip()
        print('*********************************')
        print('{:^20s}'.format(val))
        print('*********************************')
        os.chdir('/raid6/homes/ahy3nz/Trajectories/{}'.format(val))
        if os.path.isfile('/raid6/homes/ahy3nz/Trajectories/{}/ST_{}.gro'.format(val, val)):
            os.system("gmx grompp -f md_{}.mdp -c ST_{}.gro -p {}.top -o md_{} > md_grompp_{}.log 2>&1".format(val, val, val, val, val))
        else:
            failed_list.append(val)
    else:
        pass
print('-------Failed List-------')
for i, element in enumerate(failed_list):
    print('{:^20s}'.format(element))
print('-------------------------')
