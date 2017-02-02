import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
"""
    Given a inputfile of all filenames
    Remove PBCs, analyze it, and print out the temperature over time
    Assumes certain file organization and naming 
    """

parser = OptionParser()
parser.add_option("-f", action = "store", type = "string", dest = "filename")
(options, args) = parser.parse_args()

inputfile = open(options.filename,'r')
filenames = inputfile.readlines()
for i, val in enumerate(filenames):
    if "#" not in val and val.rstrip():
        val = val.strip()
        print('*********************************')
        print('{:^20s}'.format(val))
        print('*********************************')
        os.chdir('/Users/ahy3nz/Trajectories/{}'.format(val))
        os.system("echo '0' | gmx trjconv -f md_{}.xtc -s md_{}.tpr -pbc mol -dt 20 -o nopbc".format(val,val))
        os.system("echo '0' | gmx trjconv -f md_{}.xtc -s md_{}.tpr -pbc mol -dt 20 -b 80000 -e 100000 \
                -o last20".format(val,val))
        os.system('cp ~/Programs/Analysis/Bilayers/AnalyzeBilayer.py .')
        os.system('python AnalyzeBilayer.py -f nopbc.xtc -c md_{}.gro -o full'.format(val))
        os.system('python AnalyzeBilayer.py -f last20.xtc -c md_{}.gro -o last20'.format(val))
        os.system('cp ~/Programs/Utilities/Plotting/PlotXVG.py .')
        #os.system("echo '12' | gmx energy -f ST_{}.edr -o temp".format(val))
        #os.system('python PlotXVG.py -f temp.xvg')
    else:
        pass

