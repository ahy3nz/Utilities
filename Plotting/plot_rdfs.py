import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages
import itertools

""" Plot RDFs for a couple MSIBI iterations"""
final_step = 24
first_step = 0

#beadtypes = ['P3', 'Nda', 'Na', 'C1', 'Qa', 'Q0']
beadtypes = ['Na', 'C1', 'Qa', 'Q0']
#outpdf = PdfPages(('MSIBIresults.pdf'))
outpdf = PdfPages(('MSIBIsmall.pdf'))

for i,j in itertools.combinations_with_replacement(beadtypes, 2):
    if 'C1' not in i or 'C1' not in j:
        target_A = np.loadtxt('{}-{}-state_A.txt'.format(i,j))
        first_A = np.loadtxt("pair_{}-{}-state_state-7.500-step{}.txt".format(i,j, first_step))
        final_A = np.loadtxt("pair_{}-{}-state_state-7.500-step{}.txt".format(i,j, final_step))
        
        #target_B = np.loadtxt("{}-{}-state_B.txt".format(i,j))
        #first_B = np.loadtxt("pair_{}-{}-state_state-2.500-step0.txt".format(i,j))
        #final_B = np.loadtxt("pair_{}-{}-state_state-2.500-step49.txt".format(i,j))
        
        fig, axarray = plt.subplots(1,1)
        axarray.set_title("Bulk fluid (900K) {}-{}".format(i,j))
        axarray.plot(target_A[:,0], target_A[:,1], color='black', 
                label="$g^{AA}_{target}$")
        axarray.plot(first_A[:,0], first_A[:,1], color='red', marker='s', 
                label="$g^{{CG}}_{{{}}}$".format(first_step))
        axarray.plot(final_A[:,0], final_A[:,1], color='blue', marker='o', 
                label="$g^{{CG}}_{{{}}}$".format(final_step))
        axarray.legend()
        axarray.set_xlabel("Distance (nm)")
        
        #axarray[1].set_title("Bilayer (305K) {}-{}".format(i,j))
        #axarray[1].plot(target_B[:,0], target_B[:,1], color='black', label="$g^{AA}_{target}$")
        #axarray[1].plot(first_B[:,0], first_B[:,1], color='red', marker='s', label="$g^{CG}_{original}$")
        ##axarray[1].plot(final_B[:,0], final_B[:,1], color='blue', marker='o', label="50th iteration CG")
        #axarray[1].legend()
        #axarray[1].set_xlabel("Distance (nm)")
        
        plt.tight_layout()
        outpdf.savefig(fig)
        plt.savefig('{}-{}.svg'.format(i,j), transparent=True)
outpdf.close()
