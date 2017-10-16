import pdb
import itertools
import os
import msibi_utils
import msibi_utils.animate_rdf
import subprocess


# Try to animate all the rdfs
atom_types = ['Na' ,  
           'C2', 'C1' , 
              'Qa' ,  'Q0' ]
exclusions =[ ('C1', 'C1')]

for atomtype_i, atomtype_j in itertools.combinations_with_replacement(atom_types,2):
    if (atomtype_i, atomtype_j) not in exclusions:
        p = subprocess.Popen("cp rdfs/{0}-{1}-state_A.txt rdfs/{0}-{1}-state-7.500.txt".format(
            atomtype_i, atomtype_j), shell=True, stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT)
        p.wait()

        msibi_utils.animate_rdf.animate_pair_at_state(atomtype_i, atomtype_j, "state-7.500", 159, "./rdfs",
            use_agg=True, to_angstrom =1, to_kcalpermol=1, n_skip=20)
