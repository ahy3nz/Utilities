import pdb
import itertools
import os
import msibi_utils
import msibi_utils.parse_logfile
import msibi_utils.plot_fit
import subprocess

#logfile_info = msibi_utils.parse_logfile.parse_logfile('optimization_9-4.log')
msibi_utils.plot_fit.plot_all_fits('optimization_9-4.log', ylims=(-1, 1))

