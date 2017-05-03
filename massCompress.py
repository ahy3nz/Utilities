import os
import subprocess
from optparse import OptionParser
"""
Compress into tar ball
"""

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "somebilayer", dest = "filename")

(options, args) = parser.parse_args()

os.system('mkdir -p {}'.format(options.filename[:-4]))
with open(options.filename) as f:
    filenames = list(f)
for line in filenames:
    line=line.rstrip()
    print(line)
    command = 'mv {} {}'.format(line, options.filename[:-4])
    print(command)
    p = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
tar_ball = options.filename[:-4]+'.tar'
os.system('tar -cvf {} {}'.format(tar_ball, options.filename[:-4]))
print("Check if tarball was successful, then remove the folder: {}".format(options.filename[:-4]))

