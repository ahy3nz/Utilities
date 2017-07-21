import gsd
import gsd.hoomd
from optparse import OptionParser

""" HOOMD trjconv

Parameters
----------
input_file : str
    input gsd trajectory
output_file : str
    output gsd trajectory
begin : int
    beginning frame
end : int
    end frame
step : int
    step size for frames

    """

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "trajectory.gsd", dest = "input_file")
parser.add_option("-o", action="store", type="string", default = "output_trajectory.gsd", dest = "output_file")
parser.add_option("-b", action="store", type="int", default = "0", dest = "begin")
parser.add_option("-e", action="store", type="int", default = "-1", dest = "end")
parser.add_option("--dt", action="store", type="int", default = "5", dest = "step")
(options, args) = parser.parse_args()
# Load
#f = gsd.fl.open(name='trajectory.gsd',mode='rb')
#f = gsd.fl.GSDFile(name='trajectory.gsd',mode='rb')
input_file = gsd.hoomd.open(name = options.input_file, mode='rb')
output_file = gsd.hoomd.open(name = options.output_file, mode='wb')
for snapshot_frame in input_file[options.begin : options.end : options.step ]:
    output_file.append(snapshot_frame)
