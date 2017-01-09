import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", action = "store", type = "string", dest = "filename",
    default="Stage1_Weak0_pullx.xvg")
(options, args) = parser.parse_args()

filename = options.filename
xvgfile = open(filename, 'r')
xvglines = xvgfile.readlines()
data = list()
legend = []
for i, line in enumerate(xvglines):
    if '@' in line and 'title' in line:
        first_apostrophe = line.find('\"')
        second_apostrophe = line.rfind('\"')
        title = line[first_apostrophe+1: second_apostrophe]
        #Find and replace superscripts and \N
        title = title.replace('\\S', '$^{')
        title = title.replace('\\s', '$_{')
        title = title.replace('\\N', '}$')
    if '@' in line and 'xaxis' in line:
        first_apostrophe = line.find('\"')
        second_apostrophe = line.rfind('\"')
        xlabel = line[first_apostrophe+1: second_apostrophe]
        xlabel = xlabel.replace('\\S', '$^{')
        xlabel = xlabel.replace('\\s', '$_{')
        xlabel = xlabel.replace('\\N', '}$')
    if '@' in line and 'yaxis' in line:
        first_apostrophe = line.find('\"')
        second_apostrophe = line.rfind('\"')
        ylabel = line[first_apostrophe+1: second_apostrophe]
        ylabel = ylabel.replace('\\S', '$^{')
        ylabel = ylabel.replace('\\s', '$_{')
        ylabel = ylabel.replace('\\N', '}$')
    if '@' in line and 'legend' in line and 's' in line:
        first_apostrophe = line.find('\"')
        second_apostrophe = line.rfind('\"')
        legend_entry = line[first_apostrophe+1: second_apostrophe]
        legend_entry = legend_entry.replace('\\S', '$^{')
        legend_entry = legend_entry.replace('\\s', '$_{')
        legend_entry  = legend_entry.replace('\\N', '}$')
        legend.append(legend_entry)
    if '#' not in line and '@' not in line:
        items = line.split()
        data.append((items))

    else:
        pass

time = len(data)
n_columns = len(data[0]) - 1
for i in range(n_columns):
    xlist = list()
    ylist = list()
    for j in range(time):
        ylist.append(data[j][i+1])
        xlist.append(data[j][0])
    plt.plot(xlist, ylist, label=legend[i])
plt.xlabel(xlabel)
plt.title(title)
plt.ylabel(ylabel)
plt.legend()
#plt.plot(data[:][0], data[:][1])

plt.savefig('{}.png'.format(filename[:-4]))
