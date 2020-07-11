import os 
import numpy as np
from Utilities.bilayer import *

def determine_components(filename = ""):
    """ Based on the filename return the number of components"""
    #simulation = filename.split("/")[-1]
    prefix = filename.split("_")
    n_components = len(prefix[0])
    return n_components


def read_data_file(filename = "default.dat"):
    """ Read a data file of bilayer structure data
    
    parameters
    ---------
    filename : str
    n_conpoments : int

    Returns
    -------
    bilayer_mixture : bilayer

    Notes
    -----
    col0: APL
    col1: APT
    col2: height
    col3: angle
    col4: oh offset (if applicable)
    col5: ffa offset (if applicable)
    etc...
    col4+n_components-1: interdigitation
    """

    try:
        data = np.loadtxt(filename, ndmin = 2)
    except ValueError:
        print("Error reading file: {}".format(filename))
    filename = filename.split("/")
    filename = filename[-1]
    name = filename[:-4]

    n_components = determine_components(filename = filename)

    # Gather properties from each line

    # Exception if only one data entry
    if(data.ndim == 1):
        apl = data[0, 0]
        apt = data[0, 1]
        height = data[0, 2]
        tilt_angle = data[0, 3]

        idig_col = 4 + n_components - 1
        try:
            idig = data[0, idig_col]
        except IndexError:
            print("Error reading file: {}".format(filename))

        # Create a matrix of all the offsets, where each row is a single simulation
        # And each column is a component
        offsets = []
        for i in np.arange(4, idig_col):
            offsets.append(data[0, i])


    else:
        apl = data[:, 0]
        apt = data[:, 1]
        height = data[:, 2]
        tilt_angle = data[:, 3]

        idig_col = 4 + n_components - 1
        idig = data[:, idig_col]

        # Create a matrix of all the offsets, where each row is a single simulation
        # And each column is a component
        offsets = []
        for i in np.arange(4, idig_col):
            offsets.append(data[:, i])

    # Compute average and standard errors
    apl_avg = np.mean(apl)
    apl_std = np.std(apl)
    apt_avg = np.mean(apt)
    apt_std = np.std(apt)
    height_avg = np.mean(height)
    height_std = np.std(height)
    tilt_angle_avg = np.mean(tilt_angle)
    tilt_angle_std = np.std(tilt_angle)
    idig_avg = np.mean(idig)
    idig_std = np.std(idig)
    offsets_avg = np.mean(offsets, axis=1)
    offsets_std = np.std(offsets, axis=1)

    # Save as a bilayer object
    bilayer_mixture = bilayer(name = name, apl = apl_avg, apl_std = apl_std,
            apt = apt_avg, apt_std = apt_std, height = height_avg, height_std = height_std,
            tilt_angle = tilt_angle_avg, tilt_angle_std = tilt_angle_std, 
            idig = idig_avg, idig_std = idig_std, offsets = offsets_avg, 
            offsets_std = offsets_std, n_components = n_components)

    return bilayer_mixture

def collect_data_files(curr_path = ""):
    """ Given a directory, obtain all the .dat files

    Parameters
    =======--
    curr_path : str
        Path to search for data files

    Returns
    -------
    data_files : list
        List of paths to data files 
        """

    data_files = [f for f in os.listdir(curr_path) if '.dat' in f]
    return data_files
