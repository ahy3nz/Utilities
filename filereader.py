import os 
import numpy as np

def _determine_components(filename = filename):
    """ Based on the filename return the number of components"""
    prefix = filename.split("_")
    n_components = len(prefix)


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

    data = np.loadtxt(filename)
    name = filename[:-4]

    n_components = _determine_components(filename = filename)

    # Gather properties from each line
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
    apt_std = np.std(apl)
    apt_avg = np.mean(apt)
    apt_std = np.std(apt)
    height_avg = np.mean(height)
    height_std = np.std(height)
    tilt_angle_avg = np.mean(tilt_angle)
    tilt_angle_std = np.std(tilt_angle)
    idig_avg = np.mean(idig)
    idig_std = np.mean(idig)
    offsets_avg = np.mean(offsets, axis=1)
    offsets_std = np.std(offsets, axis=1)

    # Save as a bilayer object
    bilayer_mixture = bilayer(name = name, apl = apl_avg, apl_std = apl_std,
            apt = apt_avg, apt_std = apt_std, height = height_avg, height_std = height_std,
            tilt_angle = tilt_angle_avg, tilt_angle_std = tilt_angle_std, 
            idig = idig_avg, idig_std = idig_std, offsets = offsets_avg, 
            offsets_std = offsets_std)

    return bilayer_mixture
