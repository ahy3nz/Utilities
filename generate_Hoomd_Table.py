import os
import numpy as np
import sys
"""
Convert Martini (or any gromacs) LJ parameters 
to tabulated form for HOOMD simulation.
table_energy [=] 0.1 kcal/mol
table_dist [=] 6 angstrom
"""
# These are conversion factors to convert from gromacs 
# units (kJ mol-1, nm, kj mol-1  nm-1) to table units
ENERGY_GMX_TABLE = 0.4184
DISTANCE_GMX_TABLE = 0.6
FORCE_GMX_TABLE = 4.184/6

def convert_unit(self, value=value, tag=None):
    """ Convert a gromacs unit to tabular unit

    Parameters
    ---------
    value : float
        value to be converted
    tag : str
        specify fundamental unit

    Returns
    -------
    converted_unit : float
        Value, converted to tabular units

"""
    if tag == "energy":
        converted_unit = value * ENERGY_GMX_TABLE
    elif tag == "force":
        converted_unit = value * FORCE_GMX_TABLE
    elif tag == "distance":
        converted_unit = value * DISTANCE_GMX_TABLE
    else:
        sys.exit("Specify a tag (energy, force, distance) for unit conversion")
    return converted_unit

def calc_LJ_energy(self, r = 1, C12 = 0, C6 = 0):
    """ Compute energy according to LJ potential

    Parameters
    ----------
    r : float
        distance (nm)
    C12 : float
        C12 parameter (see section 4.1.1 in gmx manual)
    C6 : float
        C6 parameter (see section 4.1.1 in gmx manual)

    Returns
    -------
    energy : float
        energy (see section 4.1.1 in gmx manual)

    Notes
    -----
    This LJ energy calculation is based off gromacs units.
    See section 4.1.1 in the gromacs manual
    """
    r12 = r ** 12
    r6 = r ** 6
    energy = (C12/r12) - (C12/r6)
    return energy

def calc_LJ_force(self, r = 1, C12 = 0, C6 = 0):
    """ Compute force according to derivative of LJ potential

    Parameters
    ----------
    r : float
        distance (nm)
    C12 : float
        C12 parameter (see section 4.1.1 in gmx manual)
    C6 : float
        C6 parameter (see section 4.1.1 in gmx manual)

    Returns
    -------
    force : float
        force (see section 4.1.1 in gmx manual)

    Notes
    -----
    This LJ force calculation is based off gromacs units.
    See section 4.1.1 in the gromacs manual
    """
    r13 = r ** 13
    r7 = r ** 7
    force = abs((12*C12/r13) - (6*C6/r7))
    return force



