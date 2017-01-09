from __future__ import print_function
import mdtraj as mdtraj
import sys
import scipy.integrate as integrate
import os
from optparse import OptionParser
import pdb
import ipdb
import numpy as np
import matplotlib
import collections
from collections import OrderedDict
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def calc_APL(traj, n_lipid):
    ''' 
    Input: Trajectory and number of lipids
    Compute areas by looking at x and y unit cell lengths
    Return: array of area per lipids (n_frame x 1) [Angstrom]
    '''
    area = 100 * traj.unitcell_lengths[:, 0] * traj.unitcell_lengths[:, 1]
    areaavg = np.mean(area)
    areastd = np.std(area)/(len(area)**0.5)
    apl_list = np.eye(traj.n_frames, 1)
    apl_list[:,0] = area[:]/(n_lipid/2)
    apl_avg = areaavg/(n_lipid/2)
    apl_std = areastd/(n_lipid/2)
    return (apl_avg, apl_std, apl_list)

def get_lipids(topol):
    ''' Input a topology object
        Iterate through each atom
        If the atom's residue isn't water (so a lipid), add it to the dictionary
        lipid_dict is a dictionary mapping residue indices to a list of respective atom indices
        headgroup_dict is a dictionary mapping moleculetypes (DSPC, alc12, etc) to a list of 
        indices of the headgroup
    '''
    # Dictionary of resname keys that map to atom list values
    lipid_dict= OrderedDict()
    headgroup_dict = OrderedDict()
    for i in topol.select('all'):
        atom_i = topol.atom(i)
        if not atom_i.residue.is_water:
            residue_i = atom_i.residue
            resname = atom_i.residue.name
            # If the lipid_dict already has the residue key, append it
            if residue_i.index in lipid_dict:
                lipid_dict[residue_i.index].append(i)
            # If the lipid_dict doesn't have the residue key, make a list and append it
            else:
                lipid_dict[residue_i.index] = list()
                lipid_dict[residue_i.index].append(i)
            # Figure out head groups
            if 'DSPC' in resname:
                if 'DSPC' in headgroup_dict:
                    if 'P' in atom_i.name or 'OM' in atom_i.name or 'OA' in atom_i.name:
                        headgroup_dict['DSPC'].append(i)
                else:
                    if 'P' in atom_i.name or 'OM' in atom_i.name or 'OA' in atom_i.name:
                        headgroup_dict['DSPC'] = list()
                        headgroup_dict['DSPC'].append(i)
            elif 'DPPC' in resname:
                if 'DPPC' in headgroup_dict:
                    if 'P' in atom_i.name or 'OM' in atom_i.name or 'OA' in atom_i.name:
                        headgroup_dict['DPPC'].append(i)
                else:
                    if 'P' in atom_i.name or 'OM' in atom_i.name or 'OA' in atom_i.name:
                        headgroup_dict['DPPC'] = list()
                        headgroup_dict['DPPC'].append(i)
            elif 'ISIS' in resname:
                if 'ISIS' in headgroup_dict:
                    if 'O' in atom_i.name or 'OE' in atom_i.name:
                        headgroup_dict['ISIS'].append(i)
                else:
                    if 'O' in atom_i.name or 'OE' in atom_i.name:
                        headgroup_dict['ISIS'] = list()
                        headgroup_dict['ISIS'].append(i)
            elif 'SS' in resname:
                print("SS headgroups not included")
            elif 'acd16' in resname:
                if 'acd16' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['acd16'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['acd16'] = list()
                        headgroup_dict['acd16'].append(i)

            elif 'acd22' in resname:
                if 'acd22' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['acd22'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['acd22'] = list()
                        headgroup_dict['acd22'].append(i)
            elif 'alc12' in resname:
                if 'alc12' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc12'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc12'] = list()
                        headgroup_dict['alc12'].append(i)
            elif 'alc14' in resname:
                if 'alc14' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc14'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc14'] = list()
                        headgroup_dict['alc14'].append(i)
            elif 'alc16' in resname:
                if 'alc16' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc16'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc16'] = list()
                        headgroup_dict['alc16'].append(i)
            elif 'alc18' in resname:
                if 'alc18' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc18'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc18'] = list()
                        headgroup_dict['alc18'].append(i)
            elif 'alc20' in resname:
                if 'alc20' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc20'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc20'] = list()
                        headgroup_dict['alc20'].append(i)
            elif 'alc22' in resname:
                if 'alc22' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc22'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc22'] = list()
                        headgroup_dict['alc22'].append(i)
            elif 'alc24' in resname:
                if 'alc24' in headgroup_dict:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc24'].append(i)
                else:
                    if 'CH' not in atom_i.name:
                        headgroup_dict['alc24'] = list()
                        headgroup_dict['alc24'].append(i)

    return lipid_dict, headgroup_dict

def get_lipid_tails(topol, lipid_dict):
    ''' Input topology the lipid dictionary
        Look at the atoms in the lipid dictionary
        For each atom, get the index, (shifted to zero), residue name and residue index
        For that particular residue name, figure out if it belongs in a tail and add it
        If there are multiple tails for that residue, denote differences with 'a' and 'b'
        So Lipida and Lipidb are different keys in the lipid_tails dict
        But respective values correspond to that tail
        Return lipid_tails, a dictionary mapping each lipid tail to its atoms
    '''

    # Get an atom
    # Get that atom's residue
    # Shift atom indices by the index of the first atom in the residue 0
    # Based on the residue, check if that atom's shifted index falls within the tail range
    lipid_tails = OrderedDict()
    for lipid in lipid_dict.keys():
        lipid_atoms = lipid_dict[lipid]
        for atom_index in lipid_atoms:
            shifted_index = atom_index - lipid_atoms[0]
            atom_i = topol.atom(atom_index)
            resname = atom_i.residue.name
            resindex = atom_i.residue.index
            # This might need improvement, right now hard coding lipid tail definitions
            # Looking at 12 carbons after the headgroup
            if 'DSPC' in resname:
                if 14 == shifted_index or 16 <= shifted_index <= 26:
                    #if (resname + str(resindex) + 'a') in lipid_tails:
                    if (str(resindex) + 'a') in lipid_tails:
                        #lipid_tails[(resname + str(resindex) + 'a')].append(atom_index)
                        lipid_tails[(str(resindex) + 'a')].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex) + 'a')] = list()
                        #lipid_tails[(resname + str(resindex) + 'a')].append(atom_index)
                        lipid_tails[(str(resindex) + 'a')] = list()
                        lipid_tails[(str(resindex) + 'a')].append(atom_index)


                elif 35 == shifted_index or 37 <= shifted_index <= 47:
                    #if (resname + str(resindex) + 'b') in lipid_tails:
                    if (str(resindex) + 'b') in lipid_tails:
                        #lipid_tails[(resname + str(resindex) + 'b')].append(atom_index)
                        lipid_tails[(str(resindex) + 'b')].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex) + 'b')] = list()
                        #lipid_tails[(resname + str(resindex) + 'b')].append(atom_index)
                        lipid_tails[(str(resindex) + 'b')] = list()
                        lipid_tails[(str(resindex) + 'b')].append(atom_index)


            elif 'DPPC' in resname:
                if 14 == shifted_index or 16 <= shifted_index <= 26:
                   # if (resname + str(resindex) + 'a') in lipid_tails:
                    if (str(resindex) + 'a') in lipid_tails:
                        #lipid_tails[(resname + str(resindex) + 'a')].append(atom_index)
                        lipid_tails[(str(resindex) + 'a')].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex) + 'a')] = list()
                        #lipid_tails[(resname + str(resindex) + 'a')].append(atom_index)
                        lipid_tails[(str(resindex) + 'a')] = list()
                        lipid_tails[(str(resindex) + 'a')].append(atom_index)


                elif 33 == shifted_index or 35 <= shifted_index <= 45:
                    #if (resname + str(resindex) + 'b') in lipid_tails:
                    if ( str(resindex) + 'b') in lipid_tails:
                        #lipid_tails[(resname + str(resindex) + 'b')].append(atom_index)
                        lipid_tails[(str(resindex) + 'b')].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex) + 'b')] = list()
                        #lipid_tails[(resname + str(resindex) + 'b')].append(atom_index)
                        lipid_tails[(str(resindex) + 'b')] = list()
                        lipid_tails[(str(resindex) + 'b')].append(atom_index)

            elif 'ISIS' in resname:
                if 6 <= shifted_index <= 17:
                    #if (resname + str(resindex) + 'a') in lipid_tails:
                    if ( str(resindex) + 'a') in lipid_tails:
                        #lipid_tails[(resname + str(resindex) + 'a')].append(atom_index)
                        lipid_tails[(str(resindex) + 'a')].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex) + 'a')] = list()
                        #lipid_tails[(resname + str(resindex) + 'a')].append(atom_index)
                        lipid_tails[( str(resindex) + 'a')] = list()
                        lipid_tails[( str(resindex) + 'a')].append(atom_index)

                elif 19 == shifted_index or 21 <= shifted_index <= 31:
                    #if (resname + str(resindex) + 'b') in lipid_tails:
                    if (str(resindex) + 'b') in lipid_tails:
                        lipid_tails[(str(resindex) + 'b')].append(atom_index)
                        #lipid_tails[(resname + str(resindex) + 'b')].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex) + 'b')] = list()
                        #lipid_tails[(resname + str(resindex) + 'b')].append(atom_index)
                        lipid_tails[(str(resindex) + 'b')] = list()
                        lipid_tails[(str(resindex) + 'b')].append(atom_index)

            elif 'SS' in resname:
                if 5 <= shifted_index <= 16:
                   # if (resname + str(resindex) + 'a') in lipid_tails:
                    if ( str(resindex) + 'a') in lipid_tails:
                        #lipid_tails[(resname + str(resindex) + 'a')].append(atom_index)
                        lipid_tails[(str(resindex) + 'a')].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex) + 'a')] = list()
                        #lipid_tails[(resname + str(resindex) + 'a')].append(atom_index)
                        lipid_tails[(str(resindex) + 'a')] = list()
                        lipid_tails[(str(resindex) + 'a')].append(atom_index)

                elif 18 == shifted_index or 20 <= shifted_index <= 30:
                    #if (resname + str(resindex) + 'b') in lipid_tails:
                    if ( str(resindex) + 'b') in lipid_tails:
                        #lipid_tails[(resname + str(resindex) + 'b')].append(atom_index)
                        lipid_tails[(str(resindex) + 'b')].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex) + 'b')] = list()
                        #lipid_tails[(resname + str(resindex) + 'b')].append(atom_index)
                        lipid_tails[(str(resindex) + 'b')] = list()
                        lipid_tails[(str(resindex) + 'b')].append(atom_index)


            elif 'acd16' in resname:
                #if shifted_index == 16 or 4 <= shifted_index <= 14:
                if 0 <= shifted_index <= 14:
                    #if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))] = list()
                        lipid_tails[(str(resindex))].append(atom_index)

            elif 'acd22' in resname:
                #if shifted_index == 22 or 10 <= shifted_index <= 20:
                if 0 <= shifted_index <= 20:
                   # if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))] = list()
                        lipid_tails[(str(resindex))].append(atom_index)


            elif 'alc12' in resname:
                if 0 <= shifted_index <= 11:
                    #if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex) )].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))] = list()
                        lipid_tails[(str(resindex))].append(atom_index)

            elif 'alc14' in resname:
                if 0 <= shifted_index <= 13:
                    #if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))] = list()
                        lipid_tails[(str(resindex))].append(atom_index)


            elif 'alc16' in resname:
                if 0 <= shifted_index <= 15:
                    #if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))] = list()
                        lipid_tails[(str(resindex))].append(atom_index)

            elif 'alc18' in resname:
                if 0 <= shifted_index <= 17:
                    #if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex) )].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[( str(resindex))] = list()
                        lipid_tails[( str(resindex))].append(atom_index)

            elif 'alc20' in resname:
                if 0 <= shifted_index <= 19:
                    #if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[( str(resindex))] = list()
                        lipid_tails[( str(resindex))].append(atom_index)

            elif 'alc22' in resname:
                if 0 <= shifted_index <= 21:
                    #if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))] = list()
                        lipid_tails[(str(resindex))].append(atom_index)
            elif 'alc24' in resname:
                if 0 <= shifted_index <= 23:
                    #if (resname + str(resindex)) in lipid_tails:
                    if ( str(resindex)) in lipid_tails:
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))].append(atom_index)
                    else:
                        #lipid_tails[(resname + str(resindex))] = list()
                        #lipid_tails[(resname + str(resindex))].append(atom_index)
                        lipid_tails[(str(resindex))] = list()
                        lipid_tails[(str(resindex))].append(atom_index)

            else:
                print('Lipid {} not incorporated in lipid tail identification'.format(resname))
                sys.exit()
    return lipid_tails

def calc_tilt_angle(traj, topol, lipid_tails):
    ''' 
    Input: Trajectory, topology, dictionary of lipid tails with atom index values
    Compute characteristic vector using eigenvector associated with
    lowest eigenvalue of inertia tensor.
    Compute angle between charactersitic vector and lipid tail,
    adjusted for the first quadrant of  cartesian coordinate space
    Return: array of tilt angles (n_frame x n_lipid_tail)
    '''

    surface_normal = np.asarray([0, 0, 1.0])
    angle_list = []
    angle_list = np.eye(traj.n_frames, len(lipid_tails.keys()))
    index = 0
    for key in lipid_tails.keys():
        lipid_i_atoms = lipid_tails[key]
        traj_lipid_i = traj.atom_slice(lipid_i_atoms)
        director = mdtraj.geometry.order._compute_director(traj_lipid_i)
        lipid_angle = np.rad2deg(np.arccos(np.dot(director, surface_normal)))
        for i,angle in enumerate(lipid_angle):
            if angle >= 90:
                angle = 180- angle
                lipid_angle[i] = angle
        #angle_list.append(lipid_angle)
        angle_list[:,index] = lipid_angle
        index += 1
    angle_avg = np.mean(angle_list)
    angle_std = np.std(angle_list)/(len(angle_list)**0.5)
    return angle_avg, angle_std, angle_list


def calc_APT(apl_list, angle_list, n_tails_per_lipid):
    ''' Input: a matrix of area per lipids (each row is a frame               
        a matrix of tilt angels (each row is a frame, each column is a lipid)
        Return matrix of area per tail (n_frame x n_lipid_tail)
    '''
    # Each element in angle list correspond to a tail, and that element is a row of tilts per frame
    # Each element in apl list is the apl for a frame
    apt_list = angle_list
    apt_list = np.cos(np.deg2rad(angle_list[:,:]))*apl_list[:]/n_tails_per_lipid
    apt_avg = np.mean(apt_list)
    apt_std = np.std(apt_list)/(len(apt_list)**0.5)
    return apt_avg, apt_std, apt_list

def calc_mean(dataset):
    ''' Generic mean calculation 
    of a dataset. Assumes each elemtn in the dataset
    is a simtk Quantity'''

    avg = dataset[0]
    for i in range(1, len(dataset) - 1):
        avg = avg.__add__(dataset[i])
    avg = avg.__truediv__(len(dataset))
    return avg

def calc_stdev(avg, dataset):
    ''' Generic standard deviation calculation
    of a dataset. Assumes each element in the dataset
    is a simtk Quantity'''

    variance = (avg.__sub__(avg)).__pow__(2)
    for val in dataset:
        deviation =  val.__sub__(avg)
        variance = variance.__add__(deviation.__pow__(2))
    return (variance.__div__(len(dataset))).sqrt()

def read_xvg(filename):
    '''Given an xvg file, read the file
    Return the data as a list of lists
    Return the legend as a list '''

    xvgfile = open(filename, 'r')
    xvglines = xvgfile.readlines()
    data = list()
    legend = []
    for i, line in enumerate(xvglines):
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
    return data, legend

def calc_head_distance(traj, topol, head_indices):
    """
    Input: Trajectory, topology, indices of headgroup atoms
    For each frame, compute the average z-coordinate of the headgroup in the top and bot leaflet
    Compute the difference as the headgroup distance
    Return: array of headgroup distances (n_frame x 1)
    """
    mass_top = 0
    mass_bot = 0
    zcoord_top = 0
    zcoord_bot = 0
    atom_counter = 0
    for atom_j in head_indices:
        if 'CH0' in topol.atom(atom_j).name:
            mass_i = 12.01
        elif 'CH1' in topol.atom(atom_j).name:
            mass_i = 13.02
        elif 'CH2' in topol.atom(atom_j).name:
            mass_i = 14.03
        elif 'CH3' in topol.atom(atom_j).name:
            mass_i = 15.04
        elif 'C' in topol.atom(atom_j).name:
            mass_i = 12.01
        elif 'P' in topol.atom(atom_j).name:
            mass_i = 30.97
        elif 'O' in topol.atom(atom_j).name:
            mass_i = 16.00
        if atom_counter < len(head_indices)/2:
            zcoord_top += mass_i * traj.atom_slice([atom_j]).xyz[:,0,2]
            mass_top += mass_i
        else:
            zcoord_bot += mass_i * traj.atom_slice([atom_j]).xyz[:,0,2]
            mass_bot += mass_i
        atom_counter +=1
    zcoord_top = zcoord_top / mass_top
    zcoord_bot = zcoord_bot / mass_bot
    head_dist_list = 10 * abs(zcoord_top - zcoord_bot)
    head_dist_avg = np.mean(head_dist_list)
    head_dist_std = np.std(head_dist_list)
    return head_dist_avg, head_dist_std, head_dist_list

def compute_headgroup_distances(traj, topol, headgroup_dict):
    """
    Input: trajectory, topology, dictionary mapping molecule types to their headgroup indices
    For each molecule type, compute the distance between headgroups
    Return: dictionary of molecule types and their headgroup distances for each frame 
    Dictionary values are lists of (n_frame x 1)
    """
    headgroup_distance_dict = OrderedDict()
    for key in headgroup_dict.keys():
        headgroup_distance_dict[key] = calc_head_distance(traj, topol, headgroup_dict[key])
    return headgroup_distance_dict

def calc_bilayer_height(headgroup_distance_dict):
    """
    Input: Dictionary of molecule types and their headroup distances
    Calculate bilayer height by comparing DSPC or DPPC headgroup distances
    Return: bilayer height average, bilayer height std, bilayer heigh per frame list
    """
    if headgroup_distance_dict['DSPC']:
        dist_list = headgroup_distance_dict['DSPC'][2]
        dist_avg = np.mean(dist_list)
        dist_std = np.std(dist_list)
    elif headgroup_distance_dict['DPPC']:
        dist_list = headgroup_distance_dict['DPPC'][2]
        dist_avg = np.mean(dist_list)
        dist_std = np.std(dist_list)
    else:
        print ('No phosphate groups to compare')

    return (dist_avg, dist_std, dist_list)
def calc_offsets(headgroup_distance_dict):
    """ 
    Input: dictionary of molecule types and their headgroup distances
    Calculate the offsets with respect to the phosphate group
    Return: dictionary of offsets with respect to the phosphate group
    Values are (distance averaged over n_Frames, distance std over n_frame)
    """
    offset_dict = OrderedDict()
    if headgroup_distance_dict['DSPC']:
        use_DSPC = True
        use_DPPC = False
    elif headgroup_distance_dict['DPPC']:
        use_DSPC = False
        use_DPPC = True
    else:
        print('No phosphate groups to compare to')
    for key in headgroup_distance_dict.keys():
        # Determine if reference is DSPC or DPPC
        if use_DSPC:
            offset_list = (headgroup_distance_dict['DSPC'][2] - headgroup_distance_dict[key][2])/2
            offset_avg = np.mean(offset_list)
            offset_std = np.std(offset_list)
            offset_dict[key] = (offset_avg, offset_std)
        elif use_DPPC:
            offset_list = (headgroup_distance_dict['DPPC'][2] - headgroup_distance_dict[key][2])/2
            offset_avg = np.mean(offset_list)
            offset_std = np.std(offset_list)
            offset_dict[key] = (offset_avg, offset_std)
        else:
            print ('No phosphate groups to compare')

    return offset_dict

def calc_nematic_order(traj, lipid_dict):
    top_chains = []
    bot_chains = []
    for i, key in enumerate(lipid_dict.keys()):
        indices = [int(item) for item in lipid_dict[key]]
        if i <= 63:
            top_chains.append(indices)
        else:
            bot_chains.append(indices)
    s2_top = mdtraj.compute_nematic_order(traj, indices=top_chains)
    s2_bot = mdtraj.compute_nematic_order(traj, indices=bot_chains)
    s2_list = (s2_top + s2_bot)/2
    s2_ave = np.mean(s2_list)
    s2_std = np.std(s2_list)
    return s2_ave, s2_std, s2_list

def calc_density_profile(traj, topol, lipid_dict, n_bins = 50):
    """
    Input:  Trajectory, topology, lipid dictionary (mapping residue numbers to atom numbers), # bins
    Compute the density profile of lipids using provided bins
    Return: 1D histogram of densities in various bins for each frame (for top, bottom, and entire)
    """
    # First, divide box into 50 segments, based on z coordinates
    # Z end is based on the max z box dimension
    n_lipid = len(lipid_dict.keys())
    z_end = np.amax(traj.xyz[:,:,2])
    area = np.mean(traj.unitcell_lengths[:, 0] * traj.unitcell_lengths[:, 1])
    # Interval is the last z coordinate divied by bins
    z_interval = z_end/(n_bins)
    v_slice = z_interval * area
    density_profile_top = 1e-40 * np.ones((traj.n_frames, n_bins + 1))
    density_profile_bot = 1e-40 * np.ones((traj.n_frames, n_bins + 1)) 
    density_profile = 1e-40 * np.ones((traj.n_frames, n_bins + 1))
    # Generate windows
    bins = np.linspace(0, z_end, num = n_bins + 1)
    # For each leaflet, count the mass of a slice in the histogram, divided by volume of that slice
    # masses dictionary
    badcount = 0
    for i, key in enumerate(lipid_dict.keys()):
        lipid_i = lipid_dict[key]
        for atom_i in lipid_i:
        # loop through each lipid atom, get the z coordinate (probably an array over time)
            mass_i = get_mass(topol, atom_i)
            z_i = traj.atom_slice([atom_i]).xyz[:,0,2]
            # Row represents hte frame and the elemtn is the window it belongs in
            window_i = np.floor(z_i/z_interval)
            for j, bin_j in enumerate(window_i):
                if i < n_lipid/2:
                    density_profile_bot[j, int(bin_j)] += mass_i
                else:
                    density_profile_top[j, int(bin_j)] += mass_i
                density_profile[j, int(bin_j)] += mass_i
    
    # Divide by volume of slice to get the density
    density_profile_bot /= v_slice
    density_profile_top /= v_slice
    density_profile /= v_slice
    
    # Convert from g/nm3 to kg/m3, and nm to m
    density_profile_bot = 1e24 * density_profile_bot
    density_profile_top = 1e24 * density_profile_top
    density_profile = 1e24 * density_profile
    density_profile_avg = np.mean(density_profile, axis = 0)

    return density_profile, density_profile_avg, density_profile_bot, density_profile_top, bins

def get_mass(topol, atom_i):
    """
    Input: trajectory, atom index
    Mass dictionary is in units of amu
    Return: mass of that atom (g)
    """
    mass_dict = {'O': 15.99940, 'OM': 15.99940, 'OA': 15.99940, 'OE': 15.99940, 'OW': 15.99940,'N': 14.00670,
            'NT': 14.00670, 'NL': 14.00670, 'NR': 14.00670, 'NZ': 14.00670, 'NE': 14.00670, 'C': 12.01100, 
            'CH0': 12.0110, 'CH1': 13.01900, 'CH2': 14.02700, 'CH3': 15.03500, 'CH4': 16.04300, 'CH2r': 14.02700,
            'CR1': 13.01900, 'HC': 1.00800, 'H':  1.00800, 'P': 30.97380, 'CL': 35.45300, 'F': 18.99840, 
            'CL-': 35.45300}
    mass_i = 1.66054e-24 * mass_dict[topol.atom(atom_i).name]
    return mass_i

def calc_interdigitation(density_profile_top, density_profile_bot):
    """
    Input: top and bottom density profile [kg/m3], but units irrelevant since this is dimensionless
    Compute interdigitation according to "Structural Properties.." by Hartkamp (2016)
    Densities based on nm bins, convert to Angstrom
    Return: Interdigation avg, std, and array of interdigation at each frame (n_frame x 1) [A]
    """
    interdig = integrate.simps( (4*density_profile_top*density_profile_bot)/
            (density_profile_top + density_profile_bot)**2)
    interdig *= 0.1
    interdig_avg = np.mean(interdig)
    interdig_std = np.std(interdig)
    return interdig_avg, interdig_std, interdig


# CODE STARTS HERE
parser = OptionParser()
parser.add_option('-f', action="store", type="string", default = 'nopbc.xtc', dest = 'trajfile')
parser.add_option('-c', action="store", type="string", default = 'Stage5_ZCon0.gro', dest = 'grofile')
parser.add_option('-o', action='store', type='string', default = 'BilayerAnalysis', dest = 'outfilename')

(options, args) = parser.parse_args()
trajfile = options.trajfile
grofile = options.grofile
outfilename = options.outfilename

print('Loading trajectory <{}>...'.format(trajfile))
traj = mdtraj.load(trajfile, top=grofile)
topol = traj.topology


# Compute system information
print('Gathering system information <{}>...'.format(grofile))
lipid_dict, headgroup_dict = get_lipids(topol)
lipid_tails = get_lipid_tails(topol, lipid_dict)

n_lipid = len(lipid_dict.keys())
n_lipid_tails = len(lipid_tails.keys())
n_tails_per_lipid = n_lipid_tails/n_lipid


# Vectorized Calculations start here
print('Calculating area per lipid...')
apl_avg, apl_std, apl_list = calc_APL(traj,n_lipid)
print('Calculating tilt angles...')
angle_avg, angle_std, angle_list = calc_tilt_angle(traj, topol, lipid_tails)
print('Calculating area per tail...')
apt_avg, apt_std, apt_list = calc_APT(apl_list, angle_list, n_tails_per_lipid)
print('Calculating nematic order...')
s2_ave, s2_std, s2_list = calc_nematic_order(traj, lipid_dict)
print('Calculating headgroup distances...')
headgroup_distance_dict = compute_headgroup_distances(traj, topol, headgroup_dict)
print('Calculating bilayer height...')
Hpp_ave, Hpp_std, Hpp_list = calc_bilayer_height(headgroup_distance_dict)
print('Calculating component offsets...')
offset_dict = calc_offsets(headgroup_distance_dict)
print('Calculating density profile...')
density_profile, density_profile_avg, density_profile_top, density_profile_bot, bins = \
    calc_density_profile(traj, topol, lipid_dict)
print('Calculating interdigitation...')
interdig_avg, interdig_std, interdig_list = calc_interdigitation(density_profile_top, density_profile_bot)

# Printing properties

print('Outputting to <{}>...'.format(outfilename))
outfile = open((outfilename + '.txt'),'w')
outpdf = PdfPages((outfilename+'.pdf'))
outfile.write('{:<20s}: {}\n'.format('Trajectory',trajfile))
outfile.write('{:<20s}: {}\n'.format('Structure',grofile))
outfile.write('{:<20s}: {}\n'.format('# Frames',traj.n_frames))
outfile.write('{:<20s}: {}\n'.format('Lipids',n_lipid))
outfile.write('{:<20s}: {}\n'.format('Tails',n_lipid_tails))
outfile.write('{:<20s}: {} ({})\n'.format('APL (A^2)',apl_avg, apl_std))
outfile.write('{:<20s}: {} ({})\n'.format('APT (A^2)',apt_avg, apt_std))
outfile.write('{:<20s}: {} ({})\n'.format('Bilayer Height (A)',Hpp_ave, Hpp_std))
outfile.write('{:<20s}: {} ({})\n'.format('Tilt Angle', angle_avg, angle_std))
outfile.write('{:<20s}: {} ({})\n'.format('S2', s2_ave, s2_std))
outfile.write('{:<20s}: {} ({})\n'.format('Interdigitation (A)', interdig_avg, interdig_std))
for key in offset_dict.keys():
    outfile.write('{:<20s}: {} ({})\n'.format
            ((key + ' offset (A)'), offset_dict[key][0], offset_dict[key][1]))
outfile.write('{:<20s}: {} ({})\n'.format(
    'Leaflet 1 Tilt Angle', np.mean(angle_list[:, 0 :int(np.floor(n_lipid_tails/2))]),
    np.std(angle_list[:, 0 :int(np.floor(n_lipid_tails/2))])))
outfile.write('{:<20s}: {} ({})\n'.format(
    'Leaflet 2 Tilt Angle', np.mean(angle_list[:, int(np.floor(n_lipid_tails/2)):len(angle_list[0])]), 
    np.std(angle_list[:, int(np.floor(n_lipid_tails/2)):len(angle_list[0])])))

# Plotting
fig1 = plt.figure(1)
plt.subplot(3,2,1)
plt.plot(apl_list)
plt.title('APL')

plt.subplot(3,2,2)
plt.plot(np.mean(angle_list, axis=1))
plt.title('Tilt Angle ($^o$)')

plt.subplot(3,2,3)
plt.plot(np.mean(apt_list,axis=1))
plt.title('APT')

plt.subplot(3,2,4)
plt.plot(Hpp_list)
plt.title('H$_{PP}$')

plt.subplot(3,2,5)
plt.plot(s2_list)
plt.title('S2')

plt.subplot(3,2,6)
plt.plot(interdig_list)
plt.title('Interdigitation (A)')

plt.tight_layout()
outpdf.savefig(fig1)
plt.close()

fig2 = plt.figure(2)
plt.subplot(2,1,1)
plt.plot(bins,density_profile_avg)
plt.xlabel('Depth (nm)')
plt.title('Density Profile (kg m$^{-3}$)')


plt.subplot(2,1,2)
plt.hist(np.mean(angle_list[:, 0 : int(np.floor(n_lipid_tails/2))], axis = 0), bins = 50,  
        alpha = 0.5, facecolor = 'blue', normed = True)
plt.hist(np.mean(angle_list[:, int(np.floor(n_lipid_tails/2)) : len(angle_list[0])], axis = 0), bins = 50,  
        alpha = 0.5, facecolor = 'red', normed = True)
plt.title('Angle Distribution by Leaflet')
plt.xlabel('Angle ($^o$)')

plt.tight_layout()
outpdf.savefig(fig2)
plt.close()
outpdf.close()

print('**********')
print('{:^10s}'.format('Done'))
print('**********')
