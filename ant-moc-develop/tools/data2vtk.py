#!/usr/bin/env python3
"""@package vtu2xls

This package provides methods for converting h5 files into vtk files.
It relies on h5py and pyevtk.

Author: An Wang, USTB (wangan@xs.ustb.edu.cn)
"""

import sys
import collections

import numpy as np

def _get_numpy_array(hdf5_group, key):
    """A helper routine to ensure that the data is a proper NumPy array"""

    data_array = np.array(hdf5_group['{}/'.format(key)][...])
    data_array = np.atleast_1d(data_array)
    data_array = data_array.flatten()
    return data_array


def convert_fsr_to_vtk(filepath='./fsr_data.h5'):
    """This routine loads an HDF5 file and convert it to a VTK file

    Parameters
    ----------
    filepath : str
        Filepath for reaction rates HDF5 file (default is './fsr_data.h5')

    """

    # Create a h5py file handle for the file
    import h5py
    f = h5py.File(filepath, 'r')

    # Specify the data source
    domain_type = 'FSR'

    # Check that the file has an 'fsrs' attribute
    if '# fsrs' not in f.attrs:
        print('Unable to load HDF5 file "%s" since it does '
              'not contain an \'# fsrs\' attribute' % filepath)
        sys.exit()

    if domain_type not in f.keys():
        print('Unable to load HDF5 file "%s" since it does '
              'not contain domain type "%s"' % (filepath, domain_type))
        sys.exit()

    # Check that the file has FSR points or centroids
    point_types = ['Points', 'Centroids']
    point_type = ''
    for t in point_types:
        if t in f[domain_type]:
            point_type = t
            break
    if point_type in point_types:
        print('FSR point type is "%s"' % point_type)
    else:
        print('Unable to find any kind of FSR points in HDF5 file "%s"' % filepath)
        sys.exit()

    # Each of the tally type has several datasets
    tally_types = ['Fission RX',
                   'NuFission RX',
                   'Total RX',
                   'Absorption RX',
                   'Scalar Flux',
                   'Fission XS',
                   'NuFission XS',
                   'Total XS',
                   'Absorption XS']

    # Instantiate dictionary to hold FSR data
    fsr_points = {}
    fsr_data = {}
    num_fsrs = int(f.attrs['# fsrs'])

    # Iterate over all domains (e.g., fsrs, tracks) in the HDF5 file
    domain_obj = f[domain_type]
    for group_name in sorted(domain_obj):

        print('Found data for %s "%s"' % (domain_type, str(group_name)))

        # Create shortcut to HDF5 group for this domain
        group_obj = domain_obj[group_name]

        # Read FSR centroids from the file
        if group_name in point_types:
            if group_name == point_type:
                print('Importing data for %s X Y Z' % group_name)
                fsr_points['X'] = _get_numpy_array(group_obj, 'X')
                fsr_points['Y'] = _get_numpy_array(group_obj, 'Y')
                fsr_points['Z'] = _get_numpy_array(group_obj, 'Z')

        # Read FSR volumes from the file
        elif group_name == 'Volumes':
            print('Importing data for', group_name)
            fsr_data[group_name] = _get_numpy_array(domain_obj, group_name)

        # Read reaction rates from the file
        elif group_name in tally_types:
            for energy_name in sorted(group_obj):

                # Set the name of the data array
                array_name = group_name
                if 'sum' != energy_name:
                  array_name = array_name + ' ' + energy_name

                print('Importing data for', array_name)

                fsr_data[array_name] = _get_numpy_array(group_obj, energy_name)

        else:
            print('Unrecognized data "%s"' % group_name)

    # Use PyEVTK to generate the VTK file
    from pyevtk.hl import pointsToVTK

    # Write the data to a VTK file
    head_str, sep_str, tail_str = filepath.rpartition('.')
    filepath = pointsToVTK(head_str, fsr_points['X'], fsr_points['Y'], fsr_points['Z'], fsr_data)

    print('Finished writing VTK data to "%s"' % filepath)
    print('The number of FSRs is', num_fsrs)

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print('The path to an HDF5 file must be specified')
        sys.exit()

    file = sys.argv[1]

    # Read data from the file and produce a VTK file
    convert_fsr_to_vtk(file)
