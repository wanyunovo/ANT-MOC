#!/usr/bin/env python3
"""@package vtu2xls

This package provides methods for converting vtu files into xls files.

Requirements:
    xlwt

Author: An Wang, USTB (wangan@xs.ustb.edu.cn)
"""

import sys
import xml.etree.ElementTree as ET
import xlwt

def convert_vtk_to_xls(file = 'reaction_rates.vtu', z = -1):
    """Converts a VTK formatted mesh into an XLS file.
    This method takes the output of ANT-MOC as its argument and
    extracts data at a specified z-section, which is then dumped
    into an XLS file.

    Arguments:

        file: relative path to the VTK file

        z: index to the section from which we extract data
           If z is less than 0 or is larger than the number of
           z-sections in the file, the reaction rate at point
           (x,y) will be the sum of the reaction reates at
           (x,y,0) to (x,y,nz), where nz is the total number
           of z-sections. This is called flattend.

    """

    # Parses the XML file and returns root
    tree = ET.parse(file)
    root = tree.getroot()

    # Finds the extent of the mesh and the number of lattice cells
    piece = root.find('.//Piece')
    extent = piece.get('Extent').split(' ')
    num_cells_xy = int(piece.get('NumberOfCellsXY'))

    # Computes the dimensions of the mesh
    nx = int(extent[0])
    ny = int(extent[1])
    nz = int(extent[2])

    # Checks the arguments
    is_flattened = (z < 0 or z >= nz)

    print('Mesh dimensions = [%d, %d, %d]' % (nx, ny, nz))
    if is_flattened:
        print('Summing up reaction rates along the z-axis')
    else:
        print('Extracting data at the section where z = ', z)

    # Gets all of the valid indices
    cell_data = piece.find('CellData')
    valid_indices = []
    for data_array in cell_data.findall('DataArray'):
        data_name = data_array.get('Name')
        if data_name == 'Valid Indices':
            data = data_array.text.split()
            valid_indices = [int(i) for i in data]

    # Processes data arrays one by one
    # This is done by first extract an array of data from a XML node
    # and then write it down to an XLS worksheet
    work_book = xlwt.Workbook()
    for data_array in cell_data.findall('DataArray'):
        data_name = data_array.get('Name')

        # Skips cross-sections and non-volume-averaged values
        if data_name.find('Avg Scalar Flux') < 0 and data_name.find('Avg Fission RX') < 0:
            print('Skipping "%s"' % data_name)
            continue

        print('Importing data for "%s"' % data_name)

        # Creates a worksheet
        sheet = work_book.add_sheet(data_name)

        # Reads reaction rates.
        data = data_array.text.split()

        # We have to reorder the data to make sure that the starting point
        # is the upper left corner at z=0. Remember that (i,j) is the position
        # in the tally mesh and should be mapped into the Cartesian coordinate
        # system.
        # Iterates through cells on a specified z-section
        for count in range(num_cells_xy):
            # Sums up reaction rates along the z-axis
            value = 0.
            if is_flattened:
                for k in range(nz):
                    lookup_count = count + (nz - 1 - k) * num_cells_xy
                    value += float(data[lookup_count])
            else:
                lookup_count = count + (nz - 1 - z) * num_cells_xy
                value = float(data[lookup_count])

            # Computes the lattice cell position in the XLS file
            indx = valid_indices[lookup_count]
            i = int(indx % nx)
            j = int(indx // nx % ny)
            row = ny - 1 - j
            column = i
            sheet.write(row, column, value)

    # Saves the workbook
    head, sep, tail = file.rpartition('.')
    if is_flattened:
        file = head + '_flattened.xls'
    else:
        file = head + '.xls'
    work_book.save(file)

    print('Finished writing worksheets to "%s"' % file)


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print('The path to a VTK file must be specified')
        sys.exit()
    else:
        file = sys.argv[1]

    if len(sys.argv) > 2:
        z = int(sys.argv[2])
    else:
        z = -1

    convert_vtk_to_xls(file, z)
