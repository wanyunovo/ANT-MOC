#!/usr/bin/env python3

"""Extract records from antmoc log files.

Run this script with "-h" or "--help" to get help.

"""

# sys.argv
import sys

# package antmoclog.options : class LogOptions
# package antmoclog.data    : class LogDB
# package antmoclog.extract : function extractRecords
from antmoclog import options, data, extract

# Parser for command line arguments
options = options.LogOptions()

# Reset default field specs
# Comment out this block to keep the original default ".*", or
# just pass field specs by command line arguments.
options["specs"].default = [
    'File',
    'JobId',
    'CaseName',
    'Refines',
    'Azims',
    'XYSpacing',
    'Polars',
    'ZSpacing',
    'Quadrature',
    'Tolerance',
    'Modules',
    'Domains',
    'Ranks',
    'Threads',
    'ExtFSRs',
    'FSRs',
    'Chains2D',
    'Tracks3D',
    'Segments',
    'Iterations',
    'Keff',
    '.*Time',
    'Mem.*',
]

# Parse command line arguments
options.parse(sys.argv[1:])

# Check if we should print a help message
if options["help"].value:
    options.help()
    exit(1)
elif options["help-fields"].value:
    options.help_fields()
    exit(1)

# Create a LogDB object and set it up with current options
logdb = data.LogDB()
logdb.setup(options)

# Create a tiny extractor and set it up with current options
extractor = extract.TinyExtractor(logdb)
extractor.setup(options)
extractor.extract()

