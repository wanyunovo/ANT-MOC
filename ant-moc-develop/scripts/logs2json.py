#!/usr/bin/env python3

"""An example for exporting log files to json files.

This script takes options "--saveto".
"""

# sys.argv
import sys
import timeit

# package antmoclog.options : class LogOptions
# package antmoclog.data    : class LogDB
from antmoclog import options, data

# Parser for command line arguments
options = options.LogOptions()

# Reset the default value
options["saveto"].default = "antmoc_logdb/"

# Parse command line arguments
options.parse(sys.argv[1:])

# Check if we should print a help message
if options["help"].value:
    options.help()
    exit(1)

# Initialize a LogDB
logdb = data.LogDB()
logdb.setup(options)

t_start = timeit.default_timer()

# Read log files and dump them into json files
logdb.dump_to(options["saveto"].value)

t_stop = timeit.default_timer()
print("Time = {:.3f} s".format(t_stop - t_start))

