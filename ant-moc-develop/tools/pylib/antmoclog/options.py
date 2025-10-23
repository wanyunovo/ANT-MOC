#!/usr/bin/env python3

"""Class for antmoclog options.

Classes: LogOptions

Authors: An Wang, USTB (wangan.cs@gmail.com)

Date: 2020/11/16

"""

import multiprocessing
from antmocutils.baseoptions import BaseOptions, BaseOpt
import antmoclog


class LogOptions(BaseOptions):
    """Singleton class Options for antmoclog.

    This class can parse command line options. Each of the BaseOpt object can
    be accessed by either its full name or short name.

    Command line arguments are supposed to be passed to method parse(...).
    e.g.
        import sys
        options = LogOptions()
        options.parse(sys.argv[1:])

    Examples
    --------

    Instantiate LogOptions
    >>> options = LogOptions()

    >>> isinstance(options, LogOptions)
    True

    Get option value of "h"
    >>> print(options["h"].value)
    False

    Get option value of "output"
    >>> print(options["output"].value)
    antmoc.records.csv

    Add an option to the object
    >>> options.add_opts([ BaseOpt(name="test") ])

    >>> options["test"].value = "hello world"

    Print the option value of "test"
    >>> print(options["test"].value)
    hello world

    """

    # Uncomment this block to make it an singleton
    #instance = None

    #def __new__(cls):
    #    if cls.instance is None:
    #        cls.instance = super(LogOptions, cls).__new__(cls)
    #    return cls.instance

    def __init__(self):
        super().__init__()

        # Append available options
        self.add_opts([
            BaseOpt(name="help",        shortname="h", isbool=True, doc="Show this message"),

            BaseOpt(name="help-fields", shortname="",  isbool=True, doc="Show a list of available fields"),

            BaseOpt(name="nprocs",      shortname="n", default=multiprocessing.cpu_count(),
                    doc="Number of processes to run"),

            BaseOpt(name="cache",       isbool=True,
                    doc="Cache log files in memory. Not to cache files will save memory but read files everytime they are needed."),

            BaseOpt(name="filenames",   shortname="f", default=["log/**/*.log"],  sep=" ", valspec="GLOB1[ GLOB2 ...]",
                    doc="Log filename patterns"),

            BaseOpt(name="fileformat",  shortname="x", default="txt", valspec="STR",
                    doc="Log file format, which may used as the default format"),

            BaseOpt(name="saveto",      valspec="DIR",
                    doc="Save all of the log files in json to a directory"),

            BaseOpt(name="specs",       shortname="e", default=[".*"],     sep=" ", valspec="SPEC1[ SPEC2 ...]",
                    doc="Field specs used to extract log records. A spec is a string with three parts: field name, op, and value. The op and value are used to filter log files."),

            BaseOpt(name="output",      shortname="o", default="antmoc.records.csv", valspec="FILE",
                    doc="CSV file to store extracted records"),

            BaseOpt(name="sortby",      shortname="", valspec="FIELD",
                    doc="Sort the results by a field"),

            BaseOpt(name="sep",         shortname="",  default=" | ",
                    doc="Separator between field values used for formatted printing"),

            BaseOpt(name="truncate",    shortname="",  isbool=True,
                    doc="Truncate the output file"),

            BaseOpt(name="summary",     shortname="",  isbool=True,
                    doc="Print a summary for the query rather than list the results"),
        ])

    def help_fields(self):
        return antmoclog.fields.help_fields()


if __name__ == "__main__":
    import doctest
    doctest.testmod()

