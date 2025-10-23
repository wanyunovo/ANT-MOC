#!/bin/sh

# add this script to .gitconfig or GIT_DIR/config, where
# GIT_DIR is usually .git
#
# [diff "HDF5"]
#     textconv = /path/to/h5dump-wrapper.sh
#
# you can also take advantage of git config
# git config diff.HDF5.textconv /path/to/h5dump-wrapper.sh
#
# diff is called by git with 7 parameters:
# path old-file old-hex old-mode new-file new-hex new-mode

h5dump "$1" | tail -n +2
