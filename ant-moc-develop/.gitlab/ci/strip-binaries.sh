#!/usr/bin/env bash

find -L "${@:-"build"}" -type f -exec readlink -f '{}' \; | \
  xargs file -i | \
  grep 'charset=binary' | \
  grep 'x-executable\|x-archive\|x-sharedlib' | \
  awk -F: '{print $1}' | \
  xargs strip -s
