#!/bin/bash

# Script to set up TERRA_PATH in order to find the files in ../src

cwd=`pwd`
export TERRA_PATH="${cwd}/../src/?.t;${cwd}/../src/?.rg"
