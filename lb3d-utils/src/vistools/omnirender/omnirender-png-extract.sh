#!/bin/bash
# Use this script to extract the command line options used to generate a png file from the file.
# This requires -m to have been used when rendering.

strings $1 | grep "<argv>" | sed -e "s/.*<argv>\(.*\)<\/argv>.*/\1/"

