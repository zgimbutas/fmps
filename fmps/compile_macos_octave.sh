#!/bin/sh
#
#  Build the OpenMP MEX file for GNU Octave on macOS (Homebrew gfortran
#  and libgomp). Run from fmps/. Use it from octave-cli: the octave-gui
#  interpreter thread has too small a stack for the FMM tree code.
#
make TARGET=octave-macos-openmp -f makefile.mwrap -j8 clean
make TARGET=octave-macos-openmp -f makefile.mwrap -j8 || exit 1
make TARGET=octave-macos-openmp -f makefile.mwrap -j8 clean
