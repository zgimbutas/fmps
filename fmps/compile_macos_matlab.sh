#!/bin/sh
#
#  Build the OpenMP MEX file for MATLAB on Apple Silicon (maca64).
#  Run from fmps/. MATLAB's bin goes first on PATH so that its mex is
#  used (TeX Live also installs a mex), and the matlab-maca64-openmp
#  target links MATLAB's own libomp: Homebrew libgomp SIGILLs in
#  MATLAB's thread teardown. Uses the newest /Applications/MATLAB_*.app,
#  as mwrap.inc does.
#
MATLAB_APP=$(ls -d /Applications/MATLAB_*.app 2>/dev/null | sort | tail -1)
if [ -z "$MATLAB_APP" ]; then
  echo "compile_macos_matlab.sh: no /Applications/MATLAB_*.app found" >&2
  exit 1
fi
PATH=.:$MATLAB_APP/bin:$PATH
make TARGET=matlab-maca64-openmp -f makefile.mwrap -j8 clean
make TARGET=matlab-maca64-openmp -f makefile.mwrap -j8 || exit 1
make TARGET=matlab-maca64-openmp -f makefile.mwrap -j8 clean
