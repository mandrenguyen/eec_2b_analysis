#!/bin/bash

# Source before using ROOT/hadd on files containing RooUnfold objects.
case $- in
  *u*) _roounfold_had_nounset=1 ;;
  *) _roounfold_had_nounset=0 ;;
esac
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK=${SCRIPT_DIR}
ROOUNFOLD_DIR=${WORK}/RooUnfold_build_test
if [ ! -f "${ROOUNFOLD_DIR}/build/libRooUnfold.so" ] && \
   [ -f "${WORK}/../RooUnfold_build_test/build/libRooUnfold.so" ]; then
  ROOUNFOLD_DIR=${WORK}/../RooUnfold_build_test
fi
ROOUNFOLD_INC=${ROOUNFOLD_DIR}/src/src
ROOUNFOLD_BUILD=${ROOUNFOLD_DIR}/build

if [ ! -f "${ROOUNFOLD_BUILD}/libRooUnfold.so" ]; then
  echo "ERROR: missing RooUnfold library: ${ROOUNFOLD_BUILD}/libRooUnfold.so" >&2
  return 1 2>/dev/null || exit 1
fi

unset ROOT_INCLUDE_PATH
unset CPLUS_INCLUDE_PATH
unset CPATH
# A loaded CMSSW runtime collides with the LCG view.
unset LD_LIBRARY_PATH
unset PYTHONPATH

# LCG setup scripts are not nounset-safe.
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_106a/x86_64-el9-gcc14-opt/setup.sh
if [ ${_roounfold_had_nounset} -eq 1 ]; then
  set -u
else
  set +u
fi
unset _roounfold_had_nounset

export ROOT_INCLUDE_PATH=${ROOUNFOLD_INC}:${ROOUNFOLD_BUILD}
export LD_LIBRARY_PATH=${ROOUNFOLD_BUILD}:${LD_LIBRARY_PATH:-}
export ROOUNFOLD_DIR ROOUNFOLD_INC ROOUNFOLD_BUILD
