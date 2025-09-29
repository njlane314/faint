#!/usr/bin/env bash
set -e

# Determine topdir from FAINT or script path
if [ -z "$FAINT" ]; then
  TOPDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
else
  TOPDIR="$FAINT"
fi

LIBDIR="${TOPDIR}/build/lib"
INCDIR="${TOPDIR}/include"
MACRO="${TOPDIR}/scripts/setup_faint.C"

# Ensure library path and includes are visible to ROOT
case "$(uname -s)" in
  Darwin) export DYLD_LIBRARY_PATH="${LIBDIR}:${DYLD_LIBRARY_PATH}" ;;
  *)      export LD_LIBRARY_PATH="${LIBDIR}:${LD_LIBRARY_PATH}" ;;
esac
export ROOT_INCLUDE_PATH="${INCDIR}:${ROOT_INCLUDE_PATH}"

if [ "$(uname -s)" = "Darwin" ]; then
  LIBEXT="dylib"
else
  LIBEXT="so"
fi

# Start ROOT, ensure the helper macro is loaded, then forward any user macro call
root -l -q -e "gROOT->LoadMacro(\"${MACRO}\"); setup_faint(\"${LIBDIR}/libfaint.${LIBEXT}\",\"${INCDIR}\");" "$@"
