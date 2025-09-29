#!/usr/bin/env bash
set -e
# Determine topdir from RAREXSEC or script path
if [ -z "$RAREXSEC" ]; then
  TOPDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
else
  TOPDIR="$RAREXSEC"
fi
LIBDIR="${TOPDIR}/build/lib"
INCDIR="${TOPDIR}/include"
MACRO="${TOPDIR}/scripts/setup_rarexsec.C"

# Ensure library path and includes are visible to ROOT
case "$(uname -s)" in
  Darwin) export DYLD_LIBRARY_PATH="${LIBDIR}:${DYLD_LIBRARY_PATH}" ;;
  *)      export LD_LIBRARY_PATH="${LIBDIR}:${LD_LIBRARY_PATH}" ;;
esac
export ROOT_INCLUDE_PATH="${INCDIR}:${ROOT_INCLUDE_PATH}"

libext=$([ "$(uname -s)" = Darwin ] && echo dylib || echo so)
cmd=".L ${MACRO}; setup_rarexsec(\"${LIBDIR}/librarexsec.${libext}\",\"${INCDIR}\");"

# Start ROOT, run setup, then forward any user macro call e.g. '-- -q my.C'
root -l -q -e "$cmd" "$@"
