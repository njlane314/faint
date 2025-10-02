#!/usr/bin/env bash
set -e

if [ -z "$RAREXSEC" ]; then
  TOPDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
else
  TOPDIR="$RAREXSEC"
fi

LIBDIR="${TOPDIR}/build/lib"
INCDIR="${TOPDIR}/include"
MACRO="${TOPDIR}/setup_rarexsec.C"

JSON_INC_PATH="${JSON_INC:-${NLOHMANN_JSON_INC:-}}"
if [ -n "${JSON_INC_PATH}" ]; then
  export ROOT_INCLUDE_PATH="${JSON_INC_PATH}:${ROOT_INCLUDE_PATH}"
fi

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

root -l -q -e "gROOT->LoadMacro(\"${MACRO}\"); setup_rarexsec(\"${LIBDIR}/librarexsec.${LIBEXT}\",\"${INCDIR}\");" "$@"