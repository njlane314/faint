#!/usr/bin/env bash
set -e

if [ -z "$RAREXSEC" ]; then
  SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  if [ -f "${SCRIPT_DIR}/setup_rarexsec.C" ] || [ -d "${SCRIPT_DIR}/build" ]; then
    TOPDIR="${SCRIPT_DIR}"
  else
    TOPDIR="$( cd "${SCRIPT_DIR}/.." && pwd )"
  fi
else
  TOPDIR="$RAREXSEC"
fi

lib_candidates=(
  "${TOPDIR}/build/lib"
  "${TOPDIR}/lib"
)
for candidate in "${lib_candidates[@]}"; do
  if [ -d "${candidate}" ]; then
    LIBDIR="${candidate}"
    break
  fi
done
: "${LIBDIR:=${TOPDIR}/build/lib}"

inc_candidates=(
  "${TOPDIR}/include"
  "${TOPDIR}/include/rarexsec"
)
for candidate in "${inc_candidates[@]}"; do
  if [ -d "${candidate}" ]; then
    INCDIR="${candidate}"
    break
  fi
done
: "${INCDIR:=${TOPDIR}/include}"

macro_candidates=(
  "${TOPDIR}/setup_rarexsec.C"
  "${TOPDIR}/scripts/setup_rarexsec.C"
)
for candidate in "${macro_candidates[@]}"; do
  if [ -f "${candidate}" ]; then
    MACRO="${candidate}"
    break
  fi
done

if [ -z "${MACRO:-}" ]; then
  echo "rarexsec-root: could not locate setup_rarexsec.C" >&2
  exit 1
fi

JSON_INC_PATH="${JSON_INC:-${NLOHMANN_JSON_INC:-}}"
if [ -n "${JSON_INC_PATH}" ]; then
  export ROOT_INCLUDE_PATH="${JSON_INC_PATH}:${ROOT_INCLUDE_PATH}"
fi

export LD_LIBRARY_PATH="${LIBDIR}:${LD_LIBRARY_PATH}"
export ROOT_INCLUDE_PATH="${INCDIR}:${ROOT_INCLUDE_PATH}"

LIBEXT="so"

root -l -q -e "gROOT->LoadMacro(\"${MACRO}\"); setup_rarexsec(\"${LIBDIR}/librarexsec.${LIBEXT}\",\"${INCDIR}\");" "$@"