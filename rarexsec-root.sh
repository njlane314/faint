#!/usr/bin/env bash
set -euo pipefail
if [[ -z "${RAREXSEC:-}" ]]; then
  SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  if [[ -f "${SCRIPT_DIR}/setup_rarexsec.C" ]] || [[ -d "${SCRIPT_DIR}/build" ]]; then
    TOPDIR="${SCRIPT_DIR}"
  else
    TOPDIR="$( cd "${SCRIPT_DIR}/.." && pwd )"
  fi
else
  TOPDIR="${RAREXSEC}"
fi
export RAREXSEC="${TOPDIR}"
LIBDIR=""
for d in "${TOPDIR}/build/lib" "${TOPDIR}/lib"; do
  [[ -d "$d" ]] && { LIBDIR="$d"; break; }
done
: "${LIBDIR:=${TOPDIR}/build/lib}"
INCDIR=""
for d in "${TOPDIR}/include" "${TOPDIR}/include/rarexsec"; do
  [[ -d "$d" ]] && { INCDIR="$d"; break; }
done
: "${INCDIR:=${TOPDIR}/include}"
MACRO=""
for f in "${TOPDIR}/setup_rarexsec.C" "${TOPDIR}/scripts/setup_rarexsec.C"; do
  [[ -f "$f" ]] && { MACRO="$f"; break; }
done
if [[ -z "$MACRO" ]]; then
  echo "rarexsec-root: could not locate setup_rarexsec.C" >&2
  exit 1
fi
uname_s="$(uname -s || true)"
LIBEXT=$([[ "$uname_s" == "Darwin" ]] && echo "dylib" || echo "so")
LIBPATH="${LIBDIR}/librarexsec.${LIBEXT}"
if [[ ! -f "$LIBPATH" ]]; then
  ALT=$([[ "$LIBEXT" == "so" ]] && echo "dylib" || echo "so")
  [[ -f "${LIBDIR}/librarexsec.${ALT}" ]] && LIBPATH="${LIBDIR}/librarexsec.${ALT}"
fi
export LD_LIBRARY_PATH="${LIBDIR}:${LD_LIBRARY_PATH:-}"
if [[ "$uname_s" == "Darwin" ]]; then
  export DYLD_LIBRARY_PATH="${LIBDIR}:${DYLD_LIBRARY_PATH:-}"
fi
export ROOT_INCLUDE_PATH="${INCDIR}:${ROOT_INCLUDE_PATH:-}"
JSON_INC_PATH="${JSON_INC:-${NLOHMANN_JSON_INC:-}}"
[[ -n "${JSON_INC_PATH}" ]] && export ROOT_INCLUDE_PATH="${JSON_INC_PATH}:${ROOT_INCLUDE_PATH}"
CALL=""
while (( "$#" )); do
  case "$1" in
    -h|--help) echo "Usage: $(basename "$0") -c|--call MacroName"; exit 0 ;;
    -c|--call) CALL="${2:-}"; shift 2 ;;
    --) shift; break ;;
    *) echo "Usage: $(basename "$0") -c|--call MacroName" >&2; exit 2 ;;
  esac
done
if [[ -z "$CALL" ]]; then
  echo "Usage: $(basename "$0") -c|--call MacroName" >&2
  exit 2
fi
BOOT_CODE="gROOT->LoadMacro(\\\"${MACRO}\\\"); setup_rarexsec(\\\"${LIBPATH}\\\",\\\"${INCDIR}\\\");"
EXEC_CODE="rx_call(\\\"${CALL}\\\");"
root -l -b -q -e "${BOOT_CODE} ${EXEC_CODE}"