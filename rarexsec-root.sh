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
escape_cpp_string() {
  local str="${1//\\/\\\\}"
  printf '%s' "${str//\"/\\\"}"
}

RAREXSEC_CALL="${CALL}"
export RAREXSEC_CALL
ESC_MACRO="$(escape_cpp_string "${MACRO}")"
ESC_LIB="$(escape_cpp_string "${LIBPATH}")"
ESC_INC="$(escape_cpp_string "${INCDIR}")"
TMP_MACRO="$(mktemp "${TMPDIR:-/tmp}/rarexsec-root-XXXX.C")"
trap 'rm -f "${TMP_MACRO}"; unset RAREXSEC_CALL' EXIT
cat >"${TMP_MACRO}" <<EOF
void rarexsec_root_entry() {
  gROOT->LoadMacro("${ESC_MACRO}");
  setup_rarexsec("${ESC_LIB}","${ESC_INC}");
  const char* call = gSystem->Getenv("RAREXSEC_CALL");
  if (!call || !*call) {
    ::Error("rarexsec_root_entry", "RAREXSEC_CALL is not set");
    return;
  }
  rx_call(call);
}
rarexsec_root_entry();
EOF
root -l -b -q "${TMP_MACRO}"
