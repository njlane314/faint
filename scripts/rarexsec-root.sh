#!/usr/bin/env bash

set -euo pipefail

TOPDIR="${RAREXSEC:-$(cd "$(dirname "$0")"/.. && pwd)}"
export RAREXSEC="$TOPDIR"

LIBDIR="$TOPDIR/build/lib"; [[ -d "$LIBDIR" ]] || LIBDIR="$TOPDIR/lib"
INCDIR="$TOPDIR/include"
SETUP="$TOPDIR/scripts/setup_rarexsec.C"; [[ -f "$SETUP" ]] || SETUP="$TOPDIR/setup_rarexsec.C"
LIB="$LIBDIR/librarexsec.so"
MACRO="$TOPDIR/analysis/main.C"

: "${RAREXSEC_CFG:=$TOPDIR/data/samples.json}"
: "${RAREXSEC_CONFIG:=$RAREXSEC_CFG}"
: "${RAREXSEC_TREE:=events}"
: "${RAREXSEC_BEAMLINE:=numi-fhc}"
: "${RAREXSEC_PERIODS:=run1}"

export RAREXSEC_CFG RAREXSEC_CONFIG RAREXSEC_TREE RAREXSEC_BEAMLINE RAREXSEC_PERIODS

export LD_LIBRARY_PATH="$LIBDIR:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$INCDIR:${ROOT_INCLUDE_PATH:-}"

root -l -b -q -e "gROOT->LoadMacro(\"$SETUP\"); setup_rarexsec(\"$LIB\",\"$INCDIR\"); gROOT->LoadMacro(\"$MACRO\"); main();"
