#!/usr/bin/env python3
"""Small helper that mirrors the historical dk2nu example script.

The original version failed with ``AttributeError: Dk2NuFlux`` when the
ROOT dictionary for dk2nu had not been loaded yet.  A number of grid
nodes ship only the bare ROOT python module so we need to explicitly
load the dk2nu shared library before accessing ``ROOT.Dk2NuFlux``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

try:
    import ROOT  # type: ignore
except Exception as exc:  # pragma: no cover - import error is fatal
    sys.stderr.write(
        "Failed to import ROOT. Ensure the ROOT environment is initialised.\n"
    )
    raise


_DK2NU_LIB_CANDIDATES: tuple[str, ...] = (
    "libdk2nuTree",
    "libdk2nuTree_v2",
    "libdk2nuTree_cpp",
    "libdk2nu",
)


def ensure_dk2nu_dictionary_loaded() -> str | None:
    """Attempt to load the dk2nu ROOT dictionary.

    Returns the library name that succeeded, or ``None`` if the
    dictionary was already available.  Raises ``RuntimeError`` with a
    detailed explanation if the class remains unavailable.
    """

    if hasattr(ROOT, "Dk2NuFlux"):
        return None

    tried: list[str] = []
    for lib in _DK2NU_LIB_CANDIDATES:
        status = ROOT.gSystem.Load(lib)
        tried.append(f"{lib} (status={status})")
        if status >= 0 and hasattr(ROOT, "Dk2NuFlux"):
            return lib

    raise RuntimeError(
        "Could not load the dk2nu dictionary. Tried: "
        + ", ".join(tried)
        + ".\n"
        "Make sure the dk2nu UPS product is setup (e.g. `setup dk2nu`) or"
        " the dk2nu shared libraries are on LD_LIBRARY_PATH."
    )


def write_pattern_file(output: Path, patterns: Iterable[str]) -> None:
    output.write_text("\n".join(str(p) for p in patterns) + "\n")


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Instantiate ROOT.Dk2NuFlux")
    parser.add_argument(
        "--fhc-pattern",
        default="/pnfs/uboone/persistent/users/bnayak/flux_files/"
        "uboone_beamsim_g4.10.4/me000z200i/run*/files/"
        "g4numi_minervame_me000z200i_*.root",
        help="Glob that resolves to forward horn current files.",
    )
    parser.add_argument(
        "--rhc-pattern",
        default="/pnfs/uboone/persistent/users/bnayak/flux_files/"
        "uboone_beamsim_g4.10.4/me000z-200i/run*/files/"
        "g4numi_minervame_me000z-200i_*.root",
        help="Glob that resolves to reverse horn current files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="dk2nu_g4104_fhc_rhc.root",
        help="Output ROOT filename.",
    )
    parser.add_argument(
        "--filelist",
        default="dk2nu_g4104_fhc_rhc.files",
        help="File containing the horn-polarity patterns.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    filelist_path = Path(args.filelist)
    write_pattern_file(filelist_path, (args.fhc_pattern, args.rhc_pattern))
    print(f"Wrote horn-polarity patterns to {filelist_path}")
    print(f"  FHC: {args.fhc_pattern}")
    print(f"  RHC: {args.rhc_pattern}")
    print(f"Output will be stored in {args.output}")

    try:
        loaded_lib = ensure_dk2nu_dictionary_loaded()
        if loaded_lib:
            print(f"Loaded dk2nu dictionary from {loaded_lib}")
        flux = ROOT.Dk2NuFlux(True, str(filelist_path), str(args.output))
    except RuntimeError as exc:
        sys.stderr.write(str(exc) + "\n")
        return 2
    except AttributeError:
        sys.stderr.write(
            "ROOT does not provide Dk2NuFlux even after loading candidate"
            " libraries. Check your dk2nu installation.\n"
        )
        return 2

    # The variable is unused but keeping a reference ensures PyROOT keeps the
    # underlying C++ instance alive for the remainder of the script.
    _ = flux
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
