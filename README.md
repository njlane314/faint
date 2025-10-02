# Framework for Analysis of Infrequent cross‑secTions

The rarexsec project provides reusable building blocks for neutrino cross-section
studies together with helper scripts and example ROOT macros.

## Prerequisites

Before building, ensure the following tools are available on your system:

- A C++17-capable compiler such as GCC or Clang.
- [ROOT 6](https://root.cern/) with the `root-config`, `rootcling` (or
  `rootcint`), and `root` executables on your `PATH`.
- The [nlohmann/json](https://github.com/nlohmann/json) single-header
  library.  Many ROOT builds, including those distributed through CVMFS,
  already provide it.  If it is installed in a non-standard location, set the
  `NLOHMANN_JSON_INC` environment variable (or make variable) to the directory
  that contains the `nlohmann/json.hpp` header so the build system can locate
  it.
- Standard development utilities: `make`, `bash`, and `git`.

You can quickly check that ROOT is discoverable with `root-config --version`.

## Build instructions

1. Clone the repository and enter it:
   ```bash
   git clone https://github.com/<your-org>/rarexsec.git
   cd rarexsec
   ```
2. Build the project libraries from the provided `build/Makefile`:
   ```bash
   cd build
   make -j"$(nproc)"
   ```
   The compiled artifacts are placed under `build/lib` and `build/bin`.  Add
   `-j"$(nproc)"` if you want to compile in parallel.
3. (Optional) Install the headers, libraries, and helper scripts into a prefix:
   ```bash
   make install PREFIX=/desired/install/location
   ```

If ROOT is not detected, the build will fall back to producing the core shared
and static libraries only; dictionary generation is skipped in that case.

## Running the example ROOT macro

After building, the rarexsec libraries must be discoverable by ROOT.  The
`rarexsec-root.sh` wrapper sets up the include and library paths and
executes any macro you pass to it.  From the repository root run:

```bash
./rarexsec-root.sh -b -q macros/example_macro.C
```

This command loads the generated libraries, configures ROOT include paths, and
runs the `example_macro()` entry point defined in `macros/example_macro.C`.  Drop
`-b -q` to open an interactive ROOT session preloaded with the rarexsec API.

Alternatively, you can call the setup macro manually from ROOT by providing the
absolute library and include directories:

```bash
root -l -q -e 'setup_rarexsec("$PWD/build/lib/librarexsec.so","$PWD/include")' macros/example_macro.C
```

Replace the library suffix with `.dylib` on macOS.  Once the environment is
configured, any custom macro can use the headers under `include/rarexsec/` and
link against the libraries in `build/lib/`.
