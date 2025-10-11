# rarexsec

## Build instructions

1. Clone the repository and enter it:
   ```bash
   git clone https://github.com/<your-org>/rarexsec.git
   cd rarexsec
   ```
2. Build the project libraries from the repository root:
   ```bash
   make -j"$(nproc)"
   ```
   The compiled artifacts are placed under `build/lib` and `build/bin`.  Add
   `-j"$(nproc)"` if you want to compile in parallel.
3. (Optional) Install the headers, libraries, and helper scripts into a prefix:
   ```bash
   make install PREFIX=/desired/install/location
   ```

## Running the example ROOT macro

After building, the rarexsec libraries must be discoverable by ROOT.  The `scripts/rarexsec-root.sh` wrapper sets up the include and library paths and executes any macro you pass to it.  From the repository root run:

```bash
./scripts/rarexsec-root.sh -b -q macros/inspect_simulation_samples.C
```

This command loads the generated libraries, configures ROOT include paths, and runs the `inspect_simulation_samples()` entry point defined in `macros/inspect_simulation_samples.C`.

Alternatively, you can call the setup macro manually from ROOT by providing the absolute library and include directories:

```bash
root -l -q -e 'gROOT->LoadMacro("scripts/setup_rarexsec.C"); setup_rarexsec("$PWD/build/lib/librarexsec.so","$PWD/include")' macros/example_macro.C
```

### Macro multiple entry points

Define the functions together and call the one you need in a single command:

```bash
# macros/analysis_tools.C
void scan_nominal() { /* ... */ }
void scan_systematics() { /* ... */ }

./scripts/rarexsec-root.sh -b -q 'macros/analysis_tools.C(\"scan_systematics()\")'
```
